classdef Kraft_Ambience_Extractor < audioPlugin
    
    % myBasicSourcePlugin is a template for a basic source plugin. Use this
    % template to create your own basic source plugin.
    
    properties
        % Use this section to initialize properties that the end-user
        % interacts with.
        
        Width = 0.5;
        IOBufferSize;
        windowSize = 2048;
        overlapRatio = 0.75;
        hopSize;
        
        inputChanNum = 2;
        outputChanNum = 11;
        
        inputVectorSize = -1;
        inputVectorNumPerHop = -1;
        loadedVectorNumThisHop = -1;
        inputBuffer;
        outputBuffer;
        
        
        windowFuncInput;
        windowFuncOutput;
        
        Decorrelation_Strength = 1;
        Output_Gain = 1;
        
        RLR;
        RS_Fr;
        sigmaLR;
        sigmaBackFront;
        sigmaLoHi;
        H_A_L;
        H_A_R;
        H_A_Rr;
        H_A_Fr;
        H_A_Lo;
        H_A_Hi;
        
    end
    properties (Access = private)
        % Use this section to initialize properties that the end-user does
        % not interact with directly.
    end
    properties (Constant)
        PluginInterface = ...
            audioPluginInterface(...
            audioPluginParameter('Decorrelation_Strength', ...
            'Mapping', {'lin',0, 1}), ...
            audioPluginParameter('Width', ...
            'Mapping', {'lin', 0, 1}), ...
            audioPluginParameter('Output_Gain', ...
            'Mapping', {'lin', 0, 2}), ...
            'InputChannels',2, ...
            'OutputChannels',11, ...
            'PluginName','Kraft Ambience Extractor (2 in / 11 out)');
    end
    methods
        %------------------------------------------------------------------
        % Construct
        %------------------------------------------------------------------
        function plugin = Kraft_Ambience_Extractor
            plugin.hopSize = floor(plugin.windowSize*(1-plugin.overlapRatio));
            %initializeProcessingBuffer(plugin);
            plugin.windowFuncInput = repmat(hann(2048),1,plugin.inputChanNum);
            plugin.windowFuncOutput = repmat(hann(2048),1,plugin.outputChanNum);
            
            plugin.RLR = 2*rand(plugin.windowSize,1)-1;
            plugin.RS_Fr = 2*rand(plugin.windowSize,1)-1;
            sigmas = coder.load('kraftFilterDataSigma_improved.dat');
            plugin.sigmaLR = sigmas(:,1);
            plugin.sigmaBackFront = sigmas(:,2);
            plugin.sigmaLoHi = sigmas(:,3);
            
            %initializeFilters(plugin);
            reset(plugin);
        end
        
        function out = process(plugin,in)
            % This section contains instructions to produce the output
            % audio signal. Use plugin.MyProperty to access a property of
            % your plugin. Use getSamplesPerFrame(plugin) to get the frame
            % size used by the environment.
            
            if(plugin.inputVectorSize < 0)
                plugin.inputVectorSize = size(in,1);
                plugin.inputVectorNumPerHop = plugin.hopSize/plugin.inputVectorSize;
                plugin.loadedVectorNumThisHop = 0;
            end
            
            
            plugin.inputBuffer = [plugin.inputBuffer((plugin.inputVectorSize+1):end,:); in];
            plugin.loadedVectorNumThisHop = plugin.loadedVectorNumThisHop + 1;
            
            if (plugin.loadedVectorNumThisHop == plugin.inputVectorNumPerHop)
                %   process the signal here
                processedSig = Kraft_ambience_extraction(plugin, plugin.windowFuncInput.*plugin.inputBuffer);
                %   overlap and add
                plugin.outputBuffer = plugin.outputBuffer + plugin.windowFuncOutput.*processedSig;
            end
            
            %   output buffering
            out = plugin.outputBuffer(1:plugin.inputVectorSize,:);
            %   the earliest hop size of signal is discarded from buffer
            plugin.outputBuffer = [plugin.outputBuffer(plugin.inputVectorSize+1:end,:); zeros(plugin.inputVectorSize,plugin.outputChanNum)];

            if (plugin.loadedVectorNumThisHop == plugin.inputVectorNumPerHop)
                %   clear the loading count
                plugin.loadedVectorNumThisHop = 0;
            end
            
            
            
            % mid-side processing
            out(:,1:2) = midSideProcessing(plugin, out(:,1:2));
            out(:,3:4) = midSideProcessing(plugin, out(:,3:4));
            out(:,5:6) = midSideProcessing(plugin, out(:,5:6));
            out(:,7:8) = midSideProcessing(plugin, out(:,7:8));
        end
        function reset(plugin)
            % This section contains instructions to reset the plugin
            % between uses, or when the environment sample rate changes.
            plugin.inputVectorSize = -1;
            initializeProcessingBuffer(plugin);
            initializeFilters(plugin);
        end

        function set.Decorrelation_Strength(plugin, val)
            % This section contains instructions to execute if the
            % specified property is modified. Properties associated with
            % parameters are updated automatically. Use the set method to
            % execute more complicated instructions.
            plugin.Decorrelation_Strength = val;
            initializeFilters(plugin);
        end
        
        function out = midSideProcessing(plugin, in)
            mid = (in(:,1) + in(:,2)) / 2;
            side = (in(:,1) - in(:,2)) / 2;
            mid = mid * (1-plugin.Width)* 2;
            side = side * plugin.Width * 2;
            out = [mid + side, mid-side];
        end
        
        function out = Kraft_ambience_extraction(plugin, in)
            IN = fft(in, 2*plugin.windowSize);
            IN = IN(1:plugin.windowSize, :);
            IN_POWER = abs(IN).^2;
            gl = sqrt(IN_POWER(:,1)./sum(IN_POWER,2));
            gr = sqrt(IN_POWER(:,2)./sum(IN_POWER,2));
            AMBINECE_SIG = (gl.*IN(:,2) - gr.*IN(:,1))./(gl.*plugin.H_A_R - gr.*plugin.H_A_L);
            DIRECT_SIG = (plugin.H_A_R.*XL - plugin.H_A_L.*XR)./(gl.*plugin.H_A_R - gr.*plugin.H_A_L);
            
            OUT = [AMBINECE_SIG.*plugin.H_A_L,...
                   AMBINECE_SIG.*plugin.H_A_R]; 
            OUT = [OUT, OUT, OUT, OUT];
            
            OUT(:, 1:2) = [OUT(:,1).*plugin.H_A_Fr, OUT(:,2).*plugin.H_A_Fr];
            OUT(:, 5:6) = OUT(:, 1:2);
            OUT(:, 3:4) = [OUT(:,1).*plugin.H_A_Rr, OUT(:,2).*plugin.H_A_Rr];
            OUT(:, 7:8) = OUT(:, 3:4);
            
%             OUT(:, 1:2) = [OUT(:,1).*plugin.H_A_Lo, OUT(:,2).*plugin.H_A_Lo];
            OUT(:, 3:4) = [OUT(:,3).*plugin.H_A_Lo, OUT(:,4).*plugin.H_A_Lo];
%             OUT(:, 5:6) = [OUT(:,5).*plugin.H_A_Hi, OUT(:,6).*plugin.H_A_Hi];
            OUT(:, 7:8) = [OUT(:,7).*plugin.H_A_Hi, OUT(:,8).*plugin.H_A_Hi];
            
            OUT(:, 9:10) = [DIRECT_SIG.*gl,DIRECT_SIG.*gr];
            
            out = real(ifft(OUT, 2*plugin.windowSize));
            out = out(1:plugin.windowSize, :);
            %   setup output level
            %   multiple by two because we lose half power when doing fft
            %   and ifft with discarding the mirrored part
            out = out*plugin.Output_Gain*2;
        end
        
        function initializeProcessingBuffer(plugin)
            plugin.inputBuffer = zeros(plugin.windowSize, plugin.inputChanNum);
            plugin.outputBuffer = zeros(plugin.windowSize, plugin.outputChanNum);
        end
        
        function initializeFilters(plugin)
            decorrelation_strength_filter = plugin.sigmaLR*plugin.Decorrelation_Strength;

            RL =(1/pi)*atan(decorrelation_strength_filter.*plugin.RLR) + (1/2);
            RR = 1.- RL;

            plugin.H_A_L = RL.*complex(ones(size(RL)), zeros(size(RL)));
            plugin.H_A_R = RR.*complex(zeros(size(RR)), ones(size(RR)));
            

            RLr =(1/pi)*atan(decorrelation_strength_filter.*plugin.RS_Fr) + (1/2);
            RRr = 1.- RLr;
            
            
            plugin.H_A_Rr = RRr.*plugin.sigmaBackFront;
            plugin.H_A_Fr = 1-plugin.H_A_Rr;
            
            plugin.H_A_Rr = 2*plugin.H_A_Rr;
            plugin.H_A_Fr = 2*plugin.H_A_Fr;
            
            plugin.H_A_Lo = plugin.sigmaLoHi;
            plugin.H_A_Hi = 1-plugin.H_A_Lo;
        end
    end
end