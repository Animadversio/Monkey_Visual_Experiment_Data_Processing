classdef OptimizeChannelWeightAverageLayer < nnet.layer.RegressionLayer
    % OptimizeChannelWeightAverageLayer   Implementation
    % Modified from OptimizeChannelAverageLayer
    %   Copyright 2019 Binxu Wang binxu.wang@wustl.edu
    
    properties
        % LearnableParameters   Learnable parameters for the layer
        %   This layer has no learnable parameters.
        LearnableParameters = nnet.internal.cnn.layer.learnable.PredictionLearnableParameter.empty();
    end
    
    properties (SetAccess = private)
        % HasSizeDetermined   Specifies if all size parameters are set
        %   For pooling layers, this is always true.
        HasSizeDetermined = true
        
        % ChannelNumbers   An array of channel numbers to optimize for
        Channels
        % 
        WeightsChannel
        % ChannelWeights   An array of Weights to put on the channels. 
        Weights
        % Mathematically equivalent to optimize the cos of angle between the
        % feature vector and Weights vector. 
        nWeights % how many group of weights
    end
    
    methods
        function this = OptimizeChannelWeightAverageLayer(channels, varargin)
            this.Name = 'optimizechannelweightaverage';
            this.WeightsChannel = channels;
            % Set channels
            if nargin == 2
                this.Weights = varargin{1};
                if size(this.Weights,1)==length(channels)
                    this.nWeights = size(this.Weights,2);
                    this.Weights = this.Weights';
                elseif size(this.Weights,2)==length(channels)
                    this.nWeights = size(this.Weights,1);
                else
                    error("If Weights provided, length of weight should match that of channels");
                end
            else
                this.Weights = ones(numel(channels),1);
                this.nWeights = 1;
            end
            this.Channels = 1:this.nWeights;
        end
        
        
        function G = forwardLoss(this, X, ~)
            % Returns a dummy value 1, as a gpuArray if X is a gpuArray.
            % A value is necessary for interaction with forward and
            % backpropagation of the internal network; this value is not
            % actually used.
            
            %G = ones(1, 'like', X);
            G = zeros([1 this.nWeights], 'like', X);
            for w = 1:this.nWeights
            for i = 1:numel(this.WeightsChannel)
                G(w) = G(w) + mean(X(:,:,this.WeightsChannel(i),w) * this.Weights(i),'all'); % dX(:,:,this.Channels(i),i) = 1;
            end
            end
            disp(G)
        end
        
        function dX = backwardLoss(this, X, labels)
            dX = zeros(size(X), 'like', X);
            numChannels = size(this.Channels, 2);
            for w = 1:this.nWeights
            for i = 1:numel(this.WeightsChannel)
                dX(:,:,this.WeightsChannel(i),w) = this.Weights(i); % dX(:,:,this.Channels(i),i) = 1;
            end
            end
        end
        
        function tf = isValidInputSize(~, ~) %#ok<STOUT>
            error(message('nnet_cnn:internal:cnn:visualize:layer:OptimizeChannelAverageLayer:IsValidInputSizeProhibited'));
        end
    end
end