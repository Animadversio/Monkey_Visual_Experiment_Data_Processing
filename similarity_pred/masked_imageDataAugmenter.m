classdef masked_imageDataAugmenter < handle
    
    
    
    methods
        
        function B = augment(self,A)
            %augment Augment input image data.
            %
            %   augmentedImage = augment(augmenter,A) performs image
            %   augmentation on the input A. A can be a numeric image or a
            %   cell array of numeric and categorical images. When A is a cell array of
            %   images, the augment function performs identical
            %   augmentations on all the images and returns a cell array B
            %   of augmented images. Images in a cell array can be of
            %   different sizes and types.
            if iscell(A)
                B = cell(size(A));
                self.AffineTransforms = zeros(3,3,length(A)); % 3 x 3 x batchSize matrix
                
                % Select Rand values before the for loop so that same
                % augmentations are applied to all the images
                self.selectRandValues();
                
                for img = 1:numel(A)
                    if (isnumeric(A{img}) || iscategorical(A{img}) || islogical(A{img})) && (ndims(A{img}) < 4)
                        [B{img},tform] = augmentSingleImage(self,A{img});
                        self.AffineTransforms(:,:,img) = tform;
                    else
                        error(message('nnet_cnn:imageDataAugmenter:invalidImageCellarray',img));
                    end
                end
                
            elseif (isnumeric(A) || iscategorical(A) || islogical(A)) && (ndims(A) < 4)
                self.selectRandValues();
                [B,tform] = augmentSingleImage(self,A);
                self.AffineTransforms = tform;
                
            else
                error(message('nnet_cnn:imageDataAugmenter:invalidImage'));
            end
            
        end
    end
    
    methods(Hidden)
        
        function [B,tform] = augmentSingleImage(self,A)
            
            if (isnumeric(A) || iscategorical(A) || islogical(A)) && (ndims(A) < 4)
                
                interp = 'linear';
                imgIsCategorical = false;
                imgIsLogical = false;
                fillValue = manageFillValue(A,self.FillValue);
                
                if(iscategorical(A))
                    % Convert image to numeric for warp and save
                    % current categories for converting the result
                    % back to categorical
                    imgIsCategorical = true;
                    interp = 'nearest';
                    if iscategorical(A)
                        cats = categories(A);
                        fillValue = str2double(cats{end}) + 1;
                    end
                    A = double(A);
                end
                
                if islogical(A)
                    interp = 'nearest';
                    A = uint8(A);
                    imgIsLogical = true;
                end
                
                % tform could be different for different images as
                % centered rotation is performed and that depends
                % on size of the image
                tform = self.makeAffineTransform(size(A));
                
                B = nnet.internal.cnnhost.warpImage2D(A,tform(1:3,1:2),interp,fillValue);
                
                if imgIsCategorical
                    % If original image was categorical, covert the
                    % result to categorical
                    B = categorical(B, 1:numel(cats), cats);
                end
                
                if imgIsLogical
                    B = logical(B);
                end
            else
                error(message('nnet_cnn:imageDataAugmenter:invalidImageCellarray',img));
            end
            
        end
    end
end