function scores = scorePopulationVec_direction(reprMat, targetVec, meanVec, stdVec, maskVec, score_mode)
% Similar to `scorePopulationResponse_vec` but not the same structure. 
% Signature
%    scores = scorePopulationVec_direction(reprMat, targetVec, ...
%                    meanVec, stdVec, maskVec, score_mode)
% 
% reprMat: shape [ChanNum, imgN]
% 	 We assume the acttsr 's  activation have been baseline subtracted. 
% targetVec: shape [ChanNum, ]
%    We assume the targetVec 's a direction vector in the normalized activation space.
%    i.e. population activation with baseline subtracted and centralized using mean and std. 
% meanVec: shape [ChanNum, ] 
% stdVec: shape [ChanNum, ]
% maskVec: boolean shape [ChanNum, ] True for valid units
%
% score_mode: str for scoring mode, it contains a criterion like "corr", "MSE", "dot", "L1"
%    And it can contain areal qualification like "_V1", "_V4", "_V1V4"
%    Score will be computed in the intersected map between maskMat and score_mode. 
% array_layout: it would be "Alfa" or "Beto_new" correspond to different mappings of channel number to brain areas.
% 
% Return: 
%  scores: vector of shape [imgN]
if nargin <= 6, array_layout = "Alfa"; end

reprMat_norm = (reprMat - meanVec) ./ stdVec;
% targmat_norm = (targetMat - meanMat) ./ stdMat;
% chanmsk = parse_mode2maskVec(score_mode,array_layout,spikeID); % parse out the areal qualification in score_mode
finalmsk = maskVec; % do intersection. 
imgN = size(reprMat,2);
% Apply the mask and get vector format population activation. 
% vec_reprs_norm = reshape(acttsr_norm(repmat(finalmsk,1,1,imgN)),[],imgN); % shape [selected unitN, imgN]
vec_reprs_norm = reprMat_norm(finalmsk,:);
vec_targs_norm = reshape(targetVec(finalmsk),[],1); % shape [selected unitN, 1]
if(sum(finalmsk,'all')==0),scores=nan(1,imgN); return; end
if contains(score_mode,'dot')
scores = (vec_targs_norm)'*(vec_reprs_norm);
elseif contains(score_mode,'corr')
scores = corr(vec_reprs_norm, vec_targs_norm);
elseif contains(score_mode,'cosine')
scores = 1-pdist2(vec_targs_norm', vec_reprs_norm', 'cosine')';

elseif contains(score_mode,'MSE')
% scores = -nanmean((vec_targs_norm - vec_reprs_norm).^2,1);
error("MSE is not a suitable distance for direction")
elseif contains(score_mode,'L1')
% scores = -nanmean(abs(vec_targs_norm - vec_reprs_norm),1);
error("L1 is not a suitable distance for direction")
else
error("Method not implemented")
end
end