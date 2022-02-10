function scores = scorePopulationResponse_vec(acttsr, targetMat, meanMat, stdMat, maskMat, score_mode,array_layout)
% acttsr: shape [64, MAXUNUM+1, imgN]
% 	 We assume the acttsr 's  activation have been baseline subtracted. 
% targetMat: shape [64, MAXUNUM+1, ]
%    We assume the targetMat 's  activation have been baseline subtracted. 
% meanMat: shape [64, MAXUNUM+1, ]
% stdMat: shape [64, MAXUNUM+1, ]
% maskMat: boolean shape [64, MAXUNUM+1, ] True for valid units
% score_mode: str for scoring mode, it contains a criterion like "corr", "MSE", "dot", "L1"
%    And it can contain areal qualification like "_V1", "_V4", "_V1V4"
%    Score will be computed in the intersected map between maskMat and score_mode. 
% array_layout: it would be "Alfa" or "Beto_new" correspond to different mappings of channel number to brain areas.
% 
% Return: 
%  scores: vector of shape [imgN]
if nargin == 6, array_layout = "Alfa"; end

acttsr_norm = (acttsr - meanMat) ./ stdMat;
targmat_norm = (targetMat - meanMat) ./ stdMat;
chanmsk = parse_mode2mask(score_mode,array_layout); % parse out the areal qualification in score_mode
finalmsk = maskMat & chanmsk; % do intersection. 
imgN = size(acttsr_norm,3);
% Apply the mask and get vector format population activation. 
vec_reprs_norm = reshape(acttsr_norm(repmat(finalmsk,1,1,imgN)),[],imgN); % shape [selected unitN, imgN]
vec_targs_norm = targmat_norm(finalmsk); % shape [selected unitN, 1]
if(sum(finalmsk,'all')==0),scores=nan(1,imgN); return; end
if contains(score_mode,'dot')
scores = (vec_targs_norm)'*(vec_reprs_norm);
elseif contains(score_mode,'corr')
scores = corr(vec_reprs_norm, vec_targs_norm);
elseif contains(score_mode,'MSE')
scores = -nanmean((vec_targs_norm - vec_reprs_norm).^2,1);
elseif contains(score_mode,'L1')
scores = -nanmean(abs(vec_targs_norm - vec_reprs_norm),1);
elseif contains(score_mode,'cosine')
scores = 1-pdist2(vec_targs_norm', vec_reprs_norm', 'cosine');
% fprintf("Using dot product actually....")
else
error("Method not implemented")
end
% 
% if strcmp(score_mode,'cosine')
% score = nansum(actmat_norm(:,2:end).*targetMat(:,2:end).*maskMat(:,2:end),'all');
% fprintf("Using dot product actually....")
% elseif strcmp(score_mode,'dot')
% % score = nansum(actmat_norm(:,2:end).*targetMat(:,2:end).*maskMat(:,2:end),'all');
% score = nansum(actmat_norm(:,2:end).*targmat_norm(:,2:end).*maskMat(:,2:end),'all');
% elseif strcmp(score_mode,'dot_IT')
% maskMat = (maskMat&([1:64]'<33));
% % score = nansum(actmat_norm(:,2:end).*targetMat(:,2:end).*maskMat(:,2:end),'all'); % First time! 
% score = nansum(actmat_norm(:,2:end).*targmat_norm(:,2:end).*maskMat(:,2:end),'all');
% elseif strcmp(score_mode,'dot_chan')
% meanVec = nansum(meanMat(:,2:end),2);
% stdVec = sqrt(nansum(stdMat(:,2:end).^2,2));
% targetVec = nansum(targetMat(:,2:end),2);
% actvec_norm = (nansum(actmat(:,2:end),2) - meanVec) ./ stdVec;
% score = nansum(actvec_norm.*targetVec,'all');
% elseif strcmp(score_mode,'corr_IT')
% maskMat = (maskMat&([1:64]'<33));
% score = corr(actmat_norm(maskMat),targmat_norm(maskMat), 'Rows', 'complete');
% elseif strcmp(score_mode,'corr_V1')
% maskMat = (maskMat&(([1:64]'>32)&[1:64]'<49));
% score = corr(actmat_norm(maskMat),targmat_norm(maskMat), 'Rows', 'complete');
% elseif strcmp(score_mode,'corr_V4')
% maskMat = (maskMat&([1:64]'>48));
% score = corr(actmat_norm(maskMat),targmat_norm(maskMat), 'Rows', 'complete');
% elseif strcmp(score_mode,'corr_ITV4')
% maskMat = (maskMat&(([1:64]'<33)|[1:64]'>48));
% score = corr(actmat_norm(maskMat),targmat_norm(maskMat), 'Rows', 'complete');
% elseif strcmp(score_mode,'corr_V1V4')
% maskMat = (maskMat&([1:64]'>32));
% score = corr(actmat_norm(maskMat),targmat_norm(maskMat), 'Rows', 'complete');
% elseif strcmp(score_mode,'corr')
% % score = corr(actmat_norm(maskMat),targetMat(maskMat), 'Rows', 'complete'); % older version
% score = corr(actmat_norm(maskMat),targmat_norm(maskMat), 'Rows', 'complete');
% elseif strcmp(score_mode,'L1_V1V4')
% maskMat = (maskMat&([1:64]'>32));
% score = -nanmean(abs(actmat_norm(maskMat)-targmat_norm(maskMat)), 'all'); % negative to keep maximizing
% elseif strcmp(score_mode,'L1_IT')
% maskMat = (maskMat&([1:64]'<33));
% score = -nanmean(abs(actmat_norm(maskMat)-targmat_norm(maskMat)), 'all'); % negative to keep maximizing
% elseif strcmp(score_mode,'L1')
% score = -nanmean(abs(actmat_norm(maskMat)-targmat_norm(maskMat)), 'all'); % negative to keep maximizing
% elseif strcmp(score_mode,'MSE_V1V4')
% maskMat = (maskMat&([1:64]'>32));
% score = -nanmean((actmat_norm(maskMat)-targmat_norm(maskMat)).^2, 'all'); % negative to keep maximizing
% elseif strcmp(score_mode,'MSE_IT')
% maskMat = (maskMat&([1:64]'<33));
% score = -nanmean((actmat_norm(maskMat)-targmat_norm(maskMat)).^2, 'all'); % negative to keep maximizing
% elseif strcmp(score_mode,'MSE_ITV4')
% maskMat = (maskMat&(([1:64]'<33)|([1:64]'>48)));
% score = -nanmean((actmat_norm(maskMat)-targmat_norm(maskMat)).^2, 'all'); % negative to keep maximizing
% elseif strcmp(score_mode,'MSE')
% score = -nanmean((actmat_norm(maskMat)-targmat_norm(maskMat)).^2, 'all'); % negative to keep maximizing
% else
% error("Method not implemented")
% end
end