function score = scorePopulationResponse(actmat, targetMat, meanMat, stdMat, maskMat, score_mode)
actmat_norm = (actmat - meanMat) ./ stdMat;
targmat_norm = (targetMat - meanMat) ./ stdMat;
if strcmp(score_mode,'cosine')
score = nansum(actmat_norm(:,2:end).*targetMat(:,2:end).*maskMat(:,2:end),'all');
fprintf("Using dot product actually....")
elseif strcmp(score_mode,'dot')
% score = nansum(actmat_norm(:,2:end).*targetMat(:,2:end).*maskMat(:,2:end),'all');
score = nansum(actmat_norm(:,2:end).*targmat_norm(:,2:end).*maskMat(:,2:end),'all');
elseif strcmp(score_mode,'dot_IT')
maskMat = (maskMat&([1:64]'<33));
% score = nansum(actmat_norm(:,2:end).*targetMat(:,2:end).*maskMat(:,2:end),'all'); % First time! 
score = nansum(actmat_norm(:,2:end).*targmat_norm(:,2:end).*maskMat(:,2:end),'all');
elseif strcmp(score_mode,'dot_chan')
meanVec = nansum(meanMat(:,2:end),2);
stdVec = sqrt(nansum(stdMat(:,2:end).^2,2));
targetVec = nansum(targetMat(:,2:end),2);
actvec_norm = (nansum(actmat(:,2:end),2) - meanVec) ./ stdVec;
score = nansum(actvec_norm.*targetVec,'all');
elseif strcmp(score_mode,'corr_IT')
maskMat = (maskMat&([1:64]'<33));
score = corr(actmat_norm(maskMat),targmat_norm(maskMat), 'Rows', 'complete');
elseif strcmp(score_mode,'corr_V1')
maskMat = (maskMat&(([1:64]'>32)&[1:64]'<49));
score = corr(actmat_norm(maskMat),targmat_norm(maskMat), 'Rows', 'complete');
elseif strcmp(score_mode,'corr_V4')
maskMat = (maskMat&([1:64]'>48));
score = corr(actmat_norm(maskMat),targmat_norm(maskMat), 'Rows', 'complete');
elseif strcmp(score_mode,'corr_ITV4')
maskMat = (maskMat&(([1:64]'<33)|[1:64]'>48));
score = corr(actmat_norm(maskMat),targmat_norm(maskMat), 'Rows', 'complete');
elseif strcmp(score_mode,'corr_V1V4')
maskMat = (maskMat&([1:64]'>32));
score = corr(actmat_norm(maskMat),targmat_norm(maskMat), 'Rows', 'complete');
elseif strcmp(score_mode,'corr')
% score = corr(actmat_norm(maskMat),targetMat(maskMat), 'Rows', 'complete'); % older version
score = corr(actmat_norm(maskMat),targmat_norm(maskMat), 'Rows', 'complete');
elseif strcmp(score_mode,'L1_V1V4')
maskMat = (maskMat&([1:64]'>32));
score = -nanmean(abs(actmat_norm(maskMat)-targmat_norm(maskMat)), 'all'); % negative to keep maximizing
elseif strcmp(score_mode,'L1_IT')
maskMat = (maskMat&([1:64]'<33));
score = -nanmean(abs(actmat_norm(maskMat)-targmat_norm(maskMat)), 'all'); % negative to keep maximizing
elseif strcmp(score_mode,'L1')
score = -nanmean(abs(actmat_norm(maskMat)-targmat_norm(maskMat)), 'all'); % negative to keep maximizing
elseif strcmp(score_mode,'MSE_V1V4')
maskMat = (maskMat&([1:64]'>32));
score = -nanmean((actmat_norm(maskMat)-targmat_norm(maskMat)).^2, 'all'); % negative to keep maximizing
elseif strcmp(score_mode,'MSE_IT')
maskMat = (maskMat&([1:64]'<33));
score = -nanmean((actmat_norm(maskMat)-targmat_norm(maskMat)).^2, 'all'); % negative to keep maximizing
elseif strcmp(score_mode,'MSE_ITV4')
maskMat = (maskMat&(([1:64]'<33)|([1:64]'>48)));
score = -nanmean((actmat_norm(maskMat)-targmat_norm(maskMat)).^2, 'all'); % negative to keep maximizing
elseif strcmp(score_mode,'MSE')
score = -nanmean((actmat_norm(maskMat)-targmat_norm(maskMat)).^2, 'all'); % negative to keep maximizing
else
error("Method not implemented")
end
end