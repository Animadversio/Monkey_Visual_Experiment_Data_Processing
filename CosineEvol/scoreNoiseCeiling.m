function [ceil_mean, ceil_sem] = scoreNoiseCeiling(D,idx,score_method)
if nargin==2
score_method=idx;
for idx = 1:numel(D.pics_unique)
[ceil_mean_s, ceil_sem_s] = scoreNoiseCeiling(D,idx,score_method);
ceil_mean(idx) = ceil_mean_s;
ceil_sem(idx) = ceil_sem_s;
end
else
fprintf("Image %s\n",D.pics_unique{idx})
MAXUNUM = 4;
scores = [];
for i = 1:1000
trialrsp = D.responseTensor(:,:,idx)+D.responseStdTensor(:,:,idx).* randn(64,MAXUNUM+1);
score = scorePopulationResponse(trialrsp,D.responseTensor(:,:,idx),D.mean_mat,D.std_mat,D.F_P_mat<1E-6,score_method);
scores = [scores, score];
end
fprintf("%s Mean %.3f Std %.3f Sem %.3f (i.i.d. Gauss Noise)\n",score_method,mean(scores),std(scores),sem(scores))
ceil_mean = mean(scores);
ceil_sem = sem(scores);
end
end
% scorePopulationResponse(D.responseTensor(:,:,10),D.responseTensor(:,:,10),D.mean_mat,D.std_mat,D.F_P_mat<0.001,"corr_V1V4")

%%
% score_method = "corr_V1V4";
% for idx = 1:numel(D.pics_unique)
%     fprintf("%s\n",D.pics_unique{idx})
%     for score_method = ["corr_V1V4","MSE_V1V4"]
%     MAXUNUM = 4;
%     scores = [];
%     for i = 1:1000
%     trialrsp = D.responseTensor(:,:,idx)+D.responseStdTensor(:,:,idx).* randn(64,MAXUNUM+1);
%     score = scorePopulationResponse(trialrsp,D.responseTensor(:,:,idx),D.mean_mat,D.std_mat,D.F_P_mat<0.001,score_method);
%     scores = [scores, score];
%     end
%     fprintf("%s Mean %.3f Std %.3f Sem %.3f (i.i.d. Gauss Noise)\n",score_method,mean(scores),std(scores),sem(scores))
%     end
% end