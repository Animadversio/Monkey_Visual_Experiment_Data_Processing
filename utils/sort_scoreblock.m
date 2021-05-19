function [score_m,score_s,blockvec] = sort_scoreblock(blockarr,scorearr)
% sort an array of scores according to the block array labels. compute the
% mean and std for each block. 
% really useful function to summarize multiple evolution trajectories into
% a mean one. 
blockvec = min(blockarr):max(blockarr);
score_m = [];score_s = [];
for blocki = min(blockarr):max(blockarr)
    score_m(blocki) = mean(scorearr(blockarr==blocki));
    score_s(blocki) = sem(scorearr(blockarr==blocki));
end
end