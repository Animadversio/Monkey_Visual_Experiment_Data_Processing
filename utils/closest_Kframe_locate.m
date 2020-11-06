function [matchfrid_col, matchTON_col, matchTOFF_col] = closest_Kframe_locate(mvnms, imnms_col, mvpath, impath, halfSep, K, frameMax)
% Adapt from Frame Locate 
% Input format: 
%  mvnms: cell array or string array of movie names. Char is not OK.  Without .avi. 
%  imnms_col: cell array the same length as mvnms
%  mvpath: Stimuli folder to load in movies
%  impath: Stimuli folder to load in images
%  K: closest K frames to search for, default 2
%  frameMax: First `frameMax` frames will be loaded from the movie and
%     compared. 
%
% Return :
%  matchfrid_col: a cell array of matched frame index for each img in that movie. 
%  matchTON_col/matchTOFF_col: same structure as `matchfrid_col` but it's the timing from movie onset in double. 
if nargin <= 6, frameMax = 48; end % max num of frame played is 48.
if nargin <= 5, K = 2; end % by default find 2 most similar frames
if nargin <= 4, halfSep = true; end
assert(length(imnms_col)==length(mvnms), ...
	"Img name collection should have the same number of entry as movie names.")

D = torchImDist();
matchfrid_col = {};
matchTON_col = {};
matchTOFF_col = {};
for iMv = 1:numel(mvnms)
% Movie frames
fprintf("Loading Movie Frames... ")
vid = VideoReader(fullfile(mvpath, mvnms{iMv}+".avi"));
frameDur = 1000 / vid.FrameRate; 
frames = {};
for i = 1:frameMax 
frames{end+1} = vid.readFrame();
end
frames_tsr = cell2mat(reshape(frames,1,1,1,[]));

% Key static images 
fprintf("Loading Static Images... ")
keyimgs = {};
for iImg = 1:numel(imnms_col{iMv})
    nm = imnms_col{iMv}(iImg);
    keyimgs{end+1} = imread(fullfile(impath, nm+".jpg"));
end
keyimgs_tsr = cell2mat(reshape(keyimgs,1,1,1,[]));

% image distance matrix
fprintf("Compute Static Frame Distance Matrix... ")
% imfrdistmat = []; % Static N - by - Frame N matrix 
% for i = 1:numel(keyimgs)
% imfrdistmat(i,:) = D.distance(keyimgs_tsr(:,:,:,i), frames_tsr);
% end
imfrdistmat = D.distmat2(keyimgs_tsr, frames_tsr, 30);
% Closest frame id
if halfSep 
sepidx = 25;
[mindist1, matchfrid1] = mink(imfrdistmat(:,1:sepidx),1,2); 
[mindist2, matchfrid2] = mink(imfrdistmat(:,sepidx+1:end),1,2); 
matchfrid = [matchfrid1, matchfrid2 + sepidx];
mindist = [mindist1, mindist2]; % No need to sort 
else
[mindist, matchfrid] = mink(imfrdistmat,K,2); % 2 closest frames in the movie. Along the 2nd dimension. 
fprintf("Sort frame idx...")
matchfrid = sort(matchfrid,2); % Sort the min idx in ascending order 
end
% Translate frame id into onset and offset timing of that frame
matchTON = (matchfrid - 1) * frameDur;
matchTOFF = (matchfrid) * frameDur;
matchfrid_col{end+1} = matchfrid;
matchTON_col{end+1} = matchTON;
matchTOFF_col{end+1} = matchTOFF;
fprintf("Done\n")
end
if numel(mvnms)==1
matchfrid_col = matchfrid_col{1};
matchTON_col = matchTON_col{1};
matchTOFF_col = matchTOFF_col{1};
end
fprintf("Done\n")
% usually we assume that's the same for each movie
clear D
end