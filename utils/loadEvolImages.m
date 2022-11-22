function [im_gen, imnames] = loadEvolImages(stimfolder, thread, block, suffix)
% Util function to load in evolved image of certain block x thread from a
% folder
% block: negative `block` variable means counting from the last block. e.g. -1 means
%        the last block,-2 means the penultimate block
% thread: Python numbering convention here, 0 for first thread, 1 for 2nd
%        thread
if nargin<=3, suffix=".bmp";end
if block <= 0 % negative block number means counting from the last one
allimgs = string(ls(stimfolder+"\block*thread*gen*"));
lastgeni = regexp(allimgs(end), "block(\d*)_thread\d*_gen", 'tokens');
lastgeni = str2num(lastgeni{1});
blockid = lastgeni + block + 1;
else
blockid = block;
end
imnames = string(ls(fullfile(stimfolder,compose("block%03d_thread%03d_gen*",blockid,thread))));
im_gen = cellfun(@(impath)imread(fullfile(stimfolder,impath)),imnames,'Uni',0);
end