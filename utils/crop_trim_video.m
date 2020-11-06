% movpath = "N:\Stimuli\natural_animal_movies\Der _Hamster_-Film720.mp4";
movpath = "N:\Stimuli\natural_animal_movies\Movies for Cats to Watch Birds - The Ultimate Movie & Video for Your Cat.mp4";
vid = VideoReader(movpath);

%%
cropWdw = [340,150,450,450];
%%
vid.CurrentTime = 246.7;
sampframe = vid.readFrame();
figure;
[J,rect2] = imcrop(sampframe);
%%
cropWdw = [rect2(1:2), mean(rect2(3:4)), mean(rect2(3:4))];%[340,150,450,450];
%%
vid.CurrentTime = 140;
frames = {};
for i = 1:120
    img = vid.readFrame();
    frames{end+1} = imcrop(img, cropWdw); % (cropWdw(2)+1:cropWdw(2)+cropWdw(4),cropWdw(1)+1:cropWdw(1)+cropWdw(3),:)
end
%%
figure;
montage(frames)
%%
vidroot = "N:\Stimuli\natural_animal_movies";
vidnew = VideoWriter(fullfile(vidroot,"bird1.avi"));
vidnew.FrameRate = 30;vidnew.open()
for i = 1:numel(frames)
    vidnew.writeVideo(frames{i});
%     img = vid.readFrame();
%     frames{end+1} = imcrop(img, cropWdw); % (cropWdw(2)+1:cropWdw(2)+cropWdw(4),cropWdw(1)+1:cropWdw(1)+cropWdw(3),:)
end
close(vidnew);

%%

vid.CurrentTime = 201.5;cropWdw = [116.5100  248.5100  500  500];
frames = {};
for i = 1:120
    img = vid.readFrame();
    frames{end+1} = imcrop(img, cropWdw);% (cropWdw(2)+1:cropWdw(2)+cropWdw(4),cropWdw(1)+1:cropWdw(1)+cropWdw(3),:)
end
vidroot = "N:\Stimuli\natural_animal_movies";
vidnew = VideoWriter(fullfile(vidroot,"bird2.avi"));
vidnew.FrameRate = 30;vidnew.open()
for i = 1:numel(frames)
    vidnew.writeVideo(frames{i});
end
close(vidnew);
winopen(fullfile(vidroot,"bird2.avi"))
%%
vidroot = "N:\Stimuli\natural_animal_movies"; 
frames = cropTrimVideo(vid, fullfile(vidroot,"bird1.avi"), 140, [340,150,450,450]);
%%
frames = cropTrimVideo(vid, fullfile(vidroot,"bird2.avi"), 201.5, [116.5100  248.5100  500  500]);
frames = cropTrimVideo(vid, fullfile(vidroot,"bird3.avi"), 502.5, [400  220  500  500]);
frames = cropTrimVideo(vid, fullfile(vidroot,"bird4.avi"), 600.5, [400  220  500  500]);
frames = cropTrimVideo(vid, fullfile(vidroot,"bird5.avi"), 803.5, [350  220  500  500]);
%%
function frames = cropTrimVideo(vid, newvidpath, startTime, cropWdw)
vid.CurrentTime = startTime;cropWdw = cropWdw;
frames = {};
for i = 1:120
    img = vid.readFrame();
    frames{end+1} = imcrop(img, cropWdw);% (cropWdw(2)+1:cropWdw(2)+cropWdw(4),cropWdw(1)+1:cropWdw(1)+cropWdw(3),:)
end

vidnew = VideoWriter(newvidpath);
vidnew.FrameRate = 30;vidnew.open()
for i = 1:numel(frames)
    vidnew.writeVideo(frames{i});
end
close(vidnew);
% winopen(newvidpath)
pathparts = split(newvidpath,"\");
nameparts = split(pathparts{end},".");
mkdir(fullfile(pathparts{1:end-1},nameparts{1}))
for i = 1:numel(frames)
    impath = fullfile(pathparts{1:end-1},nameparts{1},compose("%s_%03d.jpg",nameparts{1},i));
    imwrite(frames{i},impath)
end
end
