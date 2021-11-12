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
frames = cropTrimVideo(vid, fullfile(vidroot,"bird2.avi"), 201.5, [116.5100  220  500  500]);
frames = cropTrimVideo(vid, fullfile(vidroot,"bird3.avi"), 502.5, [400  220  500  500]);
frames = cropTrimVideo(vid, fullfile(vidroot,"bird4.avi"), 600.5, [400  220  500  500]);
frames = cropTrimVideo(vid, fullfile(vidroot,"bird5.avi"), 803.5, [350  220  500  500]);
%%
vidroot = "N:\Stimuli\natural_animal_movies"; 
frames = cropTrimVideo(vid, fullfile(vidroot,"bird6.avi"), 2000, [280,250,450,450],true);
%%
frames = cropTrimVideo(vid, fullfile(vidroot,"bird7.avi"), 1300, [550  100  600  600],true);
frames = cropTrimVideo(vid, fullfile(vidroot,"bird8.avi"), 1501.5, [600  220  500  500],true);
frames = cropTrimVideo(vid, fullfile(vidroot,"bird9.avi"), 1600, [400  250  450  450],true);
frames = cropTrimVideo(vid, fullfile(vidroot,"bird10.avi"), 1905, [1  220  500  500],true);
frames = cropTrimVideo(vid, fullfile(vidroot,"bird11.avi"), 2539, [1  120  600  600],true);
frames = cropTrimVideo(vid, fullfile(vidroot,"bird12.avi"), 2250, [350  220  500  500],true);
%%
%%
movpath = "C:\Users\Ponce lab\Documents\ml2a-monk\freeviewMovie\theyre babies godsends i - tiny vs the hierarchy.avi";
newviddir = "C:\Users\Ponce lab\Documents\ml2a-monk\freeviewMovie";
vid = VideoReader(movpath);
%%
frames = TrimVideo(vid, fullfile(newviddir,"monkeybaby_1.avi"), 30, 150, true);
frames = TrimVideo(vid, fullfile(newviddir,"monkeybaby_2.avi"), 180, 300, true);
frames = TrimVideo(vid, fullfile(newviddir,"monkeybaby_3.avi"), 330, 450, true);
frames = TrimVideo(vid, fullfile(newviddir,"monkeybaby_4.avi"), 480, 600, true);
frames = TrimVideo(vid, fullfile(newviddir,"monkeybaby_5.avi"), 630, 750, true);
%%
frames = TrimVideo(vid, fullfile(newviddir,"monkeybaby_6.avi"), 550, 565, true);
frames = TrimVideo(vid, fullfile(newviddir,"monkeybaby_7.avi"), 680, 695, true);
frames = TrimVideo(vid, fullfile(newviddir,"monkeybaby_8.avi"), 120, 135, true);
%%
frames = TrimVideo(vid, fullfile(newviddir,"monkeybaby_9.avi"), 150, 165, true);
frames = TrimVideo(vid, fullfile(newviddir,"monkeybaby_10.avi"), 360, 375, true);
frames = TrimVideo(vid, fullfile(newviddir,"monkeybaby_11.avi"), 500, 515, true);
frames = TrimVideo(vid, fullfile(newviddir,"monkeybaby_12.avi"), 530, 545, true);
%%
function frames = cropTrimVideo(vid, newvidpath, startTime, cropWdw, doFrames)
if nargin == 4, doFrames = true; end
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
winopen(newvidpath)
if doFrames
pathparts = split(newvidpath,"\");
nameparts = split(pathparts{end},".");
mkdir(fullfile(pathparts{1:end-1},nameparts{1}))
for i = 1:numel(frames)
    impath = fullfile(pathparts{1:end-1},nameparts{1},compose("%s_%03d.jpg",nameparts{1},i));
    imwrite(frames{i},impath)
end
end
end

function frames = TrimVideo(vid, newvidpath, startTime, endTime, doFrames)
if nargin == 4, doFrames = true; end
vid.CurrentTime = startTime;
framesN = round((endTime - startTime) * vid.FrameRate);
frames = {};
for i = 1:framesN
    img = vid.readFrame();
    frames{end+1} = img;%imcrop(img, cropWdw);% (cropWdw(2)+1:cropWdw(2)+cropWdw(4),cropWdw(1)+1:cropWdw(1)+cropWdw(3),:)
end

vidnew = VideoWriter(newvidpath);
vidnew.FrameRate = vid.FrameRate;vidnew.open()
for i = 1:numel(frames)
    vidnew.writeVideo(frames{i});
end
close(vidnew);
winopen(newvidpath)
if doFrames
pathparts = split(newvidpath,"\");
nameparts = split(pathparts{end},".");
mkdir(fullfile(pathparts{1:end-1},nameparts{1}))
for i = 1:numel(frames)
    impath = fullfile(pathparts{1:end-1},nameparts{1},compose("%s_%03d.jpg",nameparts{1},i));
    imwrite(frames{i},impath)
end
end
end
