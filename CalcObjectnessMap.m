objpath = "E:\DL_Projects\Vision\objectness-release-v2.2";
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
Animal = "Beto";
load(fullfile(mat_dir,Animal+"_ImageRepr.mat"), 'ReprStats');
load(fullfile(mat_dir,Animal+"_Evol_Stats.mat"), 'EStats');
%%
for Expi = 1:45
    imgsize = EStats(Expi).evol.imgsize;
    imsize_pix = imgsize * 40;
    imgpos = EStats(Expi).evol.imgpos;
    prefchan = EStats(Expi).evol.pref_chan;
    img = ReprStats(Expi).Evol.BestBlockAvgImg;%Manif.BestImg;
    boxes = runObjectness(img,100);
    [objHeatMap, rawmap] = computeObjectnessHeatMap(img,boxes(boxes(:,5)>0.7,:));
    figure(2);
    subplot(131);imshow(img),drawBoxes(boxes(boxes(:,5)>0.7,:));
    subplot(132);imshow(img)
    subplot(133);imshow(objHeatMap)
    title(compose("Exp %d prefchan %d %.1f deg [%.1f %.1f]\nObjectiveness: Top 8 range [%.3f %.3f]",...
	Expi,prefchan,imgsize,imgpos(1),imgpos(2),max(boxes(:,5)),min(boxes(:,5))))
    pause;
end