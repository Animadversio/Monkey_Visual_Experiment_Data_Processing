% Create Animation for Dual evolution in Reduced Dimension Evolution
Animal = "Beto";
MatStats_path = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(MatStats_path, Animal+"_HessBGEvolStats.mat"), 'HEStats')
%%
Animal = "Alfa";
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics"; 
load(fullfile(mat_dir, Animal + "_BigGAN_FC6_Evol_Stats.mat"), 'BFEStats'); 

%% 
result_dir = "E:\OneDrive - Washington University in St. Louis\BigGAN_Evol_Movies";
ExpType = "BigGAN_FC6";
EStats = BFEStats;

for Expi = 13:length(EStats)
fprintf("Processing BigGAN-FC6 Evolution Exp %d\n",Expi)
thread_n = EStats(Expi).evol.thread_num;
% ui = EStats(Expi).evol.unit_in_pref_chan(1);
% assert(all(ui==EStats(Expi).evol.unit_in_pref_chan))
% the following snippet is to get rid of 0 unit (null channel)
prefchan_id = find((EStats(Expi).units.spikeID == EStats(Expi).evol.pref_chan(1))); % Note unit U will be included here. 
unit_in_pref_chan = EStats(Expi).evol.unit_in_pref_chan; % this numbering corresponds to A,B,C... U is not included. 
% chid = find((EStats(Expi).units.unit_num_arr == unit_in_pref_chan(1)) & (EStats(Expi).units.spikeID == EStats(Expi).evol.pref_chan(1))); 
ui = unit_in_pref_chan(1); %find(prefchan_id==chid);
Window = 50:200;
% Find the image name and score for the best image in each block.
% Compute this before hand for faster video generation 
% Sort a image name list for each gen each thread. 
imgColl = repmat("", EStats(Expi).evol.block_n-1,thread_n);
sortImgColl = cell(EStats(Expi).evol.block_n-1,thread_n);
sortScoreColl = cell(EStats(Expi).evol.block_n-1,thread_n);
maxscoreColl = zeros(EStats(Expi).evol.block_n-1,thread_n);
minscoreColl = zeros(EStats(Expi).evol.block_n-1,thread_n);
meanscoreColl = zeros(EStats(Expi).evol.block_n-1,thread_n);
for thr_i = 1:thread_n
for blocki = 1:EStats(Expi).evol.block_n-1
    gen_scores = squeeze(mean(EStats(Expi).evol.psth{thr_i,blocki}(ui,Window,:),[1,2]));
    [maxScore, maxIdx] = max(gen_scores); % Idx in local generation
    maximgidx = EStats(Expi).evol.idx_seq{thr_i,blocki}(maxIdx); % get the trial index array for this block, get the trial index for block max
    imgfullfn = ls(fullfile(EStats(Expi).meta.stimuli, EStats(Expi).imageName(maximgidx)+"*"));
    assert(~isempty(imgfullfn), "Image not found %s",fullfile(EStats(Expi).meta.stimuli, EStats(Expi).imageName(maximgidx)+"*"))
%     imgColl(blocki,thr_i) = fullfile(EStats(Expi).meta.stimuli, imgfullfn); % this is a temporary solution. the suffix may not be .bmp % FIXED @ Jan.29
    maxscoreColl(blocki,thr_i) = max(gen_scores);
    minscoreColl(blocki,thr_i) = min(gen_scores);
    meanscoreColl(blocki,thr_i) = mean(gen_scores);
    suffix = imgfullfn(end-3:end);
    [sortScore, sortIdx] = sort(gen_scores,'Descend');
    sortimgidx = EStats(Expi).evol.idx_seq{thr_i,blocki}(sortIdx); % get the trial index array for this block, get the trial index for block max
    sortImgColl{blocki,thr_i} = arrayfun(@(idx)fullfile(EStats(Expi).meta.stimuli, ...
        string(EStats(Expi).imageName{idx})+suffix), sortimgidx); % full image path for sorted image name list for each block each thread
    sortScoreColl{blocki,thr_i} = sortScore;
end
end
%% Generation averaged psth and sem
% evol_stim_fr = cellfun(@(psth)mean(psth,3),EStats(Expi).evol.psth,'UniformOutput',false);
% evol_stim_fr = cell2mat(reshape(evol_stim_fr',1,1,[],thread_n));
% evol_stim_sem = cellfun(@(psth)std(psth,0,3)/sqrt(size(psth,3)),EStats(Expi).evol.psth,'UniformOutput',false);
% evol_stim_sem = cell2mat(reshape(evol_stim_sem',1,1,[],thread_n));
% % Scalor score for evolution
% score_avg = cellfun(@(psth)mean(psth(:,51:200,:),'all') - mean(psth(:,1:50,:),'all'),EStats(Expi).evol.psth);
% score_sem = cellfun(@(psth)std(squeeze(mean(psth(:,51:200,:),[1,2])))...
%     /sqrt(size(psth,3)),EStats(Expi).evol.psth); % - mean(psth(:,1:50,:),[1,2])
%% Generate Movies with Color coded frame
v = VideoWriter(fullfile(savepath,compose('%s_BGFC6Evol_Exp%02d_imbatch_frame',Animal,Expi)));
v.FrameRate = 2;
open(v);
% Set up figure canvas, get the object for use.
h4=figure(6);set(6,'position',[263          34        1611         849]);clf;
T = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
ax1 = nexttile(1); ax2 = nexttile(2);
stimparts = split(EStats(Expi).meta.stimuli,"\");
% expday = datetime(EStats.meta.expControlFN(1:6),'InputFormat','yyMMdd');
% fdrnm = compose("%s-Chan%02d", stimparts{end-1}, pref_chan(1));
ST = suptitle(compose("%s BigGAN FC6 Evol Exp %02d pref chan %s", ...
    stimparts{end-1}, Expi, EStats(Expi).units.unit_name_arr(EStats(Expi).units.pref_chan_id(1))));
% Plot the 2 batch of images
for blocki = 1:EStats(Expi).evol.block_n-1
    set(h4,"CurrentAxes",ax1)  % change the image shown
    sortScore = sortScoreColl{blocki,1};
    imglist1 = score_frame_image_arr(sortImgColl{blocki,1}, sortScoreColl{blocki,1}, [minscoreColl(blocki,1), maxscoreColl(blocki,1)], parula); 
    montage(imglist1)
    title(compose("Gen%d Rate Max %.1f Mean %.1f Min %.1f", blocki, maxscoreColl(blocki,1), ...
        meanscoreColl(blocki,1), minscoreColl(blocki,1) )) ; 
    if thread_n==2,set(h4,"CurrentAxes",ax2)
    imglist2 = score_frame_image_arr(sortImgColl{blocki,2}, sortScoreColl{blocki,2}, [minscoreColl(blocki,2), maxscoreColl(blocki,2)], parula); 
    montage(imglist2)
    title(compose("Gen%d Rate Max %.1f Mean %.1f Min %.1f", blocki, maxscoreColl(blocki,2), ...
        meanscoreColl(blocki,2), minscoreColl(blocki,2) )) ; 
    end
    drawnow;
%     pause(0.1)
    Fs = getframe(h4);
    writeVideo(v,Fs);
end
close(v);
%% Generate Movies
savepath = result_dir; mkdir(savepath) % fullfile(result_dir, compose("%s_Evol_"))
color_seq = EStats(Expi).color_seq;
v = VideoWriter(fullfile(savepath,compose('%s_BGFC6Evol_Exp%02d_imbatch',Animal,Expi)));
v.FrameRate = 2;
open(v);
% Set up figure canvas, get the object for use.
h4=figure(5);set(5,'position',[263          34        1611         849]);clf;
T = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
ax1 = nexttile(1); ax2 = nexttile(2);
stimparts = split(EStats(Expi).meta.stimuli,"\");
ST = suptitle(compose("%s BigGAN FC6 Evol Exp %02d pref chan %s", ...
    stimparts{end-1}, Expi, EStats(Expi).units.unit_name_arr(EStats(Expi).units.pref_chan_id(1))));
% Plot the 2 batch of images
for blocki = 1:EStats(Expi).evol.block_n-1
    set(h4,"CurrentAxes",ax1)  % change the image shown
    montage(sortImgColl{blocki,1}); 
    title(compose("Gen%d Rate Max %.1f Mean %.1f Min %.1f", blocki, maxscoreColl(blocki,1), ...
        meanscoreColl(blocki,1), minscoreColl(blocki,1) )) ; 
    if thread_n==2,set(h4,"CurrentAxes",ax2)
    montage(sortImgColl{blocki,2}); 
    title(compose("Gen%d Rate Max %.1f Mean %.1f Min %.1f", blocki, maxscoreColl(blocki,2), ...
        meanscoreColl(blocki,2), minscoreColl(blocki,2) )) ; 
    end
    drawnow;
%     pause(0.1)
    Fs = getframe(h4);
    writeVideo(v,Fs);
end
close(v);
end
