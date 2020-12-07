%% Find out what feature in the PSTH of Movie and Static image response are the most similar
%  And what can assist decoding the most. 
%  Can we decode frame identity from instantaneous firing rate? 



%% Data computation PSTH correlation
MovPSTHWdw = [-20:100];
ImgPSTHWdw = [80:200];
kerL = 21;
wdw = meta_mov.rasterWindow; 
matchONidx = int32(matchTON) - wdw(1) + 1; % Timing for each frame

PSTH_corrmat = [];
for chid = 1:numel(meta_cmb.spikeID) % loop through static image units
iCh = meta_cmb.spikeID(chid);
iU = meta_cmb.unitID(chid);
static_PSTH = cellfun(@(psth) psth(chid, ImgPSTHWdw), psth_mean(:,5:9), 'Uni', 0);
chid2 = find(meta_mov.spikeID==iCh & meta_mov.unitID==iU); % find(meta_mov.spikeID==9)
movie_frPSTH = cell(numel(psthmov_mean), size(matchONidx, 1));
for iMv = 1:numel(psthmov_mean)
    for iTic = 1:size(matchONidx, 1)
        iTons = double(matchONidx(iTic,:));
        for iTon = iTons
           movie_frPSTH{iMv, iTic}(end+1,:) = psthmov_mean{iMv}(chid2, iTon+MovPSTHWdw);
        end
    end
end
static_smth_PSTH = cellfun(@(P)downsample(smoothdata(P,2,'sgolay',kerL)',int32(kerL/2)),static_PSTH,'uni',0);
movie_smth_PSTH = cellfun(@(P)downsample(smoothdata(P,2,'sgolay',kerL)',int32(kerL/2)),movie_frPSTH,'uni',0);
PSTH_corrmat1 = cellfun(@(psthM,psthS)corr(psthM(:,2), psthS),movie_smth_PSTH,static_smth_PSTH,'uni',1);
PSTH_corrmat2 = cellfun(@(psthM,psthS)corr(psthM(:,1), psthS),movie_smth_PSTH,static_smth_PSTH,'uni',1);
PSTH_corrmat = cellfun(@(psthM,psthS)corr(mean(psthM,2), psthS),movie_smth_PSTH,static_smth_PSTH,'uni',1);
fprintf("Chan%02d%s Median corr1 %.3f corr2 %.3f corr mean %.3f\n",iCh,char(64+iU),...
                        median(PSTH_corrmat2(:,2:5),'all'),...
                        median(PSTH_corrmat1(:,2:5),'all'),median(PSTH_corrmat(:,2:5),'all'))
% pause
end

%% Data computation PSTH correlation
wdw = meta_mov.rasterWindow; 
matchONidx = int32(matchTON) - wdw(1) + 1; % Timing for each frame

coher_synop = [];

% MovPSTHWdw = [-20:100];
ImgPSTHWdw = [80:200];
kerL = 30; DSr = 30; %(kerL+1)/2
WdwL = 120;
MovDelays = -120:20:120;
kerWids = [5, 11, 21, 31, 51, 81, 121];
for i = 1:numel(MovDelays)
for j = 1:numel(kerWids)
MovDelay = MovDelays(i);
kerL = kerWids(j);DSr = (kerL+1)/2;
MovPSTHWdw = [MovDelay:MovDelay+WdwL];

cc_col = [];
for chid = 1:numel(meta_cmb.spikeID) % loop through static image units
iCh = meta_cmb.spikeID(chid);
iU = meta_cmb.unitID(chid);
static_PSTH = cellfun(@(psth) psth(chid, ImgPSTHWdw), psth_mean(:,5:9), 'Uni', 0);
chid2 = find(meta_mov.spikeID==iCh & meta_mov.unitID==iU); % find(meta_mov.spikeID==9)
movie_frPSTH = cell(numel(psthmov_mean), size(matchONidx, 1));
for iMv = 1:numel(psthmov_mean)
    for iTic = 1:size(matchONidx, 1)
        iTons = double(matchONidx(iTic,:));
        for iTon = iTons
           movie_frPSTH{iMv, iTic}(end+1,:) = psthmov_mean{iMv}(chid2, iTon+MovPSTHWdw);
        end
    end
end
% static_smth_PSTH = cellfun(@(P)downsample(smoothdata(P,2,'sgolay',kerL)',DSr),static_PSTH,'uni',0);
% movie_smth_PSTH = cellfun(@(P)downsample(smoothdata(P,2,'sgolay',kerL)',DSr),movie_frPSTH,'uni',0);
static_smth_PSTH = cellfun(@(P)downsample(smoothdata(P,2,'movmean',kerL)',DSr),static_PSTH,'uni',0);
movie_smth_PSTH = cellfun(@(P)downsample(smoothdata(P,2,'movmean',kerL)',DSr),movie_frPSTH,'uni',0);
staticCatFeatVec = cat(1, static_smth_PSTH{:,2:5});
movieCatFeatVec = cat(1, movie_smth_PSTH{:,2:5});
[cc,pp] = corr(staticCatFeatVec,movieCatFeatVec); % first entry is the 2nd occurence, 2nd entry is the first
% fprintf("Chan%02d%s all Feat Correlation %.3f (%.1e) %.3f (%.1e)\n",iCh,char(64+iU),...
%                         cc(2),pp(2),cc(1),pp(1))
cc_col(chid, :) = cc;
% PSTH_corrmat1 = cellfun(@(psthM,psthS)corr(psthM(2,:)', psthS'),movie_smth_PSTH,static_smth_PSTH,'uni',1);
% PSTH_corrmat2 = cellfun(@(psthM,psthS)corr(psthM(1,:)', psthS'),movie_smth_PSTH,static_smth_PSTH,'uni',1);
% PSTH_corrmat = cellfun(@(psthM,psthS)corr(mean(psthM',2), psthS'),movie_smth_PSTH,static_smth_PSTH,'uni',1);
% fprintf("Chan%02d%s Median corr1 %.3f corr2 %.3f corr mean %.3f\n",iCh,iU,...
%                         median(PSTH_corrmat2(:,2:5),'all'),...
%                         median(PSTH_corrmat1(:,2:5),'all'),median(PSTH_corrmat(:,2:5),'all'))
end
coher_synop(i,j,1) = median(cc_col(:,1),'all');
coher_synop(i,j,2) = median(cc_col(:,2),'all');
end
end
%% 
save(fullfile(figdir,'FeatWdw_CoherMat.mat'),'coher_synop','kerWids','MovDelays','WdwL')
%%
figure(6);T=tiledlayout(2,1,'Padd','compact');
title(T,"Median Coherence of PSTH as a Function of PSTH Feature Window")
nexttile(1)
imagesc(coher_synop(:,:,2)')
yticks(1:numel(kerWids));yticklabels(kerWids);
ylabel("MovMean Kernel Width (ms)")
xticks(1:numel(MovDelays));xticklabels(MovDelays);
xlabel("120 Window Delay from frame onset (ms)")
colorbar();
title("Occurence in First Half")
nexttile(2)
imagesc(coher_synop(:,:,1)')
yticks(1:numel(kerWids));yticklabels(kerWids);
ylabel("MovMean Kernel Width (ms)")
xticks(1:numel(MovDelays));xticklabels(MovDelays);
xlabel("120 Window Delay from frame onset (ms)")
title("Occurence in Second Half")
colorbar();
saveas(6,fullfile(figdir, "Coher_WindowFeat_Map_Overall.png"))
savefig(6,fullfile(figdir, "Coher_WindowFeat_Map_Overall.fig"))
%%
