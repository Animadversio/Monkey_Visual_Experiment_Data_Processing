%% Code to process the Style experiment. (Selectivity)
%%
Animal = "Both"; Set_Path;
expftr = contains(ExpRecord.Exp_collection, "Style");
% ExpRecord.Expi>=44 & ...%ExpRecord.Expi<=40 & 
%     contains(ExpRecord.expControlFN,"selectivity") & ...
[meta_new,rasters_new,~,Trials_new] = Project_Manifold_Beto_loadRaw(find(expftr),Animal);
%%
result_dir = "C:\Users\ponce\OneDrive - Washington University in St. Louis\Style_Tune";
Expi = 1;
%%
Triali = 4;
Trials = Trials_new{Triali};
rasters = rasters_new{Triali};
meta = meta_new{Triali};
if contains(meta.ephysFN,"Beto"), Animal = "Beto"; elseif contains(meta.ephysFN,"Alfa"), Animal = "Alfa"; else keyboard; end
savedir = fullfile(result_dir,compose("%s_PilotExp%d",Animal,Expi));
mkdir(savedir)
imgnm_uniq = string(unique(Trials.imageName));
unit_name_arr = generate_unit_labels(meta.spikeID);
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);
%
epocn_arr = arrayfun(@(nm)str2num(string(regexp(nm, 'epoch(\d*)_', 'tokens'))), imgnm_uniq);
GANmask = contains(imgnm_uniq,"real_A");
pixmask = contains(imgnm_uniq,"real_B");
photmsk = contains(imgnm_uniq,"fake_B");
imgposs = cell2mat(Trials.XY);
uniq_pos = unique(imgposs,'rows');
posmask = {};
for iPos = 1:size(uniq_pos)
    posmask{iPos} = (imgposs(:,1) == uniq_pos(iPos,1) & imgposs(:,2) == uniq_pos(iPos,2));
end
epoc_list = epocn_arr(GANmask);
% %% not taking position into account
% idx_col = cellfun(@(nm) find(contains(Trials.imageName, nm)), imgnm_uniq, 'UniformOutput', false);
% nam_mat = reshape(imgnm_uniq,3,[])';
% nam_mat = nam_mat(:,[2,1,3]);
% idx_mat = reshape(idx_col,3,[])';
% idx_mat = idx_mat(:,[2,1,3]); % real_A  fake_B  real_B
% Take the image position into account 
idx_col = cellfun(@(msk) cellfun(@(nm) find(contains(Trials.imageName, nm) & msk), ...
                imgnm_uniq, 'UniformOutput', false),...
                posmask, 'UniformOutput', false); 
idx_col = cat(2,idx_col{:});
idx_mat = reshape(idx_col,3,[],size(uniq_pos,1));
idx_mat = permute(idx_mat([2,1,3],:,:),[2,1,3]);
nam_mat = cellfun(@(idx)unique(Trials.imageName(idx)), idx_mat);
% Take position into account. Select the best location to do comparison
score_mat_all = cell2mat(cellfun(@(idx) reshape(mean(rasters(:,Window,idx),[2,3]),1,1,1,[]),idx_mat,'UniformOutput',false));
% 89  3  2 88, [image id, style id, position id, channel] 
zscore_mat_all = zscore(score_mat_all,0,[1,2,3]); % z score the response within each channel. 
zscore_mat_bestpos = squeeze(max(zscore_mat_all,[],3)); % Use the mean score to the best location. discard the other. 
%% Pool over position
StylStats = repmat(struct(),size(rasters,1),1);
Window = 51:200;
for iCh = 1:size(rasters,1)
score_mat = cellfun(@(idx) mean(rasters(iCh,Window,idx),[2,3]),idx_mat);
[~,P12,~,STATS12] = ttest(score_mat(:,1),score_mat(:,2));
[~,P13,~,STATS13] = ttest(score_mat(:,1),score_mat(:,3));
[~,P23,~,STATS23] = ttest(score_mat(:,2),score_mat(:,3));
statstr = sprintf("GAN-pix: p=%.1E(t=%.1f)\nGAN-photo: p=%.1E(t=%.1f)\npix-photo: p=%.1E(t=%.1f)",...
    P12,STATS12.tstat,P13,STATS13.tstat,P23,STATS23.tstat);
StylStats(iCh).p_GAN_pix = P12; StylStats(iCh).t_GAN_pix = STATS12.tstat;
StylStats(iCh).p_GAN_phot = P13; StylStats(iCh).t_GAN_phot = STATS13.tstat;
StylStats(iCh).p_pix_phot = P23; StylStats(iCh).t_pix_phot = STATS23.tstat;
StylStats(iCh).chan = meta.spikeID(iCh);
StylStats(iCh).unit = unit_num_arr(iCh);
StylStats(iCh).uname = unit_name_arr(iCh);

figure(3);
plot(score_mat')
xticks([1,2,3])
xticklabels(["GAN", "pix2pix", "photo"])%,'Fontsize',16
% xtickangle(350)
title(sprintf("%s Style Exp Chan %d\n%s",Animal,meta.spikeID(iCh),statstr))
colororder(jet(size(idx_mat,1))) % Blue to Red is the direction of 1st to last gen
figure(1);
plot(score_mat)
legend(["GAN", "pix2pix", "photo"])%,'Fontsize',16
title(sprintf("%s Style Exp Chan %d\n%s",Animal,meta.spikeID(iCh),statstr))
drawnow
pause(0.5)
saveas(1,fullfile(savedir, compose("%s_Exp%d_epoc_chan%s.png",Animal,Expi,unit_name_arr(iCh))))
saveas(3,fullfile(savedir, compose("%s_Exp%d_styl_chan%s.png",Animal,Expi,unit_name_arr(iCh))))
% % break
end
save(fullfile(savedir, compose("Stat.mat")),'StylStats')
StylStatTab = struct2table(StylStats);
writetable(StylStatTab, fullfile(savedir, compose("Stat.csv")))
%%
StylStatTab = readtable("C:\Users\ponce\OneDrive - Washington University in St. Louis\Style_Tune\Beto_PilotExp1\Stat.csv");
%%
t_Col = [StylStatTab.t_GAN_phot, StylStatTab.t_GAN_pix, StylStatTab.t_pix_phot]; 
figure(5);
imagesc(t_Col)
xticks([1,2,3])
xticklabels(["GAN-photo", "GAN-pix", "pix-photo"])%,'Fontsize',16
yticks(1:size(rasters,1))
yticklabels(unit_name_arr)
ytickangle(340)
%%
for iCh = 1:size(rasters,1)
score_mat = zscore_mat_bestpos(:,:,iCh);
[~,P12,~,STATS12] = ttest(score_mat(:,1),score_mat(:,2));
[~,P13,~,STATS13] = ttest(score_mat(:,1),score_mat(:,3));
[~,P23,~,STATS23] = ttest(score_mat(:,2),score_mat(:,3));
statstr = sprintf("GAN-pix: p=%.1E(t=%.1f)\nGAN-photo: p=%.1E(t=%.1f)\npix-photo: p=%.1E(t=%.1f)",...
        P12,STATS12.tstat,P13,STATS13.tstat,P23,STATS23.tstat);
StylStats(iCh).p_GAN_pix_bestpos = P12; StylStats(iCh).t_GAN_pix_bestpos = STATS12.tstat;
StylStats(iCh).p_GAN_phot_bestpos = P13; StylStats(iCh).t_GAN_phot_bestpos = STATS13.tstat;
StylStats(iCh).p_pix_phot_bestpos = P23; StylStats(iCh).t_pix_phot_bestpos = STATS23.tstat;
figure(6);
plot(score_mat')
xticks([1,2,3])
xticklabels(["GAN", "pix2pix", "photo"])%,'Fontsize',16
% xtickangle(350)
title(sprintf("%s Style Exp Chan %d\n%s",Animal,meta.spikeID(iCh),statstr))
colororder(jet(size(idx_mat,1))) % Blue to Red is the direction of 1st to last gen
figure(7);
plot(score_mat)
legend(["GAN", "pix2pix", "photo"])%,'Fontsize',16
title(sprintf("%s Style Exp Chan %d\n%s",Animal,meta.spikeID(iCh),statstr))
drawnow
pause(0.7)
% break
saveas(7,fullfile(savedir, compose("%s_Exp%d_epoc_chan%s_bestpos.png",Animal,Expi,unit_name_arr(iCh))))
saveas(6,fullfile(savedir, compose("%s_Exp%d_styl_chan%s_bestpos.png",Animal,Expi,unit_name_arr(iCh))))
end
%% there seems a trend in V1 V4 Array but not IT array. 
save(fullfile(savedir, compose("Stat.mat")),'StylStats')
StylStatTab = struct2table(StylStats);
%%
writetable(StylStatTab, fullfile(savedir, compose("Stat.csv")))
%%
t_Col = [StylStatTab.t_GAN_phot_bestpos, StylStatTab.t_GAN_pix_bestpos, StylStatTab.t_pix_phot_bestpos]; 
figure(8);
imagesc(t_Col)
xticks([1,2,3])
xticklabels(["GAN-photo", "GAN-pix", "pix-photo"])%,'Fontsize',16
yticks(1:size(rasters,1))
yticklabels(unit_name_arr)
ytickangle(340)
%%
figure(4);  
plot(squeeze(mean(zscore_mat_bestpos,1)),'Linewidth',1.2)
colorseq = zeros(size(zscore_mat_bestpos,3),3);
colorseq(meta.spikeID<=32,3) = 1; % IT color is blue
colorseq(meta.spikeID>=49,[2,3]) = 0.5; % V4 color is purple
colorseq(meta.spikeID>=33 & meta.spikeID<=48,1) = 1; % V1 color is red
colororder(colorseq)
xticks([1,2,3])
xticklabels(["GAN-photo", "GAN-pix", "pix-photo"])%,'Fontsize',16
%%
ITmsk = meta.spikeID<=32;
V4msk = meta.spikeID>=49;
V1msk = meta.spikeID>=33 & meta.spikeID<=48;
figure(10); set(10,"position",[680         176        1033         802]);
subplot(131)
plot(squeeze(mean(zscore_mat_bestpos(:,:,ITmsk),1)),'Linewidth',1.2,'Color',[0,0,1])
xticks([1,2,3]);xticklabels(["GAN", "pix2pix", "photo"])
ylabel("zscore firing rate [50,200] ms")
title("IT",'FontSize',14)
subplot(132)
plot(squeeze(mean(zscore_mat_bestpos(:,:,V4msk),1)),'Linewidth',1.2,'Color',[1,0,1])
xticks([1,2,3]);xticklabels(["GAN", "pix2pix", "photo"])
title("V4",'FontSize',14)
subplot(133)
plot(squeeze(mean(zscore_mat_bestpos(:,:,V1msk),1)),'Linewidth',1.2,'Color',[1,0,0])
xticks([1,2,3]);xticklabels(["GAN", "pix2pix", "photo"])
title("V1",'FontSize',14)
suptitle(compose("%s Pilot Style Exp",Animal))
saveas(10,fullfile(savedir, compose("%s_Exp%d_area_cmp.png",Animal,Expi)));
%%
figure(8);
imagesc(squeeze(mean(zscore_mat_bestpos,1))') % mean zscore rsp for the 3 styles 
colorbar()
xticks([1,2,3]);xticklabels(["GAN", "pix2pix", "photo"])%,'Fontsize',16
yticks(1:size(rasters,1));yticklabels(unit_name_arr);ytickangle(340)
title(sprintf("%s Style Exp\nmean zscored activity to 3 styles in 3 areas",Animal))
line([0 3.5],[0.5, 0.5]+sum(meta.spikeID<33),'color','r')
line([0 3.5],[0.5, 0.5]+sum(meta.spikeID<49),'color','r')
saveas(8,fullfile(savedir, compose("%s_Exp%d_area_cmp_pcolor.png",Animal,Expi)));

%%
figure(11);
imagesc(squeeze(rasters(40,:,:))')
caxis(prctile(rasters(59,:,:),[2,99.5],'all')')
ylabel("image Num");xlabel("Time")