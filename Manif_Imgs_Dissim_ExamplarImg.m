Animal = "Both"; Set_Path;
mat_dir = "O:\Mat_Statistics";

Animal = "Beto";
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'), 'EStats') 
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'), 'Stats') 
%%
pasu_idx_vec = reshape(Stats(12).ref.pasu_idx_grid',[],1); % reshape into one row. But transpose to make similar object adjacent
pasu_val_msk = ~cellfun(@isempty,pasu_idx_vec);
%% Find Examplar Images 
flag.doEvoRef = true;
P.topN = 3;
P.botN = 3;
for Expi = 1:numel(Stats)
bsl_VEC_ALL = [];
% Manifold images. 
si=1; ui=1;
score_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).manif.psth{si},[],1));
score_std_vec = cellfun(@(psth)std(mean(psth(ui,51:200,:),[1,2]),1,'all'),reshape(Stats(Expi).manif.psth{si},[],1));
score_bsl_VEC = cell2mat(cellfun(@(psth)squeeze(mean(psth(ui,1:45,:),[2])),reshape(Stats(Expi).manif.psth{si},[],1),'uni',0)); % single trial
bsl_VEC_ALL = [bsl_VEC_ALL; score_bsl_VEC];
[sortScore,sortId]=sort(score_vec,'Descend');
[maxScore,maxId]=max(score_vec);
manifImgnm = cellfun(@(idx)string(Stats(Expi).imageName(idx(1))),Stats(Expi).manif.idx_grid{1});
fprintf("Top N Manifold image: ")
disp(manifImgnm(sortId(1:P.topN)))
fprintf("Bottom N Manifold image: " );
disp(manifImgnm(sortId(end-P.botN+1:end)))

if flag.doEvoRef
evoref_vec = cellfun(@(psth)mean(psth(1,51:200,:),'all'),reshape(EStats(Expi).ref.psth_arr,[],1));
evoref_std_vec = cellfun(@(psth)std(mean(psth(1,51:200,:),[1,2]),1,'all'),reshape(EStats(Expi).ref.psth_arr,[],1));
evoref_bsl_VEC = cell2mat(cellfun(@(psth)squeeze(mean(psth(1,1:45,:),[2])),reshape(EStats(Expi).ref.psth_arr,[],1),'uni',0));
bsl_VEC_ALL = [bsl_VEC_ALL; evoref_bsl_VEC];
[evoref_sortScore,sortId] = sort(evoref_vec,'Descend');
[evoref_maxScore,evoref_maxId] = max(evoref_vec);
fprintf("Top N natural reference image: ")
disp(EStats(Expi).ref.imgnm(sortId(1:P.topN)))
fprintf("Bottom N natural reference image: " );
disp(EStats(Expi).ref.imgnm(sortId(end-P.botN+1:end)))
end

% Pasupathy patches
if Stats(Expi).ref.didPasu
pasu_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).ref.pasu_psths',[],1));
pasu_vec(~pasu_val_msk) = []; % isnan(pasu_vec) % get rid of non-existing pasupathy images. 
pasu_std_vec = cellfun(@(psth)std(mean(psth(ui,51:200,:),[1,2]),1,'all'),reshape(Stats(Expi).ref.pasu_psths',[],1));
pasu_std_vec(~pasu_val_msk) = [];
pasu_bsl_VEC = cell2mat(cellfun(@(psth)squeeze(mean(psth(ui,1:45,:),[2])),reshape(Stats(Expi).ref.pasu_psths',[],1),'uni',0));
bsl_VEC_ALL = [bsl_VEC_ALL; pasu_bsl_VEC];
[pasu_sortScore,sortId]=sort(pasu_vec,'Descend');
[pasu_maxScore,pasu_maxId]=max(pasu_vec);
pasuimgnm = cellfun(@(idx)string(unique(Stats(Expi).imageName(idx))),reshape(Stats(Expi).ref.pasu_idx_grid',[],1),'uni',0);
pasuimgnm(~pasu_val_msk) = [];
pasuimgnm = string(pasuimgnm);

fprintf("Top N Pasupathy image: ")
disp(pasuimgnm(sortId(1:P.topN)))
fprintf("Bottom N Pasupathy image: " );
disp(pasuimgnm(sortId(end-P.botN+1:end)))
end
% Gabor patches
if Stats(Expi).ref.didGabor
gab_vec = cellfun(@(psth)mean(psth(ui,51:200,:),'all'),reshape(Stats(Expi).ref.gab_psths,[],1));
gab_std_vec = cellfun(@(psth)std(mean(psth(ui,51:200,:),[1,2]),1,'all'),reshape(Stats(Expi).ref.gab_psths,[],1));
gab_bsl_VEC = cell2mat(cellfun(@(psth)squeeze(mean(psth(ui,1:45,:),[2])),reshape(Stats(Expi).ref.gab_psths,[],1),'uni',0));
bsl_VEC_ALL = [bsl_VEC_ALL; gab_bsl_VEC];
[gab_sortScore,sortId]=sort(gab_vec,'Descend');
[gab_maxScore,gab_maxId]=max(gab_vec);
gaborimgnm = cellfun(@(idx)unique(Stats(Expi).imageName(idx)),reshape(Stats(Expi).ref.gab_idx_grid,[],1));
gaborimgnm = string(gaborimgnm); 

fprintf("Top N Gabor image: ")
disp(gaborimgnm(sortId(1:P.topN)))
fprintf("Bottom N Gabor image: " );
disp(gaborimgnm(sortId(end-P.botN+1:end)))
end
bsl_rate = nanmean(bsl_VEC_ALL);
bsl_rate_std = nanstd(bsl_VEC_ALL,1);
bsl_rate_sem = sem(bsl_VEC_ALL,1);
% if flag.doEvoRef
% saveas(21,fullfile(figdir,compose("%s_Exp%02d_pref%02d_peak_nat.png",Animal,Expi,Stats(Expi).units.pref_chan)))
% savefig(21,fullfile(figdir,compose("%s_Exp%02d_pref%02d_peak_nat.fig",Animal,Expi,Stats(Expi).units.pref_chan)))
% saveas(21,fullfile(figdir,compose("%s_Exp%02d_pref%02d_peak_nat.pdf",Animal,Expi,Stats(Expi).units.pref_chan)))
% else
% saveas(21,fullfile(figdir,compose("%s_Exp%02d_pref%02d_peak.png",Animal,Expi,Stats(Expi).units.pref_chan)))
% savefig(21,fullfile(figdir,compose("%s_Exp%02d_pref%02d_peak.fig",Animal,Expi,Stats(Expi).units.pref_chan)))
% saveas(21,fullfile(figdir,compose("%s_Exp%02d_pref%02d_peak.pdf",Animal,Expi,Stats(Expi).units.pref_chan)))
% end	
end