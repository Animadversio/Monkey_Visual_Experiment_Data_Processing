%% Catalogue failure experiments

%% the loading failure experiments
setdiff(BigGANTabs_cmb.ephysFN, cellfun(@(M)M.ephysFN,meta_cmb','Unif',0))
% setdiff(cellfun(@(M)M.ephysFN,meta_cmb','Unif',0), BigGANTabs_cmb.ephysFN)
% {'Beto-13062022-003'}
% {'Beto-30112020-002'}

%% The plotting failure experiments (all fixed)
setdiff(arrayfun(@(i)meta_cmb{i}.ephysFN,failure_Expis','UniformOutput',false),...
        cellfun(@(M)M.ephysFN, meta_sing','UniformOutput',false))
% {'Alfa-27072020-002'}
% {'Beto-16102020-003'}
% {'Beto-23072020-002'}

%% Debug: fix the typos in stimuli path 
for Expi = 1:numel(BFEStats)
if strcmp(BFEStats(Expi).meta.ephysFN, "Beto-23072020-002")
    BFEStats(Expi).meta.stimuli = "N:\Stimuli\2020-BigGAN\2020-07-23-Beto-01\2020-07-23-15-59-38";
    meta_cmb{Expi}.stimuli = "N:\Stimuli\2020-BigGAN\2020-07-23-Beto-01\2020-07-23-15-59-38";
    EvolStat = BFEStats(Expi);
    save(fullfile(BFEStats(Expi).meta.figdir,"EvolStat.mat"),'EvolStat')
elseif strcmp(BFEStats(Expi).meta.ephysFN, "Alfa-27072020-002")
    BFEStats(Expi).meta.stimuli = "N:\Stimuli\2020-BigGAN\2020-07-27-Alfa-01\2020-07-27-09-47-40";
    meta_cmb{Expi}.stimuli = "N:\Stimuli\2020-BigGAN\2020-07-27-Alfa-01\2020-07-27-09-47-40";
    EvolStat = BFEStats(Expi);
    save(fullfile(BFEStats(Expi).meta.figdir,"EvolStat.mat"),'EvolStat')
else
    continue
end
end

%% Debug Alfa-27072020-002 experiments: Add missing stimuli

load('N:\Stimuli\2020-BigGAN\2020-07-27-Alfa-01\2020-07-27-09-47-40\block026_thread000_code.mat')
imgs = G.visualize(codes);
figure;montage(imgs)
%%
missing_ids = find(strcmp("block027_thread000_gen_gen026_001699.bmp",ids));
img = G.visualize(codes(missing_ids,:));
imwrite(img,fullfile('N:\Stimuli\2020-BigGAN\2020-07-27-Alfa-01\2020-07-27-09-47-40',...
    "block027_thread000_gen_gen026_001699.bmp"))
%% Check the stimuli exists
for imgi = 1:numel(BFEStats(92).imageName)
    imgnm = BFEStats(92).imageName{imgi};
    if ~exist(fullfile('N:\Stimuli\2020-BigGAN\2020-07-27-Alfa-01\2020-07-27-09-47-40',...
            [imgnm,'.bmp']),'file')
    elseif ~exist(fullfile('N:\Stimuli\2020-BigGAN\2020-07-27-Alfa-01',...
            [imgnm(1:end-14),'.bmp']),'file')
    else
        fprintf([imgnm,'\n'])
        keyboard
    end
end


%% Debug  Beto-16102020-003  missing trial record: manual add trial record
% BigGANTabs = readtable("Evol_BigGAN_AB_final.xlsx",'Format','auto');
Expi = find(BigGANTabs.ephysFN == "Beto-16102020-003"); % 48
assert(isempty(Trials_cmb{48}.TrialRecord)) % this experiment failed to save thre TrialRecord
%%
load(fullfile(meta_cmb{48}.stimuli,"space_opts.mat"),'space_opts');
User = struct();
User.evoConfiguration = {5   [-0.5, 0.7]  2  1  'CMAES_Hessian'
                         5   [-0.5, 0.7]  2  1  'CMAES'};
User.prefChan = [5, 5];
User.space = {space_opts.name}'; % 
User.space_cfg = {{space_opts(1).name  space_opts(1).setting}
                  {space_opts(2).name  space_opts(2).setting}};
User.space_opts = space_opts;
Trials_cmb{48}.TrialRecord.User = User;
% 003 FC6 BigGAN evol 
% FC6 CMAES_Hessian ,   ch5 2 (-0.5 0.7) 
% BigGAN CMAES,   ch5 2 (-0.5 0.7) 
% This CMAES only cut off Spectrum no rescaling!!! Use this version to save time. 
% 4 mins finish 5 blocks pretty efficient  
% Taking a nap!  
% CMAES-Hessian goes up immediately.  
% Both climb .  
% 11block 10 mins.  CMA Hessian does well on the FC6 
% ~ 21 blocks. BigGAN also climbs a bit. Not too dramatic but steady. 
% Terminate ~ 30. Very efficient here!  
% During Saving the matlab crushed....  
% So TrialRecord is missing. Need to fix this bug!!!!! 
% Finished 
%%
crp_sync_localMat_to_networkMat
%%
BFEStats_fix = Evol_BigGAN_FC6_Collect_Stats_fun(meta_cmb(48), rasters_cmb(48), Trials_cmb(48));
% BFEStats(48) units is empty
BFEStats(48) = BFEStats_fix;
%%
load('S:\Data-Ephys-MAT\Beto-16102020-003_formatted.mat')
%%
load(fullfile('N:\Stimuli\2020-BigGAN\2020-10-16-Beto-01\2020-10-16-11-21-49',"space_opts.mat"),'space_opts');
User = struct();
User.evoConfiguration = {5   [-0.5, 0.7]  2  1  'CMAES_Hessian'
                         5   [-0.5, 0.7]  2  1  'CMAES'};
User.prefChan = [5, 5];
User.space = {space_opts.name}'; % 
User.space_cfg = {{space_opts(1).name  space_opts(1).setting}
                  {space_opts(2).name  space_opts(2).setting}};
User.space_opts = space_opts;

Trials.TrialRecord.User = User;
%%
savefast('S:\Data-Ephys-MAT\Beto-16102020-003_formatted.mat','Trials','lfps','meta','rasters')
%% Debug  Alfa-16112020-001 mis aligned bhv2 and PL2 
% correct matches
%       201116_Alfa_generate_BigGAN(1)
%       Alfa-16112020-001
% but these are not the exp described by the comments.

% find(strcmp("Alfa-16112020-001",arrayfun(@(S)S.meta.ephysFN,BFEStats,'Unif',0))); % 154
BigGANTabs = readtable("Evol_BigGAN_AB_final.xlsx",'Format','auto');
Expi = find(BigGANTabs.expControlFN=="201116_Alfa_generate_BigGAN(1)");
preMeta = table2preMeta(BigGANTabs(Expi,:));
preMeta.ephysFN = 'Alfa-16112020-001';
[meta_subst, rasters_subst, lfps_subst, Trials_subst] = loadExperiments_preMeta(preMeta, false);



%%
function preMeta = table2preMeta(Tabs)
iExp = 0;
for iExp = 1:size(Tabs,1)
%     rowi = rowlist(iExp);
    preMeta(iExp).ephysFN = Tabs.ephysFN{iExp}; 
    preMeta(iExp).expControlFN = Tabs.expControlFN{iExp}; % 
    preMeta(iExp).stimuli = Tabs.stimuli{iExp} ;
    preMeta(iExp).comments = Tabs.comments{iExp};
end
end