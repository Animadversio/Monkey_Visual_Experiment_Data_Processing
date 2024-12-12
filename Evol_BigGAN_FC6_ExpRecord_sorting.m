%% Evol BigGAN related experiments filtering

rootdir = "E:\OneDrive - Harvard University\Evol_BigGAN_FC6";
mat_dir = "O:\Mat_Statistics"; 
saveroot = "O:\Evol_BigGAN_FC6_cmp"; 
%%
Animal = "Both";
ExpRecord = readExpRecord("Both");
%% Rough filter logic 
%  1) has to be a generate experiment
%  2) cosine evolution is not , movie is not
%  3) it's after 2020 07 19, before that I have not invented this exp.
ftrrows = find(contains(ExpRecord.expControlFN,"generate_") & ...
               (~ contains(ExpRecord.expControlFN,"cosine"))  & ...
               (~ contains(ExpRecord.expControlFN,"Movie"))  & ...
               cellfun(@(date)date > datetime(2020,7,19),ExpRecord.expdate));
% disp(ExpRecord(ftrrows,:))
%%
%% Additional experiments from carlos' notes
ftrrows_add = find(contains(ExpRecord.expControlFN,"generate_") & ...
               (~ contains(ExpRecord.expControlFN,"cosine"))  & ...
               (~ contains(ExpRecord.expControlFN,"Movie"))  & ...
               ...
               ((ExpRecord.expdate == datetime(2020,7,24)) | ...
               (ExpRecord.expdate == datetime(2020,9,28)) | ...
               (ExpRecord.expdate == datetime(2021,5,28))) );
%% Do some manual work to sort them. 
...

%% exp from the first filter 
BigGANTabs = readtable("Evol_BigGAN_AB.xlsx",'Format','auto');
%%
% [BigGANTabs;ExpRecord(ftrrows_add,:)]
BigGANTabs_cmb = [BigGANTabs;ExpRecord(ftrrows_add,:)];
%% combined 
writetable([BigGANTabs;ExpRecord(ftrrows_add,:)],"Evol_BigGAN_AB_final.xlsx")
%%
BigGAN_Alfa = BigGANTabs(contains(BigGANTabs.ephysFN,"Alfa"), :);
BigGAN_Beto = BigGANTabs(contains(BigGANTabs.ephysFN,"Beto"), :);
%%
% preMeta = struct();
preMeta = table2preMeta(BigGANTabs);
%%
% [meta_new, rasters_new, lfps_new, Trials_new] = loadExperiments(ftrrows([43,45,46:end]), Animal, false);%43,45
[meta_new, rasters_new, lfps_new, Trials_new] = loadExperiments_preMeta(preMeta, false);
%% Process the new experiments
preMeta = table2preMeta(ExpRecord(ftrrows_add,:));
% [meta_new, rasters_new, lfps_new, Trials_new] = loadExperiments(ftrrows([43,45,46:end]), Animal, false);%43,45
[meta_add, rasters_add, lfps_add, Trials_add] = loadExperiments_preMeta(preMeta, false);
%%
meta_cmb = [meta_new,meta_add];
rasters_cmb = [rasters_new,rasters_add];
Trials_cmb = [Trials_new,Trials_add];
%% Try to see any failed experiments 
for iExp = 1:size(meta_new)
    if isempty(Trials_cmb{iExp}) || isempty(rasters_cmb{iExp}) 
        disp(meta_new{iExp})
    end
end
%% Save the failed loading exp to disk.
eptymsk = cellfun(@isempty,Trials_new);
disp(sum(eptymsk))
writetable(BigGANTabs(eptymsk,:),"S:\Evol_BigGAN_AB_loadingfail.xlsx")
meta_new(eptymsk)=[];
rasters_new(eptymsk)=[];
Trials_new(eptymsk)=[];

%% New processing Function
% BFEStats = Evol_BigGAN_FC6_Collect_Stats_fun(meta_new(1:end), rasters_new(1:end), Trials_new(1:end));
% New Function
BFEStats = Evol_BigGAN_FC6_Collect_Stats_fun(meta_cmb, rasters_cmb, Trials_cmb);


%% Sort out single exps 
[meta_sing, rasters_sing, Trials_sing] = Evol_BigGAN_FC6_find_singular_exps(meta_cmb, rasters_cmb, Trials_cmb);
%%
savefast(fullfile("S:\Evol_BigGAN_singular_ExpsMeta.mat"),'meta_sing')
%%
BFEStats_sing = Evol_BigGAN_FC6_singular_Collect_Stats_fun(meta_sing, rasters_sing, Trials_sing);
%%
for Expi = 1:numel(BFEStats_sing)
if isempty(BFEStats_sing(Expi).evol) && Expi == 17
    continue
end
fprintf("Id%02d  %s - %s - chan%02d(U%d) %s [%.1f %.1f]\n",Expi,BFEStats_sing(Expi).Animal,...
                  BFEStats_sing(Expi).meta.ephysFN,...
                  BFEStats_sing(Expi).evol.pref_chan,...
                  BFEStats_sing(Expi).evol.unit_in_pref_chan,...
                  string(BFEStats_sing(Expi).evol.space_names),...
                  BFEStats_sing(Expi).evol.imgpos(1), BFEStats_sing(Expi).evol.imgpos(2))
    
end
%%
sing_pairing = [[ 2, 1]
                [ 9,10]
                [ 9,11]
                [ 9,12]
                [13,14]
                [15,"Alfa-29072020-003"] % pair with a thread in a paired exp
                [16,"Alfa-13082020-002"] % pair with a thread in a paired exp
                [18,"Alfa-07012021-002"]];
clearly_excluded = [4,5,6,7,8,17,19,20]; % pure fc6 no pairing
%4 Beto-13042021-002  the paired BigGAN is note saved.
%3 Beto-29072020-004  no pairing just BigGAN

% 17 is this  'Alfa-16112020-001' needs to be excluded. 
% 19 is this  'Beto-28052021-002' not a pair 
% 20 is this  'Beto-28052021-003' not a pair 

%% the loading failure experiments
setdiff(BigGANTabs_cmb.ephysFN, cellfun(@(M)M.ephysFN,meta_cmb','Unif',0))
% setdiff(cellfun(@(M)M.ephysFN,meta_cmb','Unif',0), BigGANTabs_cmb.ephysFN)
%% The plotting failure experiments 
setdiff(arrayfun(@(i)meta_cmb{i}.ephysFN,failure_Expis','UniformOutput',false),...
        cellfun(@(M)M.ephysFN, meta_sing','UniformOutput',false))
% {'Alfa-27072020-002'}
% {'Beto-16102020-003'}
% {'Beto-23072020-002'}
%% Experiments found to be baseline fluctuating. 
baseline_excl = ["Beto-18082020-002",
                 "Beto-07092020-006",
                 "Beto-14092020-002",
                 "Beto-27102020-003",
                 "Alfa-22092020-003"];
sum([BFEStats.Animal]=="Alfa" & paired_mask)
sum([BFEStats.Animal]=="Beto" & paired_mask)
paired_mask = arrayfun(@(S)~isempty(S.evol), BFEStats);
%%
% find(strcmp("Beto-23072020-002",arrayfun(@(S)S.meta.ephysFN,BFEStats,'Unif',0))); % 5
% BFEStats(Expi).meta.stimuli = "N:\Stimuli\2020-BigGAN\2020-07-23-Beto-01\2020-07-23-15-59-38";
% find(strcmp("Alfa-27072020-002",arrayfun(@(S)S.meta.ephysFN,BFEStats,'Unif',0))) % 92
% "N:\Stimuli\2020-BigGAN\2020-07-27-Alfa-01\2020-07-27-09-47-40"

%%
savefast(fullfile(mat_dir,"Both_BigGAN_FC6_Evol_Stats.mat"),'BFEStats')
savefast(fullfile(mat_dir,"Both_BigGAN_FC6_singular_Evol_Stats.mat"),'BFEStats_sing')
%%
sum(arrayfun(@(S)isempty(S.units),BFEStats))
%%
empty_exp = {};
for i = 1:numel(BFEStats)
if isempty(BFEStats(i).units)
empty_exp{end+1} =  BFEStats(i).meta.ephysFN;
end
end
empty_exp = cat(1,empty_exp{:});
%%
singular_exps = cellfun(@(M)M.ephysFN,meta_sing,'UniformOutput',false)';
%%
setdiff(empty_exp,singular_exps)
%%
setdiff(singular_exps,empty_exp)
%%
meta_cmb = [meta_new,meta_add];
for i = 1:numel(meta_cmb)
    if strcmp(meta_cmb{i}.ephysFN,'Beto-16102020-003')
        disp(meta_cmb{i})
        break
    end
end
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