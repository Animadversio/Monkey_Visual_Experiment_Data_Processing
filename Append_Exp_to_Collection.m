clearvars -except meta_new rasters_new lfps_new Trials_new ExpSpecTable_Aug ExpRecord ExpSpecTable_Aug_alfa
%% Obsolete code to read Formatted data and put similar experiments into a compiled mat file locally. 
ExpSpecTable_Aug = readtable("S:\ExpSpecTable_Augment.xlsx");
%%
% ExpSpecTable = ExpSpecTable_Aug(:, ["ephysFN","expControlFN","stimuli","comments"]);
% writetable(ExpSpecTable, "S:\ExpSpecTable.xlsx")
writetable(ExpSpecTable_Aug, "S:\ExpSpecTable_Augment.xlsx")
%% Load exp file from certain entry 
% [meta_new,rasters_new,lfps_new,Trials_new] = Project_Manifold_Beto_loadRaw([108:113]);
% Animal = "Both";Set_Path;
% expftr = (contains(ExpRecord.expControlFN,"200321") | contains(ExpRecord.expControlFN,"200322"));
% Project_Manifold_Beto_loadRaw(find(expftr),Animal,true);
Animal = "Both";Set_Path;
expftr = (contains(ExpRecord.expControlFN,"200428"));
fllist = find(expftr);
Project_Manifold_Beto_loadRaw(fllist(1:end),Animal,true);
crp_sync_localMat_to_networkMat
%% Load exp file by filtering the experimental record
find(contains(ExpSpecTable_Aug.expControlFN,'generate_parallel'))';
%expid = [2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,45,48,51,55,58,61,64,67,70,72,74,77,79,83,85,88];
expid = find(ExpSpecTable_Aug.Expi > 33 & contains(ExpSpecTable_Aug.expControlFN,'selectivity'));
[meta_new,rasters_new,lfps_new,Trials_new] = Project_Manifold_Beto_loadRaw(expid);
crp_sync_localMat_to_networkMat
%%
[meta_new2,rasters_new2,lfps_new2,Trials_new2] = Project_Manifold_Beto_loadRaw([109,113]);
%
meta_new = [meta_new, meta_new2];
rasters_new = [rasters_new, rasters_new2];
lfps_new = [lfps_new, lfps_new2];
Trials_new = [Trials_new, Trials_new2];
clear meta_new2 rasters_new2 lfps_new2 Trials_new2
%% Code for appending new experiments to the older ones
load("D:\\Manifold_Evolv_Exps.mat")
%%
for i=1:length(meta_new)
    lfps{end+1}=lfps_new{i};
    meta{end+1}=meta_new{i};
    Trials{end+1}=Trials_new{i};
    rasters{end+1}=rasters_new{i}; 
end
savefast("D:\\Manifold_Evolv_Exps.mat", 'lfps', 'meta', 'rasters', 'Trials')
%%
load("D:\\Manifold_Exps.mat")
%%
for i=1
    lfps{end+1}=lfps_new{i};
    meta{end+1}=meta_new{i};
    Trials{end+1}=Trials_new{i};
    rasters{end+1}=rasters_new{i}; 
end
%%
savefast("D:\\Manifold_Exps.mat", 'lfps', 'meta', 'rasters', 'Trials')
%%
storedStruct.meta{4}.stimuli = '\\storage1.ris.wustl.edu\crponce\ActiveStimuli\2019-Manifold\beto-191008a\backup_10_08_2019_12_14_29\PC_imgs';
storedStruct.meta{5}.stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Selectivity\2019-10-09a-beto' ;
storedStruct.meta{6}.stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191010a\backup_10_10_2019_13_10_15\PC_imgs' ;

%%
load("D:\\Manifold_Exps.mat")
meta{4}.stimuli = '\\storage1.ris.wustl.edu\crponce\ActiveStimuli\2019-Manifold\beto-191008a\backup_10_08_2019_12_14_29\PC_imgs';
meta{5}.stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Selectivity\2019-10-09a-beto' ;
meta{6}.stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191010a\backup_10_10_2019_13_10_15\PC_imgs' ;

for Expi = 1: numel(meta_new)/2
lfps{end+1}=lfps_new{2 * Expi - 1};
meta{end+1}=meta_new{2 * Expi - 1};
Trials{end+1}=Trials_new{2 * Expi - 1};
rasters{end+1}=rasters_new{2 * Expi - 1}; 
disp(length(lfps))
end
savefast("D:\\Manifold_Exps.mat", 'lfps', 'meta', 'rasters', 'Trials')
%%
load("D:\\Manifold_Evolv_Exps.mat")
for Expi = 1: numel(meta_new)/2
lfps{end+1}=lfps_new{2 * Expi};
meta{end+1}=meta_new{2 * Expi};
Trials{end+1}=Trials_new{2 * Expi};
rasters{end+1}=rasters_new{2 * Expi}; 
end
savefast("D:\\Manifold_Evolv_Exps.mat", 'lfps', 'meta', 'rasters', 'Trials')
%%
storedStruct.meta = meta;
storedStruct.Trials = Trials;
storedStruct.rasters = rasters;
storedStruct.lfps = lfps;
%% solve the broken experiment problem  Beto64chan-25112019-006 Try to concatanate 2 manifold into 1. 
lfps = cat(3, lfps_new{2}, lfps_new{3});
rasters = cat(3, rasters_new{2}, rasters_new{3});
%%
meta = meta_new{2};
assert(all(meta_new{3}.spikeID==meta_new{2}.spikeID), "The spike ID for 2 experiments matched to each other.")
for field = string(fieldnames(meta))'
    if isstr(getfield(meta,field)) 
        if strcmp(getfield(meta,field), getfield(meta_new{3},field))
            fprintf("%s field matched.\n",field);
            continue;
        end
    else
        if all(getfield(meta,field) == getfield(meta_new{3},field))
            fprintf("%s field matched.\n",field);
            continue;
        end
    end
    meta = setfield(meta, field+"_2", getfield(meta_new{3},field));
end    
meta.comments = "Combined data from 2 broken experiments " + meta.comments;
%% Combine the Trials array
Trials = Trials_new{2};
for field = ["imageName", "XY", "eyeXY", "width", "eyePupil", "eventMarkers", "block"]
    Trials = setfield(Trials, field, [getfield(Trials_new{2},field);getfield(Trials_new{3},field)]);
end
for field = ["B", "trialStart"]
    Trials = setfield(Trials, field, [getfield(Trials_new{2},field),getfield(Trials_new{3},field)]);
end
for field = ["words", "event10", "event03", "MLConfig", "TrialRecord"]
    Trials = setfield(Trials, field+"_2", getfield(Trials_new{3},field));
end
%%
Trials_new{2} = Trials;
lfps_new{2} = lfps;
meta_new{2} = meta;
rasters_new{2} = rasters;
%% Only use it once! 
Trials_new(3)=[];
lfps_new(3)=[];
meta_new(3)=[];
rasters_new(3)=[];
%%
savefast( fullfile(meta.pathMat,[meta.ephysFN '_cmb_formatted']) ,'meta','rasters','lfps','Trials')