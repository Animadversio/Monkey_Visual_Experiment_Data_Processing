function [meta,rasters,lfps,Trials] = loadExperiments(rowlist, animal, no_return, no_lfp)
if nargin == 1
animal = "Beto";
no_return = false;
no_lfp = false;
elseif nargin == 2
no_return = false;
no_lfp = true;
elseif nargin == 3
no_lfp = true;
end
switch animal 
    case "Beto"
        ExpRecord = readtable("S:\ExpSpecTable_Augment.xlsx");
    case "Alfa"
        ExpRecord = readtable("S:\Exp_Record_Alfa.xlsx");
    case "Both"
        ExpSpecTable_Aug = readtable("S:\ExpSpecTable_Augment.xlsx");
        ExpSpecTable_Aug_alfa = readtable("S:\Exp_Record_Alfa.xlsx");
        ExpRecord = [ExpSpecTable_Aug; ExpSpecTable_Aug_alfa];
end
% ExpSpecTable_Aug = readtable("S:\ExpSpecTable_Augment.xlsx");
iExp = 0;
for iExp = 1:numel(rowlist)
    rowi = rowlist(iExp);
    preMeta(iExp).ephysFN = ExpRecord.ephysFN{rowi}; 
    preMeta(iExp).expControlFN = ExpRecord.expControlFN{rowi}; % 
    preMeta(iExp).stimuli = ExpRecord.stimuli{rowi} ;
    preMeta(iExp).comments = ExpRecord.comments{rowi};
end

Project_General_copyMissingFiles(preMeta); % communicating and copying data from network to local 
meta = {}; rasters = {}; lfps = {}; Trials = {};
for iExp = 1:length(preMeta) 
    tMeta = preMeta(iExp);
    try
    [meta_,rasters_,lfps_,Trials_] = loadData(tMeta.ephysFN,'expControlFN',tMeta.expControlFN) ;
    catch err %e is an MException struct
        fprintf('Error message:\n%s\n',err.message);
        fprintf('Error trace:\n%s\n',err.getReport);
        disp(tMeta)
        %keyboard
        fileID = fopen('S:\Exp_error_log.log','w+');
        fprintf(fileID,'Error message:\n%s\n',err.message);
        fprintf(fileID,'Error trace:\n%s\n',err.getReport);
        fclose(fileID);
        continue
    end
    meta_merged = rmfield( tMeta, intersect(fieldnames(tMeta), fieldnames(meta_)) );
    names = [fieldnames(meta_merged); fieldnames(meta_)];
    meta_ = cell2struct([struct2cell(meta_merged); struct2cell(meta_)], names, 1);
    if ~no_return % if true then don't return these things, only save to disk. 
        % helpful when we don't want to take up all the memory! 
    meta{iExp} = meta_;
    rasters{iExp} = rasters_;
    if ~no_lfp
    lfps{iExp} = lfps_;
    end
    Trials{iExp} = Trials_;
    end
    clear meta_  rasters_ lpfs_ Trials_ names meta_merged tMeta
end
end
