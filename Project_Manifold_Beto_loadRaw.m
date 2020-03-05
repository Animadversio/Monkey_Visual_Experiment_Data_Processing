function [meta,rasters,lfps,Trials] = Project_Manifold_Beto_loadRaw(rowlist, animal)
if nargin == 1
animal = "Beto";
end
switch animal 
    case "Beto"
        ExpSpecTable_Aug = readtable("S:\ExpSpecTable_Augment.xlsx");
    case "Alfa"
        ExpSpecTable_Aug = readtable("Exp_Record_Alfa.xlsx");
end
% ExpSpecTable_Aug = readtable("S:\ExpSpecTable_Augment.xlsx");
iExp = 0;
for iExp = 1:numel(rowlist)
    rowi = rowlist(iExp);
    preMeta(iExp).ephysFN = ExpSpecTable_Aug.ephysFN{rowi}; 
    preMeta(iExp).expControlFN = ExpSpecTable_Aug.expControlFN{rowi}; % 
    preMeta(iExp).stimuli = ExpSpecTable_Aug.stimuli{rowi} ;
    preMeta(iExp).comments = ExpSpecTable_Aug.comments{rowi};
end

Project_General_copyMissingFiles(preMeta); % communicating and copying data from network to local 

for iExp = 1:length(preMeta) 
    
    tMeta = preMeta(iExp);
    [meta_,rasters_,lfps_,Trials_] = loadData(tMeta.ephysFN,'expControlFN',tMeta.expControlFN) ;
    meta_merged = rmfield( tMeta, intersect(fieldnames(tMeta), fieldnames(meta_)) );
    names = [fieldnames(meta_merged); fieldnames(meta_)];
    meta_ = cell2struct([struct2cell(meta_merged); struct2cell(meta_)], names, 1);

    meta{iExp} = meta_;
    rasters{iExp} = rasters_;
    lfps{iExp} = lfps_;
    Trials{iExp} = Trials_;
    clear meta_  rasters_ lpfs_ Trials_ names meta_merged tMeta
end
end
