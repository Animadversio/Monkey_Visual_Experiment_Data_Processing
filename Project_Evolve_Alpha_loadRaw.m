function [meta,rasters,lfps,Trials] = Project_Evolve_Alpha_loadRaw

% day 001 
% iExp = 1; % PC space exploration
% preMeta(iExp).ephysFN = 'Beto64chan-02102019-003'; 
% preMeta(iExp).expControlFN = '191002_Beto_selectivity_basic(1)'; % '191002_Beto_selectivity_basic(1).bhv2';
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Selectivity\2019-10-02a-beto' ;
% preMeta(iExp).comments = 'Derived from PCs of Ch 29 evolution 191002. 373 images {11*11 in PC2,PC3; PC49 PC50; RND RND space respectively, +10 last gen}';
% 
% iExp = 1; % Cma
% preMeta(iExp).ephysFN = 'Beto64chan-02102019-001'; 
% preMeta(iExp).expControlFN = '191002_Beto_generate_parallel'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191002a\backup_10_02_2019_14_33_18' ;
% preMeta(iExp).comments = 'CMA Evolution of Ch 29 191002. generate PCs for later PC space tuning';
% day 002
iExp = 1; % PC space exploration
preMeta(iExp).ephysFN = 'Beto64chan-03102019-002'; 
preMeta(iExp).expControlFN = '191003_Beto_selectivity_basic'; 
preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Selectivity\2019-10-03a-beto' ;
preMeta(iExp).comments = ['Derived from PCs of Ch 6 evolution 191003(CRP). Fewer Reps 408 images {11*11 in PC2,PC3; PC49 PC50; RND12 space respectively,',...
                          ' +10 last gen +12 gabor + 23 reference image}'];

iExp = 2; % Cma
preMeta(iExp).ephysFN = 'Beto64chan-03102019-001'; 
preMeta(iExp).expControlFN = '191003_Beto_generate_parallel'; % 
preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191003a\backup_10_03_2019_13_42_32' ;
preMeta(iExp).comments = 'CMA Evolution of Ch 6 191003(CRP). generate PCs for later PC space tuning';

% iExp = 2; preMeta(iExp).ephysFN = 'monkey64chan-11062019-002'; preMeta(iExp).expControlFN = '190611_Beto_generate_simple(1)'; preMeta(iExp).stimuli = 'n:\Stimuli\2019-06-Evolutions\beto-190611b' ;

% iExp = ; % preMeta(iExp).ephysFN = ''; % preMeta(iExp).expControlFN = ''; % preMeta(iExp).stimuli = fullfile('n:',' ' ) ;

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