function [meta,rasters,lfps,Trials] = Project_Manifold_Beto_loadRaw
iExp = 0;
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
% % day 002
% % iExp = 2; % PC space exploration
% % preMeta(iExp).ephysFN = 'Beto64chan-03102019-002'; 
% % preMeta(iExp).expControlFN = '191003_Beto_selectivity_basic'; 
% % preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Selectivity\2019-10-03a-beto' ;
% % preMeta(iExp).comments = ['Derived from PCs of Ch 6 evolution 191003(CRP). Fewer Reps 408 images {11*11 in PC2,PC3; PC49 PC50; RND12 space respectively,',...
% %                           ' +10 last gen +12 gabor + 23 reference image}'];
% % % 
% iExp = 2; % Cma
% preMeta(iExp).ephysFN = 'Beto64chan-03102019-001'; 
% preMeta(iExp).expControlFN = '191003_Beto_generate_parallel'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191003a\backup_10_03_2019_13_42_32' ;
% preMeta(iExp).comments = 'CMA Evolution of Ch 6 191003(CRP). generate PCs for later PC space tuning';
% % % day 003
% % iExp = 3; % PC space exploration
% % preMeta(iExp).ephysFN = 'Beto64chan-07102019-005'; 
% % preMeta(iExp).expControlFN = '191007_Beto_selectivity_basic(1)'; 
% % preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Selectivity\2019-10-07a-beto' ;
% % preMeta(iExp).comments = ['Derived from PCs of Ch 5 evolution 191007(CRP). Fewer Reps 373 images {11*11 in PC2,PC3; PC49 PC50; RND12 space respectively,'];
% 
% iExp = 3; % Cma
% preMeta(iExp).ephysFN = 'Beto64chan-07102019-003'; 
% preMeta(iExp).expControlFN = '191007_Beto_generate_parallel(2)'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191007a\backup_10_07_2019_14_03_46' ;
% preMeta(iExp).comments = 'CMA Evolution of Ch 5 191007(CRP). generate PCs for later PC space tuning';
% % % day 004
% % iExp = 4; % PC space exploration
% % preMeta(iExp).ephysFN = 'Beto64chan-08102019-006'; 
% % preMeta(iExp).expControlFN = '191008_Beto_selectivity_basic'; 
% % preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191008a\backup_10_08_2019_12_14_29\PC_imgs' ;
% % preMeta(iExp).comments = ['Derived from PCs of Ch 20 evolution 191008(CRP). Fewer Reps 363 images {11*11 in PC2,PC3; PC49 PC50; RND12 space respectively,'];
% % 
% iExp = 4; % Cma
% preMeta(iExp).ephysFN = 'Beto64chan-08102019-004'; 
% preMeta(iExp).expControlFN = '191008_Beto_generate_parallel(2)'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191008a\backup_10_08_2019_12_14_29' ;
% preMeta(iExp).comments = 'CMA Evolution of Ch 20 191008(CRP). generate PCs for later PC space tuning';
% % 
% % % day 005
% % iExp = 5; % PC space exploration
% % preMeta(iExp).ephysFN = 'Beto64chan-09102019-003'; 
% % preMeta(iExp).expControlFN = '191009_Beto_selectivity_basic'; 
% % preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Selectivity\2019-10-09a-beto' ;
% % preMeta(iExp).comments = ['Derived from PCs of Ch 19 evolution 191009(BXW). Fewer Reps 363 images {11*11 in PC2,PC3; PC49 PC50; RND12 space respectively,'];
% % % 
% iExp = 5; % Cma
% preMeta(iExp).ephysFN = 'Beto64chan-09102019-002'; 
% preMeta(iExp).expControlFN = '191009_Beto_generate_parallel(1)'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191009a\backup_10_09_2019_13_49_17' ;
% preMeta(iExp).comments = 'CMA Evolution of Ch 19 191009(CRP). generate PCs for later PC space tuning';
% 
% % day 006
% % iExp = 6; % PC space exploration
% % preMeta(iExp).ephysFN = 'Beto64chan-10102019-003'; 
% % preMeta(iExp).expControlFN = '191010_Beto_selectivity_basic'; 
% % preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191010a\backup_10_10_2019_13_10_15\PC_imgs' ;
% % preMeta(iExp).comments = ['Derived from PCs of Ch 13 evolution 191010(CRP). Fewer Reps 363 images {11*11 in PC2,PC3; PC49 PC50; RND12 space respectively,'];
% % 
% iExp = 6; % Cma
% preMeta(iExp).ephysFN = 'Beto64chan-10102019-002'; 
% preMeta(iExp).expControlFN = '191010_Beto_generate_parallel(1)'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191010a\backup_10_10_2019_13_10_15' ;
% preMeta(iExp).comments = 'CMA Evolution of Ch 13 191009(CRP). generate PCs for later PC space tuning';
% % day 007
% % iExp = 7; % PC space exploration
% % preMeta(iExp).ephysFN = 'Beto64chan-11102019-002'; 
% % preMeta(iExp).expControlFN = '191011_Beto_selectivity_basic'; 
% % preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191011a\backup_10_11_2019_13_17_02\PC_imgs' ;
% % preMeta(iExp).comments = ['Derived from PCs of Ch 28 evolution 191011(CRP). Fewer Reps 363 images {11*11 in PC2,PC3; PC49 PC50; RND12 space respectively,'];
% % 
% iExp = 7; % Cma
% preMeta(iExp).ephysFN = 'Beto64chan-11102019-001'; 
% preMeta(iExp).expControlFN = '191011_Beto_generate_parallel'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191011a\backup_10_11_2019_13_17_02' ;
% preMeta(iExp).comments = 'CMA Evolution of Ch 28 191011(CRP). generate PCs for later PC space tuning';
% 
% % day 008
% % iExp = 8; % PC space exploration
% % preMeta(iExp).ephysFN = 'Beto64chan-15102019-005'; 
% % preMeta(iExp).expControlFN = '191015_Beto_selectivity_basic'; 
% % preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191015a\backup_10_15_2019_14_25_42\PC_imgs' ;
% % preMeta(iExp).comments = ['Derived from PCs of Ch 15 evolution 191015(CRP BXW). 363 images {11*11 in PC2,PC3; PC49 PC50; RND12 space respectively,'];
% % 
% iExp = 8; % Cma
% preMeta(iExp).ephysFN = 'Beto64chan-15102019-003'; 
% preMeta(iExp).expControlFN = '191015_Beto_generate_parallel(2)'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191015a\backup_10_15_2019_14_25_42' ;
% preMeta(iExp).comments = 'CMA Evolution of Ch 15 191015(CRP BXW). generate PCs for later PC space tuning';
% 
% % day 009
% % iExp = 9; % PC space exploration
% % preMeta(iExp).ephysFN = 'Beto64chan-16102019-002'; 
% % preMeta(iExp).expControlFN = '191016_Beto_selectivity_basic'; 
% % preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191016a\backup_10_16_2019_14_20_03\PC_imgs' ;
% % preMeta(iExp).comments = ['Derived from PCs of Ch 17 evolution 191016(CRP). 363 images {11*11 in PC2,PC3; PC49 PC50; RND12 space respectively,'];
% % 
% iExp = 9; % Cma
% preMeta(iExp).ephysFN = 'Beto64chan-16102019-001'; 
% preMeta(iExp).expControlFN = '191016_Beto_generate_parallel'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191016a\backup_10_16_2019_14_20_03' ;
% preMeta(iExp).comments = 'CMA Evolution of Ch 17 191016(CRP). generate PCs for later PC space tuning';
% 
% % day 010
% % iExp = 9; % PC space exploration
% % preMeta(iExp).ephysFN = 'Beto64chan-17102019-007'; 
% % preMeta(iExp).expControlFN = '191017_Beto_selectivity_basic'; 
% % preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191017a\backup_10_17_2019_15_39_04\PC_imgs' ;
% % preMeta(iExp).comments = ['Derived from PCs of Ch 63 evolution 191017(CRP). 363 images {11*11 in PC2,PC3; PC49 PC50; RND12 space respectively,'];
% 
% iExp = 10; % CMA
% preMeta(iExp).ephysFN = 'Beto64chan-17102019-006'; 
% preMeta(iExp).expControlFN = '191017_Beto_generate_parallel(3)'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191017a\backup_10_17_2019_15_39_04' ;
% preMeta(iExp).comments = 'CMA Evolution of Ch 63 191017(CRP). generate PCs for later PC space tuning';

% day 011
iExp = iExp + 1; %iExp = 11; % PC space exploration
preMeta(iExp).ephysFN = 'Beto64chan-30102019-003'; 
preMeta(iExp).expControlFN = '191030_Beto_selectivity_basic'; 
preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191030a-selectivity' ;
preMeta(iExp).comments = ['Following Evolution of Ch26 Explore in the PC2-3 space only and add the pasupathy images to it. white background 4 rotations.'];

iExp = iExp + 1; %iExp = 11; % CMA
preMeta(iExp).ephysFN = 'Beto64chan-30102019-001'; 
preMeta(iExp).expControlFN = '191030_Beto_generate_parallel'; % 
preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191030a\backup_10_30_2019_10_15_31' ;
preMeta(iExp).comments = 'CMA Evolution of Ch 26 191030(CRP). generate PCs for later PC space tuning';
% 
% % % day 012
iExp = iExp + 1; %iExp = 12; % PC space exploration
preMeta(iExp).ephysFN = 'Beto64chan-31102019-004'; 
preMeta(iExp).expControlFN = '191031_Beto_selectivity_basic'; 
preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191031a-selectivity' ;
preMeta(iExp).comments = ['Derived from PCs of Ch 20 evolution 191031(CRP). 121 images from PC23 space, add the pasupathy images to it. white background 4 rotations.'];

iExp = iExp + 1;% iExp = 12; % CMA
preMeta(iExp).ephysFN = 'Beto64chan-31102019-002'; 
preMeta(iExp).expControlFN = '191031_Beto_generate_parallel'; % 
preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191031a\backup_10_31_2019_13_24_49' ;
preMeta(iExp).comments = 'CMA Evolution of Ch 20 191031(CRP). generate PCs for later PC space tuning';

% % day 013
iExp = iExp + 1; %iExp = 13; % PC space exploration
preMeta(iExp).ephysFN = 'Beto64chan-05112019-004'; 
preMeta(iExp).expControlFN = '191105_Beto_selectivity_basic'; 
preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191105a-selectivity' ;
preMeta(iExp).comments = ['Derived from PCs of Ch 40 evolution 191105(CRP). 121 images from PC23 space, add the pasupathy images to it. white background 4 rotations.'];

iExp = iExp + 1; %iExp = 13; % CMA
preMeta(iExp).ephysFN = 'Beto64chan-05112019-003'; 
preMeta(iExp).expControlFN = '191105_Beto_generate_parallel(1)'; % 
preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191105a\backup_11_05_2019_11_23_09' ;
preMeta(iExp).comments = 'CMA Evolution of Ch 40 evolution 191105(CRP). generate PCs for later PC space tuning';

% 
% % day 014
iExp = iExp + 1; %iExp = 14; % PC space exploration
preMeta(iExp).ephysFN = 'Beto64chan-06112019-005'; 
preMeta(iExp).expControlFN = '191106_Beto_selectivity_basic'; 
preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191106a-selectivity' ;
preMeta(iExp).comments = ['Derived from PCs of Ch 42 evolution 191106(CRP). 121 images from PC23 space, add the pasupathy images to it. white background 4 rotations.'];

iExp = iExp + 1; %iExp = 14; % CMA
preMeta(iExp).ephysFN = 'Beto64chan-06112019-003'; 
preMeta(iExp).expControlFN = '191106_Beto_generate_parallel'; % 
preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191106a\backup_11_06_2019_13_12_32' ;
preMeta(iExp).comments = 'CMA Evolution of Ch 42 evolution 191106(CRP). generate PCs for later PC space tuning';

iExp = iExp + 1; %iExp = 15; % PC space exploration
preMeta(iExp).ephysFN = 'Beto64chan-06112019-008'; 
preMeta(iExp).expControlFN = '191106_Beto_selectivity_basic(1)'; 
preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191106b-selectivity' ;
preMeta(iExp).comments = ['Derived from PCs of Ch 35 evolution 191106(CRP). 121 images from PC23 space, add the pasupathy images to it. white background 4 rotations.'];

iExp = iExp + 1; %iExp = 15; % CMA 
preMeta(iExp).ephysFN = 'Beto64chan-06112019-007'; 
preMeta(iExp).expControlFN = '191106_Beto_generate_parallel(1)'; % 
preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191106b\backup_11_06_2019_14_17_35' ;
preMeta(iExp).comments = 'CMA Evolution of Ch 35 evolution 191106(CRP). generate PCs for later PC space tuning';
% day 015
iExp = iExp + 1; % iExp = 16; % PC space exploration
preMeta(iExp).ephysFN = 'Beto64chan-07112019-006'; 
preMeta(iExp).expControlFN = '191107_Beto_selectivity_basic'; 
preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191107a-selectivity' ;
preMeta(iExp).comments = ['Derived from PCs of Ch 36 evolution 191107(CRP). 121 images from PC23 space, add the pasupathy images to it. white background 4 rotations.'];

iExp = iExp + 1; % iExp = 16; % CMA
preMeta(iExp).ephysFN = 'Beto64chan-07112019-004'; 
preMeta(iExp).expControlFN = '191107_Beto_generate_parallel'; % 
preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191107a\backup_11_07_2019_12_29_28' ;
preMeta(iExp).comments = 'CMA Evolution of Ch 36 evolution 191107(CRP). generate PCs for later PC space tuning';

iExp = iExp + 1;% iExp = 17; % PC space exploration
preMeta(iExp).ephysFN = 'Beto64chan-07112019-012'; 
preMeta(iExp).expControlFN = '191107_Beto_selectivity_basic(1)'; %?????
preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191107b-selectivity' ;
preMeta(iExp).comments = ['Derived from PCs of Ch 34 evolution 191107(CRP). 121 images from PC23 space, add the pasupathy images to it. white background 4 rotations.'];

iExp = iExp + 1;% iExp = 17; % CMA
preMeta(iExp).ephysFN = 'Beto64chan-07112019-010'; 
preMeta(iExp).expControlFN = '191107_Beto_generate_parallel(2)'; % 
preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191107b\backup_11_07_2019_13_50_18' ;
preMeta(iExp).comments = 'CMA Evolution of Ch 34 evolution 191107(CRP). generate PCs for later PC space tuning';

% day 016
iExp = iExp + 1;% iExp = 18; % PC space exploration
preMeta(iExp).ephysFN = 'Beto64chan-08112019-004'; 
preMeta(iExp).expControlFN = '191108_Beto_selectivity_basic'; %?????
preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191108a-selectivity' ;
preMeta(iExp).comments = ['Derived from PCs of Ch 37 evolution 191108(CRP). 121 images from PC23 space, add the pasupathy images to it. white background 4 rotations.'];

iExp = iExp + 1;% iExp = 18; % CMA
preMeta(iExp).ephysFN = 'Beto64chan-08112019-003'; 
preMeta(iExp).expControlFN = '191108_Beto_generate_parallel'; % 
preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191108a\backup_11_08_2019_10_05_21' ;
preMeta(iExp).comments = 'CMA Evolution of Ch 37 evolution 191108(CRP). generate PCs for later PC space tuning';

% day 001 Interp Exp
% iExp = 1; % PC space interpolation between the 2 ending points 
% preMeta(iExp).ephysFN = 'Beto64chan-14102019-007'; 
% preMeta(iExp).expControlFN = '191014_Beto_selectivity_basic(1)'; 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191014a\backup_10_14_2019_12_51_15\PC_imgs' ;
% preMeta(iExp).comments = ['Derived from PCs of Ch 26 single unit evolution 191014(CRP BXW). (16*11 - 30)*2= 292 images'];
% 
% iExp = 2; % Cma parallel evolution
% preMeta(iExp).ephysFN = 'Beto64chan-14102019-005'; 
% preMeta(iExp).expControlFN = '191014_Beto_generate_parallel(3)'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191014a\backup_10_14_2019_12_51_15' ;
% preMeta(iExp).comments = 'CMA parallel Evolution of Ch 26 and Ch 26 191014(CRP BXW). generate PCs for later PC space tuning';

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