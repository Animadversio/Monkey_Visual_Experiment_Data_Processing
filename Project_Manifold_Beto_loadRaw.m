function [meta,rasters,lfps,Trials] = Project_Manifold_Beto_loadRaw(rowlist)
% iExp = 0;
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

% % day 011
% iExp = iExp + 1; %iExp = 11; % PC space exploration
% preMeta(iExp).ephysFN = 'Beto64chan-30102019-003'; 
% preMeta(iExp).expControlFN = '191030_Beto_selectivity_basic'; 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191030a-selectivity' ;
% preMeta(iExp).comments = ['Following Evolution of Ch26 Explore in the PC2-3 space only and add the pasupathy images to it. white background 4 rotations.'];
% 
% iExp = iExp + 1; %iExp = 11; % CMA
% preMeta(iExp).ephysFN = 'Beto64chan-30102019-001'; 
% preMeta(iExp).expControlFN = '191030_Beto_generate_parallel'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191030a\backup_10_30_2019_10_15_31' ;
% preMeta(iExp).comments = 'CMA Evolution of Ch 26 191030(CRP). generate PCs for later PC space tuning';
% % 
% % % % day 012
% iExp = iExp + 1; %iExp = 12; % PC space exploration
% preMeta(iExp).ephysFN = 'Beto64chan-31102019-004'; 
% preMeta(iExp).expControlFN = '191031_Beto_selectivity_basic'; 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191031a-selectivity' ;
% preMeta(iExp).comments = ['Derived from PCs of Ch 20 evolution 191031(CRP). 121 images from PC23 space, add the pasupathy images to it. white background 4 rotations.'];
% 
% iExp = iExp + 1;% iExp = 12; % CMA
% preMeta(iExp).ephysFN = 'Beto64chan-31102019-002'; 
% preMeta(iExp).expControlFN = '191031_Beto_generate_parallel'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191031a\backup_10_31_2019_13_24_49' ;
% preMeta(iExp).comments = 'CMA Evolution of Ch 20 191031(CRP). generate PCs for later PC space tuning';
% 
% % % day 013
% iExp = iExp + 1; %iExp = 13; % PC space exploration
% preMeta(iExp).ephysFN = 'Beto64chan-05112019-004'; 
% preMeta(iExp).expControlFN = '191105_Beto_selectivity_basic'; 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191105a-selectivity' ;
% preMeta(iExp).comments = ['Derived from PCs of Ch 40 evolution 191105(CRP). 121 images from PC23 space, add the pasupathy images to it. white background 4 rotations.'];
% 
% iExp = iExp + 1; %iExp = 13; % CMA
% preMeta(iExp).ephysFN = 'Beto64chan-05112019-003'; 
% preMeta(iExp).expControlFN = '191105_Beto_generate_parallel(1)'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191105a\backup_11_05_2019_11_23_09' ;
% preMeta(iExp).comments = 'CMA Evolution of Ch 40 evolution 191105(CRP). generate PCs for later PC space tuning';
% 
% % 
% % % day 014
% iExp = iExp + 1; %iExp = 14; % PC space exploration
% preMeta(iExp).ephysFN = 'Beto64chan-06112019-005'; 
% preMeta(iExp).expControlFN = '191106_Beto_selectivity_basic'; 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191106a-selectivity' ;
% preMeta(iExp).comments = ['Derived from PCs of Ch 42 evolution 191106(CRP). 121 images from PC23 space, add the pasupathy images to it. white background 4 rotations.'];
% 
% iExp = iExp + 1; %iExp = 14; % CMA
% preMeta(iExp).ephysFN = 'Beto64chan-06112019-003'; 
% preMeta(iExp).expControlFN = '191106_Beto_generate_parallel'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191106a\backup_11_06_2019_13_12_32' ;
% preMeta(iExp).comments = 'CMA Evolution of Ch 42 evolution 191106(CRP). generate PCs for later PC space tuning';
% 
% iExp = iExp + 1; %iExp = 15; % PC space exploration
% preMeta(iExp).ephysFN = 'Beto64chan-06112019-008'; 
% preMeta(iExp).expControlFN = '191106_Beto_selectivity_basic(1)'; 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191106b-selectivity' ;
% preMeta(iExp).comments = ['Derived from PCs of Ch 35 evolution 191106(CRP). 121 images from PC23 space, add the pasupathy images to it. white background 4 rotations.'];
% 
% iExp = iExp + 1; %iExp = 15; % CMA 
% preMeta(iExp).ephysFN = 'Beto64chan-06112019-007'; 
% preMeta(iExp).expControlFN = '191106_Beto_generate_parallel(1)'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191106b\backup_11_06_2019_14_17_35' ;
% preMeta(iExp).comments = 'CMA Evolution of Ch 35 evolution 191106(CRP). generate PCs for later PC space tuning';
% % day 015
% iExp = iExp + 1; % iExp = 16; % PC space exploration
% preMeta(iExp).ephysFN = 'Beto64chan-07112019-006'; 
% preMeta(iExp).expControlFN = '191107_Beto_selectivity_basic'; 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191107a-selectivity' ;
% preMeta(iExp).comments = ['Derived from PCs of Ch 36 evolution 191107(CRP). 121 images from PC23 space, add the pasupathy images to it. white background 4 rotations.'];
% 
% iExp = iExp + 1; % iExp = 16; % CMA
% preMeta(iExp).ephysFN = 'Beto64chan-07112019-004'; 
% preMeta(iExp).expControlFN = '191107_Beto_generate_parallel'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191107a\backup_11_07_2019_12_29_28' ;
% preMeta(iExp).comments = 'CMA Evolution of Ch 36 evolution 191107(CRP). generate PCs for later PC space tuning';
% 
% iExp = iExp + 1;% iExp = 17; % PC space exploration
% preMeta(iExp).ephysFN = 'Beto64chan-07112019-012'; 
% preMeta(iExp).expControlFN = '191107_Beto_selectivity_basic(1)'; %?????
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191107b-selectivity' ;
% preMeta(iExp).comments = ['Derived from PCs of Ch 34 evolution 191107(CRP). 121 images from PC23 space, add the pasupathy images to it. white background 4 rotations.'];
% 
% iExp = iExp + 1;% iExp = 17; % CMA
% preMeta(iExp).ephysFN = 'Beto64chan-07112019-010'; 
% preMeta(iExp).expControlFN = '191107_Beto_generate_parallel(2)'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191107b\backup_11_07_2019_13_50_18' ;
% preMeta(iExp).comments = 'CMA Evolution of Ch 34 evolution 191107(CRP). generate PCs for later PC space tuning';
% 
% % day 016
% iExp = iExp + 1;% iExp = 18; % PC space exploration
% preMeta(iExp).ephysFN = 'Beto64chan-08112019-004'; 
% preMeta(iExp).expControlFN = '191108_Beto_selectivity_basic'; %?????
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191108a-selectivity' ;
% preMeta(iExp).comments = ['Derived from PCs of Ch 37 evolution 191108(CRP). 121 images from PC23 space, add the pasupathy images to it. white background 4 rotations.'];
% 
% iExp = iExp + 1;% iExp = 18; % CMA
% preMeta(iExp).ephysFN = 'Beto64chan-08112019-003'; 
% preMeta(iExp).expControlFN = '191108_Beto_generate_parallel'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191108a\backup_11_08_2019_10_05_21' ;
% preMeta(iExp).comments = 'CMA Evolution of Ch 37 evolution 191108(CRP). generate PCs for later PC space tuning';

% % day 017
% iExp = iExp + 1;% iExp = 19; % PC space exploration
% preMeta(iExp).ephysFN = 'Beto64chan-11112019-006'; 
% preMeta(iExp).expControlFN = '191111_Beto_selectivity_basic(1)'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-19111a-selectivity' ;
% preMeta(iExp).comments = ['Derived from PCs of Ch 49 (V4) evolution 191111(BXW). 121 images from PC23 space, add the pasupathy images to it. white background 4 rotations.'];
% 
% iExp = iExp + 1;% iExp = 19; % CMA
% preMeta(iExp).ephysFN = 'Beto64chan-11112019-003'; 
% preMeta(iExp).expControlFN = '191111_Beto_generate_parallel'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-19111a\backup_11_11_2019_12_36_58' ;
% preMeta(iExp).comments = 'CMA Evolution of Ch 49 (V4) evolution 191111(BXW). generate PCs for later PC space tuning';
% 
% iExp = iExp + 1;% iExp = 20; % PC space exploration
% preMeta(iExp).ephysFN = 'Beto64chan-11112019-011'; 
% preMeta(iExp).expControlFN = '191111_Beto_selectivity_basic(2)'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-19111b-selectivity' ;
% preMeta(iExp).comments = ['Derived from PCs of Ch 60 (V4) evolution 191111(BXW). 121 images from PC23 space, add the pasupathy images to it. white background 4 rotations.'];
% 
% iExp = iExp + 1;% iExp = 20; % CMA
% preMeta(iExp).ephysFN = 'Beto64chan-11112019-010'; 
% preMeta(iExp).expControlFN = '191111_Beto_generate_parallel(1)'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-19111a\backup_11_11_2019_12_36_58' ;
% preMeta(iExp).comments = 'CMA Evolution of Ch 60 (V4) evolution 191111(BXW). generate PCs for later PC space tuning';

% % day 018
% iExp = iExp + 1;% iExp = 21; % PC space exploration
% preMeta(iExp).ephysFN = 'Beto64chan-12112019-004'; 
% preMeta(iExp).expControlFN = '191112_Beto_selectivity_basic'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-19112a-selectivity' ;
% preMeta(iExp).comments = ['Derived from PCs of Ch 56 (V4) evolution 191112(CRP). 121 images from PC23 space, add the pasupathy images to it. white background 4 rotations.'];
% 
% iExp = iExp + 1;% iExp = 21; % CMA
% preMeta(iExp).ephysFN = 'Beto64chan-12112019-002'; 
% preMeta(iExp).expControlFN = '191112_Beto_generate_parallel'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191112a\backup_11_12_2019_12_05_57' ;
% preMeta(iExp).comments = 'CMA Evolution of Ch 56 (V4) evolution 191112(CRP). generate PCs for later PC space tuning';

% % day 019
% iExp = iExp + 1;% iExp = 22; % PC space exploration
% preMeta(iExp).ephysFN = 'Beto64chan-13112019-006'; 
% preMeta(iExp).expControlFN = '191113_Beto_selectivity_basic'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191113a-selectivity' ;
% preMeta(iExp).comments = ['Derived from PCs of Ch 54 (V4) evolution 191113(BXW). [-4.8 -6.3] 3  121 images from PC23 space, add the pasupathy images to it. white background 4 rotations.'];
% 
% iExp = iExp + 1;% iExp = 22; % RF_mapper
% preMeta(iExp).ephysFN = 'Beto64chan-13112019-005'; 
% preMeta(iExp).expControlFN = '191113_Beto_rfMapper_basic(2)'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-06-RF-mapping\2019-11-13-Beto-chan54' ;
% preMeta(iExp).comments = 'RFMapping Heatmap of Ch 54 (V4) evolution 191113(BXW). chan 54 [-4.8 -6.3] -1.5:0.5:1.5 block043_thread000_gen_gen042_001719.bmp';

% iExp = iExp + 1;% iExp = 22; % CMA
% preMeta(iExp).ephysFN = 'Beto64chan-13112019-004'; 
% preMeta(iExp).expControlFN = '191113_Beto_generate_parallel(2)'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191113a\backup_11_13_2019_11_54_12' ;
% preMeta(iExp).comments = 'CMA Evolution of Ch 54 (V4) evolution 191113(BXW). [-4.8 -6.3] 3  generate PCs for later PC space tuning';
% 
% iExp = iExp + 1;% iExp = 23; % PC space exploration
% preMeta(iExp).ephysFN = 'Beto64chan-13112019-012'; 
% preMeta(iExp).expControlFN = '191113_Beto_selectivity_basic(1)'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191113b-selectivity' ;
% preMeta(iExp).comments = ['Derived from PCs of Ch 58 MU (V4) evolution 191113(BXW). [-3 -4] 3  121 images from PC23 space, add the pasupathy images to it. white background 4 rotations.'];

% iExp = iExp + 1;% iExp = 23; % RF_mapper
% preMeta(iExp).ephysFN = 'Beto64chan-13112019-011'; 
% preMeta(iExp).expControlFN = '191113_Beto_rfMapper_basic(6)'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-06-RF-mapping\2019-11-13-Beto-chan58' ;
% preMeta(iExp).comments = 'RFMapping Heatmap of Ch 58 (V4) evolution 191113(BXW). chan MU 58 [-3 -4] 3 deg -1.5:0.5:1.5 block043_thread000_gen_gen042_001719.bmp';

% iExp = iExp + 1;% iExp = 23; % CMA
% preMeta(iExp).ephysFN = 'Beto64chan-13112019-010'; 
% preMeta(iExp).expControlFN = '191113_Beto_generate_parallel(3)'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191113b\backup_11_13_2019_13_22_16' ;
% preMeta(iExp).comments = 'CMA Evolution of Ch 58 MU (V4) evolution 191113(BXW). [-3 -4] 3  generate PCs for later PC space tuning';

% % day 020
% iExp = iExp + 1;% iExp = 24; % PC space exploration
% preMeta(iExp).ephysFN = 'Beto64chan-14112019-007'; 
% preMeta(iExp).expControlFN = '191114_Beto_selectivity_basic'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191114a-selectivity' ;
% preMeta(iExp).comments = ['Derived from PCs of Ch 50 (V4) evolution 191114(CRP). ch 50 (-5,-6) 3 -degrees  121 images from PC23 space, add the pasupathy images to it. white background 4 rotations.'];

% iExp = iExp + 1;% iExp = 24; % RF_mapper
% preMeta(iExp).ephysFN = 'Beto64chan-14112019-006'; 
% preMeta(iExp).expControlFN = '191114_Beto_rfMapper_basic(4)'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-06-RF-mapping\2019-11-14-Beto' ;
% preMeta(iExp).comments = 'RFMapping Heatmap of Ch 50 (V4) evolution 191114(CRP). chan MU ch 50 (-5,-6) -1.5:0.5:1.5 block043_thread000_gen_gen042_001719.bmp';

% iExp = iExp + 1;% iExp = 24; % CMA
% preMeta(iExp).ephysFN = 'Beto64chan-14112019-004'; 
% preMeta(iExp).expControlFN = '191114_Beto_generate_parallel'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191114a\backup_11_14_2019_12_25_14' ;
% preMeta(iExp).comments = 'CMA Evolution of Ch 50 (V4) evolution 191114(CRP). ch 50 (-5,-6) 3 -degrees  generate PCs for later PC space tuning';

% iExp = iExp + 1;% iExp = 25; % PC space exploration
% preMeta(iExp).ephysFN = 'Beto64chan-14112019-007'; 
% preMeta(iExp).expControlFN = '191114_Beto_selectivity_basic'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191114a-selectivity' ;
% preMeta(iExp).comments = ['Derived from PCs of Ch 50 (V4) evolution 191114(CRP). ch 50 (-5,-6) 3 -degrees  121 images from PC23 space, add the pasupathy images to it. white background 4 rotations.'];
% 
% iExp = iExp + 1;% iExp = 25; % CMA
% preMeta(iExp).ephysFN = 'Beto64chan-14112019-009'; 
% preMeta(iExp).expControlFN = '191114_Beto_generate_parallel(1)'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191114b\backup_11_14_2019_13_50_18' ;
% preMeta(iExp).comments = 'CMA Evolution of Ch 57 (V4) evolution 191114(CRP). ch 57 (-3,-6) 3 1, SU 2/5 generate PCs for later PC space tuning';

% % day 021
% iExp = iExp + 1;% iExp = 25; % PC space exploration
% preMeta(iExp).ephysFN = 'Beto64chan-15112019-005'; 
% preMeta(iExp).expControlFN = '191115_Beto_selectivity_basic'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191115a-selectivity' ;
% preMeta(iExp).comments = ['Derived from PCs of Ch 57 (V4) evolution 191114(CRP). ch 50 (-3,-4) 4 -degrees  121 images from PC23 space, add the pasupathy images to it. white background 4 rotations.'];

% iExp = iExp + 1;% iExp = 25; % RF_mapper
% preMeta(iExp).ephysFN = 'Beto64chan-15112019-004'; 
% preMeta(iExp).expControlFN = '191115_Beto_rfMapper_basic(1)'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-06-RF-mapping\2019-11-15-Beto-chan57' ;
% preMeta(iExp).comments = 'RFMapping Heatmap of Ch 57 (V4) evolution 191118(BXW). ch  57 (-3,-4) 4-deg, grid -2:0.5:2';

% iExp = iExp + 1;% iExp = 25; % CMA
% preMeta(iExp).ephysFN = 'Beto64chan-15112019-003'; 
% preMeta(iExp).expControlFN = '191115_Beto_generate_parallel(1)'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191115a\backup_11_15_2019_12_56_11' ;
% preMeta(iExp).comments = 'CMA Evolution of Ch 57 (V4) evolution 191114(CRP). ch 57 (-3,-4) 4 1, MU 5/5 generate PCs for later PC space tuning';

% % day 022
% iExp = iExp + 1;% iExp = 26; % PC space exploration
% preMeta(iExp).ephysFN = 'Beto64chan-18112019-006'; 
% preMeta(iExp).expControlFN = '191118_Beto_selectivity_basic(1)'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191118a-selectivity' ;
% preMeta(iExp).comments = ['Derived from PCs of Ch 55 (V4) evolution 191118(BXW). ch  55 [-2.5 -2.8] 1 3 MU  -degrees  121 images from PC23 space, add the pasupathy images to it. white background 4 rotations.'];
% 
% iExp = iExp + 1;% iExp = 26; % RF_mapper
% preMeta(iExp).ephysFN = 'Beto64chan-18112019-004'; 
% preMeta(iExp).expControlFN = '191118_Beto_rfMapper_basic(1)'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-06-RF-mapping\2019-11-18-Beto-chan55' ;
% preMeta(iExp).comments = 'RFMapping Heatmap of Ch 55 (V4) evolution 191118(BXW). ch  55 [-2.5 -2.8] 1 3 MU -1.5:0.5:1.5';
% 
% iExp = iExp + 1;% iExp = 26; % CMA
% preMeta(iExp).ephysFN = 'Beto64chan-18112019-003'; 
% preMeta(iExp).expControlFN = '191118_Beto_generate_parallel(1)'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191118a\backup_11_18_2019_13_37_22' ;
% preMeta(iExp).comments = 'CMA Evolution of Ch 55 (V4) evolution 191118(BXW). ch  55 [-2.5 -2.8] 1 3 MU  generate PCs for later PC space tuning';

% % % day 023
% iExp = iExp + 1;% iExp = 27; % PC space exploration
% preMeta(iExp).ephysFN = 'Beto64chan-19112019-008'; 
% preMeta(iExp).expControlFN = '191119_Beto_selectivity_basic'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191119a-selectivity' ;
% preMeta(iExp).comments = ['Derived from PCs of Ch 26 (IT) evolution 191119(CRP). ch 26 (-2.1 -0.6) 3 1, SU 1/5  -degrees  121 images from PC23 space, add the pasupathy images to it. white background 4 rotations.'];
% 
% iExp = iExp + 1;% iExp = 27; % RF_mapper
% preMeta(iExp).ephysFN = 'Beto64chan-19112019-007'; 
% preMeta(iExp).expControlFN = '191119_Beto_rfMapper_basic(4)'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-06-RF-mapping\2019-11-19-Beto-chan26-heatmap' ;
% preMeta(iExp).comments = 'RFMapping Heatmap of Ch 26 (IT) evolution 191119(CRP). ch 26 (-2.1 -0.6) 3 1, SU 1/5 -1.5:0.5:1.5';
% 
% iExp = iExp + 1;% iExp = 27; % CMA
% preMeta(iExp).ephysFN = 'Beto64chan-19112019-006'; 
% preMeta(iExp).expControlFN = '191119_Beto_generate_parallel(1)'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191119a\backup_11_19_2019_12_11_50' ;
% preMeta(iExp).comments = 'CMA Evolution of Ch 26 (IT) evolution 191119(CRP). ch 26 (-2.1 -0.6) 3 1, SU 1/5  generate PCs for later PC space tuning';
% 
% % % day 024
% iExp = iExp + 1;% iExp = 28; % PC space exploration
% preMeta(iExp).ephysFN = 'Beto64chan-20112019-007'; 
% preMeta(iExp).expControlFN = '191120_Beto_selectivity_basic(1)'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191120a-seletivity' ;
% preMeta(iExp).comments = ['Derived from PCs of Ch 39 (V1) evolution 191120(CRP BXW). ch 39, MU (0,0) 3 1 4/5  -degrees  121 images from PC23 space, add the pasupathy images to it. white background 4 rotations.'];

% iExp = iExp + 1;% iExp = 28; % RF_mapper
% preMeta(iExp).ephysFN = 'Beto64chan-20112019-004'; 
% preMeta(iExp).expControlFN = '191120_Beto_rfMapper_basic(1)'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-06-RF-mapping\2019-11-20a-Beto-chan39' ;
% preMeta(iExp).comments = 'RFMapping Heatmap of Ch 39 (V1) evolution 191120(CRP BXW). ch 39, MU (0,0) 3 1 4/5 -1.5:0.5:1.5';
% 
% iExp = iExp + 1;% iExp = 28; % CMA
% preMeta(iExp).ephysFN = 'Beto64chan-20112019-002'; 
% preMeta(iExp).expControlFN = '191120_Beto_generate_parallel'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191120a\backup_11_20_2019_13_31_48' ;
% preMeta(iExp).comments = 'CMA Evolution of Ch 39 (V1) evolution 191120(CRP BXW). ch 39, MU (0,0) 3 1 4/5  generate PCs for later PC space tuning';
% 
% iExp = iExp + 1;% iExp = 29; % PC space exploration
% preMeta(iExp).ephysFN = 'Beto64chan-20112019-006'; 
% preMeta(iExp).expControlFN = '191120_Beto_selectivity_basic'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191120b-seletivity' ;
% preMeta(iExp).comments = ['Derived from PCs of Ch 39 (V1) evolution 191120(CRP BXW). ch 39, MU (0,0) 1 1 4/5  -degrees  121 images from PC23 space, add the pasupathy images to it. white background 4 rotations.'];
% 
% iExp = iExp + 1;% iExp = 29; % RF_mapper
% preMeta(iExp).ephysFN = 'Beto64chan-20112019-005'; 
% preMeta(iExp).expControlFN = '191120_Beto_rfMapper_basic(2)'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-06-RF-mapping\2019-11-20a-Beto-chan39' ;
% preMeta(iExp).comments = 'RFMapping Heatmap of Ch 39 (V1) evolution 191120(CRP BXW). ch 39, MU (0,0) 1 1 4/5 -1.5:0.5:1.5';
% 
% iExp = iExp + 1;% iExp = 29; % CMA
% preMeta(iExp).ephysFN = 'Beto64chan-20112019-003'; 
% preMeta(iExp).expControlFN = '191120_Beto_generate_parallel(1)'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191120b\backup_11_20_2019_14_21_08' ;
% preMeta(iExp).comments = 'CMA Evolution of Ch 39 (V1) evolution 191120(CRP BXW). ch 39, MU (0,0) 1 1 4/5  generate PCs for later PC space tuning';

% % day 025
% iExp = iExp + 1;% iExp = 30; % PC space exploration
% preMeta(iExp).ephysFN = 'Beto64chan-21112019-008'; 
% preMeta(iExp).expControlFN = '191121_Beto_selectivity_basic(2)'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191121b-selectivity' ;
% preMeta(iExp).comments = ['Derived from PCs of Ch 45 (V1) evolution 191121(CRP). ch 45 (0,0) 3 1, hash (MU 5/5)  -degrees  121 images from PC23 space, add the pasupathy images to it. white background 4 rotations.'];

% iExp = iExp + 1;% iExp = 30; % RF_mapper
% preMeta(iExp).ephysFN = 'Beto64chan-21112019-005'; 
% preMeta(iExp).expControlFN = '191121_Beto_rfMapper_basic(1)'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-06-RF-mapping\2019-11-21a-Beto-chan45' ;
% preMeta(iExp).comments = 'RFMapping Heatmap of Ch 45 (V1) evolution 191121(CRP). ch 45 (0,0) 3 1, hash (MU 5/5) -1.5:0.5:1.5';

% iExp = iExp + 1;% iExp = 30; % CMA
% preMeta(iExp).ephysFN = 'Beto64chan-21112019-004'; 
% preMeta(iExp).expControlFN = '191121_Beto_generate_parallel(2)'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191121b\backup_11_21_2019_13_02_02' ;
% preMeta(iExp).comments = 'CMA Evolution of Ch 45 (V1) evolution 191121(CRP). ch 45 (0,0) 3 1, hash (MU 5/5)  generate PCs for later PC space tuning';

% iExp = iExp + 1;% iExp = 31; % PC space exploration
% preMeta(iExp).ephysFN = 'Beto64chan-21112019-006'; 
% preMeta(iExp).expControlFN = '191121_Beto_selectivity_basic'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191121a-selectivity' ;
% preMeta(iExp).comments = ['Derived from PCs of Ch 45 (V1) evolution 191121(CRP). ch 45 (0,0) 1 1, hash (MU 5/5)  -degrees  121 images from PC23 space, add the pasupathy images to it. white background 4 rotations.'];

% iExp = iExp + 1;% iExp = 31; % CMA
% preMeta(iExp).ephysFN = 'Beto64chan-21112019-003'; 
% preMeta(iExp).expControlFN = '191121_Beto_generate_parallel(1)'; % 
% preMeta(iExp).stimuli = '\\storage1.ris.wustl.edu\crponce\Active\Stimuli\2019-Manifold\beto-191121a\backup_11_21_2019_12_42_39' ;
% preMeta(iExp).comments = 'CMA Evolution of Ch 45 (V1) evolution 191121(CRP). ch 45 (0,0) 1 1, hash (MU 5/5)  generate PCs for later PC space tuning';


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
ExpSpecTable_Aug = readtable("S:\ExpSpecTable_Augment.xls");
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
