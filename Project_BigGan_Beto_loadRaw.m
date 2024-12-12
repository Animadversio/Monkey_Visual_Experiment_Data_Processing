function B = Project_BigGan_Beto_loadRaw(varargin)

if nargin
    B = varargin{1} ; 
    meta = B.meta ;
    rasters = B.rasters;
    if isfield(B,'lfps')
        lfps = B.lfps;
    end
    Trials = B.Trials;
end

iExp =  1; 
preMeta(iExp).ephysFN  =     'Beto-20072020-002';
preMeta(iExp).expControlFN = '200720_Beto_generate_BigGAN(1)';
preMeta(iExp).stimuli = 'n:\Stimuli\2020-Evolutions\2020-07-20-Beto-01\2020-07-20-12-40-51' ;
preMeta(iExp).comments = 'CH 20 (-0.5 1) 3 1 hash  CMAES , BigGAN_class, truncnorm.random([1,128]) * 0.7 ' ;

iExp = 2 ; 
preMeta(iExp).ephysFN  =     'Beto-20072020-003';
preMeta(iExp).expControlFN = '200720_Beto_generate_integrated';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-Evolutions\2020-07-20-Beto-02\2020-07-20-12-58-55' ) ;
preMeta(iExp).comments = 'CH 20 (-0.5, 1) 3 1 CMAES, hash, CMAES' ;

iExp = 3 ; 
preMeta(iExp).ephysFN  =     'Beto-22072020-002';
preMeta(iExp).expControlFN = '200722_Beto_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-07-22-Beto-01\2020-07-22-10-14-22' ) ;
preMeta(iExp).comments = '5 (-0.7, 1) 3 2 CMAES, fc6 and biggan' ;

iExp = 4 ; 
preMeta(iExp).ephysFN  =     'Beto-23072020-002';
preMeta(iExp).expControlFN = '200723_Beto_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-07-23-Beto-01\2020-07-23-15-59-38' ) ;
preMeta(iExp).comments = '29 (-0.8 1) 3 2 CMAES with both ' ;

iExp = 5 ; 
preMeta(iExp).ephysFN  =     'Beto-24072020-002';
preMeta(iExp).expControlFN = '200724_Beto_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-07-24-Beto-01\2020-07-24-13-33-25' ) ;
preMeta(iExp).comments = '17 (-1,0) 3 1 CMAES 5/5 hash, fc6 and biggan' ;

iExp = iExp + 1 ;  
preMeta(iExp).ephysFN  =     'Beto-27072020-002';
preMeta(iExp).expControlFN = '200727_Beto_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-07-27-Beto-01\2020-07-27-14-29-47' ) ;
preMeta(iExp).comments = '26 (-1.1,0) 3 1 CMAES with fc6 (first) and BG (second),5/5 hash' ;

iExp = iExp + 1 ;  
preMeta(iExp).ephysFN  =     'Beto-28072020-003';
preMeta(iExp).expControlFN = '200728_Beto_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-07-28-Beto-01\2020-07-28-14-54-46' ) ;
preMeta(iExp).comments = 'generate 15 (-0.5,1) 3 1  using CMAES for fc6 (first thread) and CMAES-ZOHA' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Beto-29072020-002';
preMeta(iExp).expControlFN = '200729_Beto_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-07-29-Beto-01\2020-07-29-10-28-46' ) ;
preMeta(iExp).comments = 'Big gan generate start at 1028.↵CH 28 (0, 0.5) 3 1 CMAES fc6 and bigGan class. No climbing from either at 12 blocks' ;




iExp = iExp + 1 ;  
preMeta(iExp).ephysFN  =     'Beto-29072020-003';
preMeta(iExp).expControlFN = '200729_Beto_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-07-29-Beto-02\2020-07-29-11-00-02' ) ;
preMeta(iExp).comments = '18 (-0.5, 1) 4 1  CMAES fc6 vs biggan class. Hash. ' ;

iExp = iExp + 1 ;  
preMeta(iExp).ephysFN  =     'Beto-31072020-003';
preMeta(iExp).expControlFN = '200731_Beto_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-07-31-Beto-01\2020-07-31-11-07-26' ) ;
preMeta(iExp).comments = '29 [0 1.5]  3  2  CMAES , FC6,  BigGAN class. ' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Beto-03082020-002';
preMeta(iExp).expControlFN = '200803_Beto_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-08-03-Beto-01\2020-08-03-09-42-05' ) ;
preMeta(iExp).comments = 'Generate bigGan start at 942. CH 19 hash.  ↵19 (0, 3) 4 1 CMAES fc6 so ↵19  （0， 3）' ;


% iExp = iExp + 1; 
% preMeta(iExp).ephysFN  =     '';
% preMeta(iExp).expControlFN = '';
% preMeta(iExp).stimuli = fullfile('n:','' ) ;
% preMeta(iExp).comments = '' ;


iExp = iExp + 1 ;  
preMeta(iExp).ephysFN  =     'Beto-05082020-002';
preMeta(iExp).expControlFN = '200805_Beto_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-08-05-Beto-01\2020-08-05-13-56-11' ) ;
preMeta(iExp).comments = '28 (-0.5,0) 3 1 , 5/5 Hash' ;

iExp = iExp + 1 ;  
preMeta(iExp).ephysFN  =     'Beto-05082020-004';
preMeta(iExp).expControlFN = '200805_Beto_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-08-05-Beto-02\2020-08-05-15-09-06' ) ;
preMeta(iExp).comments = '28 (-0.5,0) 3 1 , 5/5 Hash   FC6 space CMAES-Hessian, BigGAN class CMAES-Hessian  ' ;

iExp = iExp + 1 ;  
preMeta(iExp).ephysFN  =     'Beto-07082020-002';
preMeta(iExp).expControlFN = '200807_Beto_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-08-07-Beto-01\2020-08-07-09-30-49' ) ;
preMeta(iExp).comments = 'CH 20 (-0.4, 1.4) 3 1  cmaes' ;

iExp = iExp + 1 ;  
preMeta(iExp).ephysFN  =     'Beto-11082020-003';
preMeta(iExp).expControlFN = '200811_Beto_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-08-11-Beto-02\2020-08-11-10-22-40' ) ;
preMeta(iExp).comments = 'CH 19 (0, 2) 4 1 hash' ;

iExp = iExp + 1 ;  
preMeta(iExp).ephysFN  =     'Beto-13082020-003';
preMeta(iExp).expControlFN = '200813_Beto_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-08-13-Beto-01\2020-08-13-10-11-18' ) ;
preMeta(iExp).comments = 'CH 28 (-1.5, 0) 3 1 CMAES hash.' ;

iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-17082020-002';
preMeta(iExp).expControlFN = '200817_Beto_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-08-17-Beto-01\2020-08-17-09-51-16' ) ;
preMeta(iExp).comments = 'CH 26 (-2, 0) 4 1 CMAES. Hash. Fc6 vs bigGAN class.' ;

iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-18082020-002';
preMeta(iExp).expControlFN = '200818_Beto_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-08-18-Beto-01\2020-08-18-16-02-08' ) ;
preMeta(iExp).comments = '25 (-0.1,0) 3 1 CMAES, fc6 and biggan class' ;

iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-20082020-002';
preMeta(iExp).expControlFN = '200820_Beto_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-08-20-Beto-01\2020-08-20-09-38-38' ) ;
preMeta(iExp).comments = 'CH 15 (0, 1.5) hash. 4 1 CMAES  fc6 vs bigGAN.' ;

iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-20082020-003';
preMeta(iExp).expControlFN = '200820_Beto_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-08-20-Beto-02\2020-08-20-10-11-12' ) ;
preMeta(iExp).comments = 'CH 49 (0, 0) 3 1 CMAES. Hash. Fc6 vs bigGAN. ' ;

iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-24082020-002';
preMeta(iExp).expControlFN = '200824_Beto_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-08-24-Beto-01\2020-08-24-09-44-23' ) ;
preMeta(iExp).comments = 'CH 52 (0, 0) 3 1 CMAES fc6 vs BigGAN.' ;

iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-24082020-004';
preMeta(iExp).expControlFN = '200824_Beto_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-08-24-Beto-02\2020-08-24-10-37-47' ) ;
preMeta(iExp).comments = 'CH 62 (0, 0) 3 1 CMAES fc6 vs BigGAN' ;

iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-26082020-005';
preMeta(iExp).expControlFN = '200826_Beto_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-08-26-Beto-01\2020-08-26-11-04-04' ) ;
preMeta(iExp).comments = 'CH 61 (-3, -2.5) 3 1 CMAES hash 1/5' ;

iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-28082020-002';
preMeta(iExp).expControlFN = '200828_Beto_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-08-28-Beto-01\2020-08-28-09-52-32' ) ;
preMeta(iExp).comments = 'CH 1 (0, 0) 3 1 Cmaes fc6 vs bigGan. 5/5 hash' ;

iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-31082020-002';
preMeta(iExp).expControlFN = '200831_Beto_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-08-31-Beto-01\2020-08-31-10-52-29' ) ;
preMeta(iExp).comments = '13 (-1,0) 3 1, CMAES, hash' ;

iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-01092020-003';
preMeta(iExp).expControlFN = '200901_Beto_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-01-Beto-01\2020-09-01-14-40-41' ) ;
preMeta(iExp).comments = '27 (0,0) 5 1 CMAES , SU 1/5, ' ;

iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-01092020-004';
preMeta(iExp).expControlFN = '200901_Beto_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-01-Beto-02\2020-09-01-15-16-24' ) ;
preMeta(iExp).comments = '29 (0,0) 5 1 CMAES, hash 5/5' ;

iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-03092020-002';
preMeta(iExp).expControlFN = '200903_Beto_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-03-Beto-01\2020-09-03-10-17-39' ) ;
preMeta(iExp).comments = 'Generate BigGAN start at 1017. CH 4 hash (-1, 0.5) 3 1 cmaes fc6 vs bigGAN.' ;

iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-03092020-003';
preMeta(iExp).expControlFN = '200903_Beto_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-03-Beto-02\2020-09-03-10-49-00' ) ;
preMeta(iExp).comments = ' CH 60  hash (0, 0) 3 1 cmaes fc6 vs bigGAN' ;

iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-07092020-002';
preMeta(iExp).expControlFN = '200907_Beto_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-07-Beto-01\2020-09-07-10-50-26' ) ;
preMeta(iExp).comments = '5 (0,0) 3 2  CMAES , fc6 and biggan, hash 5/5' ;

% Isolation problems
%    {'29'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             }
%     {'Beto-07092020-006'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              }
%     {'200907_Beto_generate_BigGAN(2)'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 }
%     {'N:\Stimuli\2020-BigGAN\2020-09-07-Beto-02\2020-09-07-11-35-06'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  }
%     {'006 at 1135 AM ↵51 (0,0) 2 1 CMAES hash 5/5 ↵fc6 and biggan. PSTH is okaish ↵upped reward to 180 and then 200 upon 1138 bmp ↵1141 AM, so one hour so far. He should go at least another 15-20 minutes. ↵bmp ↵Oh no, noise increasing across the board! Isolation in danger... ↵Bad noise gone again, and even lower. He's moving a lot, and  I think the connector may have problems. It has become harder to plug in ↵Upped to 400 ms ↵1156 AM, has been on for one hour and 15 minutes, which is not bad. Let's finish him at 1 20' ↵12 PM he's bumpy. Should not stop during a bump or he'll learn the wrong lesson. '}
% 

iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-10092020-002';
preMeta(iExp).expControlFN = '200910_Beto_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-10-Beto-01\2020-09-10-10-37-41' ) ;
preMeta(iExp).comments = 'CH 17 (-1.5, 0) 3 1 CMAES. Hash 5/5. ' ;

iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-10092020-004';
preMeta(iExp).expControlFN = '200910_Beto_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-10-Beto-02\2020-09-10-11-29-53' ) ;
preMeta(iExp).comments = 'CH 4 (0, 1.5) hash.' ;

% Isolation problems
% {'32'                                                                                                                                                                                                       }
%     {'Beto-14092020-002'                                                                                                                                                                                        }
%     {'200914_Beto_generate_BigGAN(2)'                                                                                                                                                                           }
%     {'N:\Stimuli\2020-BigGAN\2020-09-14-Beto-01\2020-09-14-10-12-27'                                                                                                                                            }
%     {'CH 51(0, 0) 3 1 cmaes. Hash 5/5. Start at 10:13. Fc6 vs BG. Both climbing at block 6. block 15, doing the same thing again! Hes not moving or bumping but noise is taken over. 22 blocks, stop at 1030.  '}


iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-14092020-003';
preMeta(iExp).expControlFN = '200914_Beto_generate_BigGAN(3)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-14-Beto-02\2020-09-14-10-33-12' ) ;
preMeta(iExp).comments = '59, (0, 0)  3 1 CMAES. Fc6 vs BG. Hash 5/5' ;

iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-17092020-002';
preMeta(iExp).expControlFN = '200917_Beto_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-17-Beto-01\2020-09-17-09-59-59' ) ;
preMeta(iExp).comments = 'Generate BG fc6 vs BG start at 959. CH 29 (0,1)  3 1 CMAES. SU 2/5' ;

iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-17092020-004';
preMeta(iExp).expControlFN = '200917_Beto_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-17-Beto-02\2020-09-17-10-57-23' ) ;
preMeta(iExp).comments = 'CH 64 (0, 0) 3 1 CMAES. Hash 5/5.' ;

iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-23092020-002';
preMeta(iExp).expControlFN = '200923_Beto_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-23-Beto-01\2020-09-23-09-55-07' ) ;
preMeta(iExp).comments = 'BG CH 57 (0,0) 3 1 cmaes. Fc6 vs BG.' ;

iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-23092020-003';
preMeta(iExp).expControlFN = '200923_Beto_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-23-Beto-02\2020-09-23-10-17-04' ) ;
preMeta(iExp).comments = 'CH 27 (0, 0) 3 1 hash. This ch has a SU (labeled as 2) but it is sparse and seems suppressed' ;

iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-23092020-004';
preMeta(iExp).expControlFN = '200923_Beto_generate_BigGAN(3)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-23-Beto-03\2020-09-23-10-34-45' ) ;
preMeta(iExp).comments = 'CH 54 (0, 0) 3 1 cmaes. Napping. Hash 5/5' ;

iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-23092020-005';
preMeta(iExp).expControlFN = '200923_Beto_generate_BigGAN(4)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-23-Beto-04\2020-09-23-10-59-09' ) ;
preMeta(iExp).comments = 'CH 55 (0,0) 3 1 hash 5/5.' ;

iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-29092020-005';
preMeta(iExp).expControlFN = '200929_Beto_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-29-Beto-02\2020-09-29-10-38-39' ) ;
preMeta(iExp).comments = 'BG vs fc6 CH 55 hash (0,0) 3 1 cmaes. Start at 1038.' ;


iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-27102020-002';
preMeta(iExp).expControlFN = '201027_Beto_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-10-27-Beto-01\2020-10-27-10-12-38' ) ;
preMeta(iExp).comments = 'CH 5 (0,0) 3 2 hash. cmaes. Fc6 vs BG class.' ;


iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-27102020-004';
preMeta(iExp).expControlFN = '201027_Beto_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-10-27-Beto-03\2020-10-27-10-53-12' ) ;
preMeta(iExp).comments = 'CH 42, same as above, start at 1053. B10, fc6 is climbing, maybe....' ;

iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-10112020-002';
preMeta(iExp).expControlFN = '201110_Beto_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-11-10-Beto-01\2020-11-10-13-33-32' ) ;
preMeta(iExp).comments = 'Lets try CH 39 for fun. (0,0) 3 2 hash. It climbed!' ;

iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-20112020-002';
preMeta(iExp).expControlFN = '201120_Beto_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-11-20-Beto-01\2020-11-20-12-08-07' ) ;
preMeta(iExp).comments = 'CH 29 (0, 1). 3 degrees. SU (2/5)' ;

iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-24112020-002';
preMeta(iExp).expControlFN = '201124_Beto_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-11-24-Beto-02\2020-11-24-10-48-51' ) ;
preMeta(iExp).comments = 'BG start at 1048. CH 48 (0,0) 3 1 cmaes. Hash 5/5. Fc6 vs BG.' ;


iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-30112020-002';
preMeta(iExp).expControlFN = '201130_Beto_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-11-30-Beto-01\2020-11-30-10-01-08' ) ;
preMeta(iExp).comments = 'Gen BG. Start at 1001. Ch 20 (0,2) 4 1 hash. I plugged in both arrays today. ' ;

iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-04122020-003';
preMeta(iExp).expControlFN = '201204_Beto_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-12-04-Beto-02\2020-12-04-10-12-46' ) ;
preMeta(iExp).comments = 'CH 57 (0,0) 2 1 hash. Start 1012' ;

iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-08122020-002';
preMeta(iExp).expControlFN = '201208_Beto_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-12-08-Beto-01\2020-12-08-08-21-42' ) ;
preMeta(iExp).comments = 'CH 36 (0,0) 2 1 hash gen BG. Start at 821. ' ;

% % this isjust an fc6 evolution -- should use for Project_Stability or
% % Validate
% iExp = iExp + 1 ; 
% preMeta(iExp).ephysFN  =     'Beto-13042021-002';
% preMeta(iExp).expControlFN = '210413_Beto_generate_BigGAN';
% preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2021-04-13-Beto-01\2021-04-13-09-48-15' ) ;
% preMeta(iExp).comments = 'BG gen start at 948. CH 19 SU. (0,0) 3 1. cmaes, fc6' ;

% % this isjust an fc6 evolution -- should use for Project_Stability or
% % Validate
% iExp = iExp + 1 ; 
% preMeta(iExp).ephysFN  =     'Beto-28052021-002';
% preMeta(iExp).expControlFN = '210528_Beto_generate_BigGAN(1)';
% preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2021-05-28-Beto-01\2021-05-28-10-28-01' ) ;
% preMeta(iExp).comments = 'CH 29 (0,0) 3 1 hash.' ;

iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-06052022-002';
preMeta(iExp).expControlFN = '220506_Beto_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2022-05-06-Beto-01\2022-05-06-13-36-52' ) ;
preMeta(iExp).comments = '47 [0,0 ] 3.5 unit s1' ;

iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-11052022-002';
preMeta(iExp).expControlFN = '220511_Beto_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2021-EvolDecomp\2022-05-11-Beto-01\2022-05-11-13-23-29' ) ;
preMeta(iExp).comments = '42 BigGAN FC6 [0 0 ] 4 MU 4/5' ;

iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-24052022-002';
preMeta(iExp).expControlFN = '220524_Beto_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2021-EvolDecomp\2022-05-24-Beto-01\2022-05-24-15-13-01' ) ;
preMeta(iExp).comments = 'chan 38 [0  0 ]  5 deg unit 1' ;

iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-31052022-002';
preMeta(iExp).expControlFN = '220531_Beto_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2021-EvolDecomp\2022-05-31-Beto-01\2022-05-31-13-26-33' ) ;
preMeta(iExp).comments = '002 BigGAN FC6 evol [1.3  -1.2]  1 4' ;

iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-07062022-003';
preMeta(iExp).expControlFN = '220607_Beto_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2021-EvolDecomp\2022-06-07-Beto-01\2022-06-07-13-29-46' ) ;
preMeta(iExp).comments = '64  MU5/5, [0  0]  5 deg rc , FC6   CMAES Hessian , BigGAN  CMAES' ;

iExp = iExp + 1 ; 
preMeta(iExp).ephysFN  =     'Beto-13062022-004';
preMeta(iExp).expControlFN = '220613_Beto_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2021-EvolDecomp\2022-06-13-Beto-02\2022-06-13-14-12-33' ) ;
preMeta(iExp).comments = 'Ch 6 [0.0, 0.0 ]  5  1 ' ;

% iExp = iExp + 1 ; 
% preMeta(iExp).ephysFN  =     '';
% preMeta(iExp).expControlFN = '';
% preMeta(iExp).stimuli = fullfile('n:','' ) ;
% preMeta(iExp).comments = '' ;

% iExp = iExp + 1 ; 
% preMeta(iExp).ephysFN  =     '';
% preMeta(iExp).expControlFN = '';
% preMeta(iExp).stimuli = fullfile('n:','' ) ;
% preMeta(iExp).comments = '' ;

% iExp = iExp + 1 ; 
% preMeta(iExp).ephysFN  =     '';
% preMeta(iExp).expControlFN = '';
% preMeta(iExp).stimuli = fullfile('n:','' ) ;
% preMeta(iExp).comments = '' ;

% iExp = iExp + 1 ; 
% preMeta(iExp).ephysFN  =     '';
% preMeta(iExp).expControlFN = '';
% preMeta(iExp).stimuli = fullfile('n:','' ) ;
% preMeta(iExp).comments = '' ;



% iExp = iExp + 1 ; 
% preMeta(iExp).ephysFN  =     '';
% preMeta(iExp).expControlFN = '';
% preMeta(iExp).stimuli = fullfile('n:','' ) ;
% preMeta(iExp).comments = '' ;

Project_General_copyMissingFiles(preMeta);

if exist('meta','var')
    fprintf('meta exists, so this must be an appending operation\n')
    iStart = length(meta)+1;
else
    iStart = 1 ; 
end

for iExp = iStart:length(preMeta)% 1:length(preMeta)
    tMeta = preMeta(iExp);
    try 
    [meta_,rasters_,lfps_,Trials_] = loadData(tMeta.ephysFN,'expControlFN',tMeta.expControlFN) ;
    meta_merged = rmfield( tMeta, intersect(fieldnames(tMeta), fieldnames(meta_)) );
    names = [fieldnames(meta_merged); fieldnames(meta_)];
    meta_ = cell2struct([struct2cell(meta_merged); struct2cell(meta_)], names, 1);
    
    meta{iExp} = meta_;
    rasters{iExp} = single(rasters_);
%     lfps{iExp} = single(lfps_);
    lfps{iExp} = nan;
    Trials{iExp} = Trials_;
    clear meta_  rasters_ lpfs_ Trials_ names meta_merged tMeta
    catch ME
    meta{iExp} = tMeta;
    rasters{iExp} = nan;
%     lfps{iExp} = single(lfps_);
    lfps{iExp} = nan;
    Trials{iExp} = nan;
    save(compose("S:\\%s-errorMess",tMeta.ephysFN),'ME')
    end
end

B.meta =        meta;
B.rasters =     rasters;
B.Trials =      Trials ;
B.lfps =        lfps ;