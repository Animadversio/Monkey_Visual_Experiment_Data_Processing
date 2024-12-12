
function A = Project_BigGan_Alfa_loadRaw(varargin)

if nargin
    A = varargin{1} ; 
    meta = A.meta ;
    rasters = A.rasters;
    if isfield(A,'lfps')
        lfps = A.lfps;
    end
    Trials = A.Trials;
end

iExp =  1; 
preMeta(iExp).ephysFN  =     'Alfa-21072020-002';
preMeta(iExp).expControlFN = '200721_Alfa_generate_integrated(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-Evolutions\2020-07-21-Alfa-01\2020-07-21-09-51-23');
preMeta(iExp).comments = 'CH 12 (0, -0.5) 3 1 CMAES. Hash, FC6 generator' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-21072020-003';
preMeta(iExp).expControlFN = '200721_Alfa_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-07-21-Alfa-01\2020-07-21-10-06-15' ) ;
preMeta(iExp).comments = 'CH 12 (0, -0.5) 3 1 CMAES. Hash. BigGan' ;

iExp = iExp + 1;
preMeta(iExp).ephysFN  =     'Alfa-21072020-004';
preMeta(iExp).expControlFN = '200721_Alfa_generate_BigGAN(3)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-07-21-Alfa-02\2020-07-21-10-19-14' ) ;
preMeta(iExp).comments = 'CH 12 (0, -0.5) 3 1 CMAES. Hash., BigGan with monkey class' ;

iExp = iExp + 1;
preMeta(iExp).ephysFN  =     'Alfa-21072020-005';
preMeta(iExp).expControlFN = '200721_Alfa_generate_BigGAN(4)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-07-21-Alfa-03\2020-07-21-10-28-52' ) ;
preMeta(iExp).comments = 'CH 12 (0, -0.5) 3 1 CMAES. Hash., BigGan, intraclass again' ;

iExp = iExp + 1;
preMeta(iExp).ephysFN  =     'Alfa-21072020-007';
preMeta(iExp).expControlFN = '200721_Alfa_generate_integrated(3)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-Evolutions\2020-07-21-Alfa-02\2020-07-21-10-51-34' ) ;
preMeta(iExp).comments = 'Ch 1 (-4, -2.5) 4 1  SU 1/5. Very pretty unit. CMAES, FC6 generator' ;

iExp = iExp + 1;
preMeta(iExp).ephysFN  =     'Alfa-21072020-008';
preMeta(iExp).expControlFN = '200721_Alfa_generate_BigGAN(5)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-07-21-Alfa-04\2020-07-21-11-08-33' ) ;
preMeta(iExp).comments = 'Ch 1 (-4, -2.5) 4 1  SU 1/5. BigGan' ;


% Experiment that had to be stopped at 11 blocks because monkey was not
% working
% iExp = iExp + 1;
% preMeta(iExp).ephysFN  =     'Alfa-22072020-002';
% preMeta(iExp).expControlFN = '200722_Alfa_generate_BigGAN(3)';
% preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-07-22-Alfa-01\2020-07-22-10-30-01' ) ;
% preMeta(iExp).comments = 'Ch 6 (0, -3.5) 3 1 CMAES for fc6 and biggan class. Ugh hes not wanting to fixate. Added juice.. CH 6 SU 1/5. ' ;

% {'Alfa-22072020-003'                                                                                                                                                                                                                                                                                                                                                                                                           }
%     {'200722_Alfa_generate_BigGAN(4)'                                                                                                                                                                                                                                                                                                                                                                                              }
%     {'N:\Stimuli\2020-BigGAN\2020-07-22-Alfa-02\2020-07-22-11-05-32'                                                                                                                                                                                                                                                                                                                                                               }
%     {'Ch 16 (0,0) 3 1, CMAES and bigGan class . Start at 1105. I am trying one with rf in center and see if that helps his work ethic. Stop at 1121. Complete, may not keep↵7/22/20. Alfa did not work great today. Not fixating, did not want to work. He got 200 mL total over 1 hour. After return to cage, we discovered he pooped in his chair.... could be the culprit↵7/23/20, chlorhexidine. In rig at 1015. 3, 5, 6, 7,'}


iExp = iExp + 1;
preMeta(iExp).ephysFN  =     'Alfa-23072020-002';
preMeta(iExp).expControlFN = '200723_Alfa_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-07-23-Alfa-01\2020-07-23-10-33-01' ) ;
preMeta(iExp).comments = 'CH 3 (0, -3) 3 1 CMAES. Fc6 and BigGan class. CH 3 MU 1/5' ;

iExp = iExp + 1;
preMeta(iExp).ephysFN  =     'Alfa-23072020-003';
preMeta(iExp).expControlFN = '200723_Alfa_generate_BigGAN(3)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-07-23-Alfa-02\2020-07-23-11-06-08' ) ;
preMeta(iExp).comments = 'CH 6 (-0.3, -3.3) 3 1 CMAES fc6 and bigGan' ;

iExp = iExp + 1;
preMeta(iExp).ephysFN  =     'Alfa-23072020-005';
preMeta(iExp).expControlFN = '200723_Alfa_generate_BigGAN(4)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-07-23-Alfa-02b\2020-07-23-11-36-31' ) ;
preMeta(iExp).comments = 'CH 6 (0, -3) 3 1 CMAES fc6 and bigGan class. SU 2/5 this time.' ;

iExp = iExp + 1;
preMeta(iExp).ephysFN  =     'Alfa-24072020-002';
preMeta(iExp).expControlFN = '200724_Alfa_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-07-24-Alfa-01\2020-07-24-09-24-28' ) ;
preMeta(iExp).comments = 'CH 5 (-1.5, -3.8) 4 1 CMAES. SU 1/5. Start at 923' ;

iExp = iExp + 1;
preMeta(iExp).ephysFN  =     'Alfa-24072020-003';
preMeta(iExp).expControlFN = '200724_Alfa_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-07-24-Alfa-02\2020-07-24-10-06-35' ) ;
preMeta(iExp).comments = 'CH 5 (-1.5, -3.8) 4 1 CMAES. SU 1/5' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-27072020-002';
preMeta(iExp).expControlFN = '200727_Alfa_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-07-27-Alfa-01\2020-07-27-09-47-40' ) ;
preMeta(iExp).comments = 'CH 28 (-2, -2.5) 4 1 cmaes fc6 and bigGan class.' ;


iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-27072020-003';
preMeta(iExp).expControlFN = '200727_Alfa_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-07-27-Alfa-02\2020-07-27-10-18-57' ) ;
preMeta(iExp).comments = 'CH 2 (-2, -3.5) 3 1 cmaes. Fc6 and bigGan. Start at 10:18' ;


iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-28072020-002';
preMeta(iExp).expControlFN = '200728_Alfa_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-07-28-Alfa-01\2020-07-28-10-56-05' ) ;
preMeta(iExp).comments = 'CH 27 (-0.2, -0.2) 3 1 CMAES fc6 vs BigGan. SU 2/5.' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-28072020-004';
preMeta(iExp).expControlFN = '200728_Alfa_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-07-28-Alfa-02\2020-07-28-11-57-40' ) ;
preMeta(iExp).comments = 'Ch 9  (–2.2, -2.2) 4 1 CMAES. Hash,' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-29072020-003';
preMeta(iExp).expControlFN = '200729_Alfa_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-07-29-Alfa-01\2020-07-29-15-19-16' ) ;
preMeta(iExp).comments = '32  [-4 -4]  4  1  CMAES  fc6; 32  [-4 -4]  4  1  CMAES  BigGAN_class   0.7 * truncnorm' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-30072020-005';
preMeta(iExp).expControlFN = '200730_Alfa_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-07-30-Alfa-01\2020-07-30-10-51-34' ) ;
preMeta(iExp).comments = 'Ch 7  hash. (-0.5, -0.5) 3 1 CMAES fc6 and BigGan class.' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-04082020-002';
preMeta(iExp).expControlFN = '200804_Alfa_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-08-04-Alfa-01\2020-08-04-09-54-25' ) ;
preMeta(iExp).comments = 'CH 22 SU 2/5. (0, -2.5) 3 1 CMAES for both' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-04082020-004';
preMeta(iExp).expControlFN = '200804_Alfa_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-08-04-Alfa-02\2020-08-04-10-54-57' ) ;
preMeta(iExp).comments = 'CH 14 SU hash. (0, -2.5) 3 1 CMAES for both' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-06082020-002';
preMeta(iExp).expControlFN = '200806_Alfa_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-08-06-Alfa-01\2020-08-06-10-18-55' ) ;
preMeta(iExp).comments = 'CH 6 (-0.3, -3) 3 1  cmaes. Hash.' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-06082020-005';
preMeta(iExp).expControlFN = '200806_Alfa_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-08-06-Alfa-02\2020-08-06-11-41-51' ) ;
preMeta(iExp).comments = 'CH 15 (0, 0) 3 1  cmaes. Hash.' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-10082020-002';
preMeta(iExp).expControlFN = '200810_Alfa_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-08-10-Alfa-01\2020-08-10-09-32-31' ) ;
preMeta(iExp).comments = 'Gen big gan start at 932. ↵CH 6 (0, -2.6) 3 1 CMEAS SU 1/5, sounds great. Very VR' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-10082020-003';
preMeta(iExp).expControlFN = '200810_Alfa_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-08-10-Alfa-02\2020-08-10-09-59-48' ) ;
preMeta(iExp).comments = 'CH 19 (1, -3) 4 1 CMEAS hash' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-10082020-005';
preMeta(iExp).expControlFN = '200810_Alfa_generate_BigGAN(3)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-08-10-Alfa-03\2020-08-10-11-16-31' ) ;
preMeta(iExp).comments = 'Ch 6 evolving to hash (using whole ch, not separated). (-0.3, -3) 3 1 CMAES. Fc6 vs biggan' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-12082020-003';
preMeta(iExp).expControlFN = '200812_Alfa_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-08-12-Alfa-01\2020-08-12-09-51-49' ) ;
preMeta(iExp).comments = 'CH 29 (0, -2) 4 1 CMAES, fc6 vs BigGan.' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-13082020-002';
preMeta(iExp).expControlFN = '200813_Alfa_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-08-13-Alfa-01\2020-08-13-12-31-20' ) ;
preMeta(iExp).comments = '29 (0,-2) 4 1 CMAES, for fc6 and big gan' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-13082020-004';
preMeta(iExp).expControlFN = '200813_Alfa_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-08-13-Alfa-03\2020-08-13-13-11-01' ) ;
preMeta(iExp).comments = '004 56 (0,0) 3 1 CMAES, 5/5 hash' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-13082020-005';
preMeta(iExp).expControlFN = '200813_Alfa_generate_BigGAN(3)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-08-13-Alfa-04\2020-08-13-13-34-42' ) ;
preMeta(iExp).comments = '60 (0,0) 3 1 CMAES fc6 and biggan 5/5 hash' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-14082020-002';
preMeta(iExp).expControlFN = '200814_Alfa_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-08-14-Alfa-01\2020-08-14-09-54-17' ) ;
preMeta(iExp).comments = '51 (-1.8, -1.8) 4 1 CMAES fc6 vs BigGan. CH 51, SU 1/5' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-14082020-004';
preMeta(iExp).expControlFN = '200814_Alfa_generate_BigGAN(3)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-08-14-Alfa-02\2020-08-14-11-06-59' ) ;
preMeta(iExp).comments = '63 hash (0, -1) 3 1 CMAES fc6 vs biggan. ' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-19082020-002';
preMeta(iExp).expControlFN = '200819_Alfa_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-08-19-Alfa-01\2020-08-19-09-59-21' ) ;
preMeta(iExp).comments = 'CH 58 (0, -2) 4 1 CMAES. Hash. Fc6 vs bigGAN class' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-19082020-003';
preMeta(iExp).expControlFN = '200819_Alfa_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-08-19-Alfa-02\2020-08-19-10-25-54' ) ;
preMeta(iExp).comments = 'CH 55 (-1, -2.2) 4 1 CMAES fc6 vs bigGAN.' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-19082020-006';
preMeta(iExp).expControlFN = '200819_Alfa_generate_BigGAN(5)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-08-19-Alfa-03\2020-08-19-11-24-20' ) ;
preMeta(iExp).comments = 'CH 50 (0, -5.5) SU 1/5. Great SU, just very far rf. ' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-25082020-005';
preMeta(iExp).expControlFN = '200825_Alfa_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-08-25-Alfa-01\2020-08-25-10-20-44' ) ;
preMeta(iExp).comments = 'CH 60 (0, -0.5) 3 1 CMAES hash. Fc6 v bigGAN.' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-25082020-006';
preMeta(iExp).expControlFN = '200825_Alfa_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-08-25-Alfa-02\2020-08-25-10-46-06' ) ;
preMeta(iExp).comments = 'CH 59 SU 1/5. (-2.5, -2.5) 3 1 CMAES fc6 vs bigGAN' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-27082020-002';
preMeta(iExp).expControlFN = '200827_Alfa_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-08-27-Alfa-01\2020-08-27-09-41-28' ) ;
preMeta(iExp).comments = 'CH 7 (-1.5, -1.5) 3 1 CMAES hash 1/5. Fc6 vs bigGAN' ;

% thie experiment only has 13 blocks
% {'Alfa-27082020-004'                                                                                                                                                                      }
%     {'200827_Alfa_generate_BigGAN(2)'                                                                                                                                                         }
%     {'N:\Stimuli\2020-BigGAN\2020-08-27-Alfa-02\2020-08-27-10-50-27'                                                                                                                          }
%     {'52 (0,0) 3 1 CMAES. Hash. Fc6 vs bigGAN start at 1050Might be the last exp. he is slowing down. Been working for 1.5 hours. Blasting him. Stop at 1109. Only 13 blocks. May not keep.  '}


iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-01092020-002';
preMeta(iExp).expControlFN = '200901_Alfa_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-01-Alfa-01\2020-09-01-10-15-53' ) ;
preMeta(iExp).comments = 'Biggan start at 10:15. CH 3 hash (-1.5, -1.5) 3 1 CMAES fc6 vs BigGAN. ' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-01092020-003';
preMeta(iExp).expControlFN = '200901_Alfa_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-01-Alfa-02\2020-09-01-10-42-42' ) ;
preMeta(iExp).comments = 'CH 52 (0, 0) 3 1 CMAES . Hash.' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-02092020-002';
preMeta(iExp).expControlFN = '200902_Alfa_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-02-Alfa-01\2020-09-02-09-54-51' ) ;
preMeta(iExp).comments = ' CH 50  SU 1/5. (0, -5) 3 1 CMAES fc6 vs BigGAN' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-02092020-003';
preMeta(iExp).expControlFN = '200902_Alfa_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-02-Alfa-02\2020-09-02-10-17-21' ) ;
preMeta(iExp).comments = ' CH 64 (0, -1) 3 1 CMAES fc6 vs biggan.' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-02092020-005';
preMeta(iExp).expControlFN = '200902_Alfa_generate_BigGAN(3)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-02-Alfa-03\2020-09-02-10-59-57' ) ;
preMeta(iExp).comments = 'CH 49 (0, 0) 3 1 CMAES fc6 vs biggan' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-04092020-003';
preMeta(iExp).expControlFN = '200904_Alfa_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-04-Alfa-01\2020-09-04-12-06-08' ) ;
preMeta(iExp).comments = 'CH 25 (2.5, -2.5) 4 1 . SU 1/5. CMAES fc6 vs bigGAN start at 1206' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-08092020-002';
preMeta(iExp).expControlFN = '200908_Alfa_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-08-Alfa-01\2020-09-08-10-11-27' ) ;
preMeta(iExp).comments = 'CH 57 (0,0) 3 1 hash(1/5). Generate bigGAN fc6 vs BG' ;


iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-08092020-003';
preMeta(iExp).expControlFN = '200908_Alfa_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-08-Alfa-02\2020-09-08-10-31-42' ) ;
preMeta(iExp).comments = 'CH 21 (0, -2.3) 3 1 hash 1/5. Start at 1031. Fc6 vs BG.' ;


iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-08092020-005';
preMeta(iExp).expControlFN = '200908_Alfa_generate_BigGAN(3)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-08-Alfa-03\2020-09-08-11-10-09' ) ;
preMeta(iExp).comments = 'CH 62 (0, 0) 3 1 CMAES. Hash (1/5) start at 1110' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-11092020-002';
preMeta(iExp).expControlFN = '200911_Alfa_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-11-Alfa-01\2020-09-11-09-46-38' ) ;
preMeta(iExp).comments = 'CH 32 (-3.5, -3.5) 4 1 CMAES. MU 2/5' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-11092020-005';
preMeta(iExp).expControlFN = '200911_Alfa_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-11-Alfa-02\2020-09-11-10-45-27' ) ;
preMeta(iExp).comments = 'CH 30. (-2, -2, ) 3 1 Start at 1045. CMAES hash 1/5' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-16092020-002';
preMeta(iExp).expControlFN = '200916_Alfa_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-16-Alfa-01\2020-09-16-10-04-30' ) ;
preMeta(iExp).comments = 'Generate BG start at 1005. CH 10 (-0.2, -2.2) 3 1 CMAES. Hash 5/5. Fc6 vd=s BG' ;


iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-16092020-004';
preMeta(iExp).expControlFN = '200916_Alfa_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-16-Alfa-02\2020-09-16-10-53-04' ) ;
preMeta(iExp).comments = 'CH 15 (0,0) 3 1 CMAES. SU 2/5. Generate BG ' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-18092020-002';
preMeta(iExp).expControlFN = '200918_Alfa_generate_BigGAN(3)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-18-Alfa-01\2020-09-18-09-34-18' ) ;
preMeta(iExp).comments = 'CH 6 (0, -2.7) 3 1 CMAES. SU 1/5. Fc6 vs BG like normal' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-18092020-003';
preMeta(iExp).expControlFN = '200918_Alfa_generate_BigGAN(4)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-18-Alfa-02\2020-09-18-09-51-51' ) ;
preMeta(iExp).comments = 'CH 6 (0, -2.7) 5 1 CMAES. SU 1/5. Fc6 vs BG like normal. ' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-18092020-004';
preMeta(iExp).expControlFN = '200918_Alfa_generate_BigGAN(5)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-18-Alfa-03\2020-09-18-10-21-58' ) ;
preMeta(iExp).comments = '004 at 1021 AM ↵repeating above but this time with CMAES_hess' ;


iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-22092020-002';
preMeta(iExp).expControlFN = '200922_Alfa_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-22-Alfa-01\2020-09-22-09-52-55' ) ;
preMeta(iExp).comments = 'Gen BG start at 952. Fc6 vs BG. CH 11 (-2, -1) 3 1 cmaes. Hash 5/5.' ;

% one generation where baseline went nuts
% {'37'                                                                                                                                                                                                                                              }
%     {'Alfa-04092020-003'                                                                                                                                                                                                                               }
%     {'200904_Alfa_generate_BigGAN(2)'                                                                                                                                                                                                                  }
%     {'N:\Stimuli\2020-BigGAN\2020-09-04-Alfa-01\2020-09-04-12-06-08'                                                                                                                                                                                   }
%     {'CH 25 (2.5, -2.5) 4 1 . SU 1/5. CMAES fc6 vs bigGAN start at 1206. Block 12. Nothing yet. 25 blocks complete stop at 1234. Neither climbed. Fc6 image is interesting now that I look at it after the fact, may be worth looking at.  complete.  '}


iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-22092020-004';
preMeta(iExp).expControlFN = '200922_Alfa_generate_BigGAN(3)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-22-Alfa-03\2020-09-22-10-33-35' ) ;
preMeta(iExp).comments = 'CH 18 (0, -4) 41 cmaes. Hash 5/5. Fc6 vs BG. ' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-22092020-005';
preMeta(iExp).expControlFN = '200922_Alfa_generate_BigGAN(4)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-22-Alfa-04\2020-09-22-10-57-28' ) ;
preMeta(iExp).comments = 'CH 61 (0, 0) 3 1 cmaes hash 5/5. Fc6 vs BG' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-22092020-006';
preMeta(iExp).expControlFN = '200922_Alfa_generate_BigGAN(5)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-22-Alfa-05\2020-09-22-11-16-22' ) ;
preMeta(iExp).comments = ' CH 54(0, 0) 3 1. cmaes hash 5/5 fc vs BG' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-24092020-002';
preMeta(iExp).expControlFN = '200924_Alfa_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-24-Alfa-01\2020-09-24-09-28-54' ) ;
preMeta(iExp).comments = 'CH 24 (-1 -3.5) 4 1 cmaes fc6 vs BG. Hash 5/5. ' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-24092020-003';
preMeta(iExp).expControlFN = '200924_Alfa_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-24-Alfa-02\2020-09-24-09-46-34' ) ;
preMeta(iExp).comments = 'CH 4 (-2, -2.5) 4 1 cmaes. Fc6 vs BG.' ;


iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-24092020-004';
preMeta(iExp).expControlFN = '200924_Alfa_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-24-Alfa-03\2020-09-24-10-08-59' ) ;
preMeta(iExp).comments = 'CH 14 (0, -2) 3 1 cmaes fc6 vs bg' ;


iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-28092020-002';
preMeta(iExp).expControlFN = '200928_Alfa_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-28-Alfa-01\2020-09-28-11-33-00' ) ;
preMeta(iExp).comments = 'CH 34  hash(0,0) 2 1 cmaes fc6 vs BG' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-28092020-004';
preMeta(iExp).expControlFN = '200928_Alfa_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-28-Alfa-02\2020-09-28-11-52-42' ) ;
preMeta(iExp).comments = ' CH 48(0,0) 2 1 cmaes. Hash. Fc6 vs BG. ' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-28092020-006';
preMeta(iExp).expControlFN = '200928_Alfa_generate_BigGAN(3)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-28-Alfa-03\2020-09-28-12-12-31' ) ;
preMeta(iExp).comments = 'BG CH 38 (0,0) 2 1 hash. Cmaes fc6 vs BG' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-28092020-008';
preMeta(iExp).expControlFN = '200928_Alfa_generate_BigGAN(4)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-09-28-Alfa-04\2020-09-28-12-35-24' ) ;
preMeta(iExp).comments = 'Bg vs fc6 CH 45 (0,0) 2 1 cmaes. Hash.' ;


% iExp = iExp + 1; 
% preMeta(iExp).ephysFN  =     '';
% preMeta(iExp).expControlFN = '';
% preMeta(iExp).stimuli = fullfile('n:','' ) ;
% preMeta(iExp).comments = '' ;


iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-21102020-002';
preMeta(iExp).expControlFN = '201021_Alfa_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-10-21-Alfa-01\2020-10-21-10-43-19' ) ;
preMeta(iExp).comments = '002 Generate BigGAN started 10:43 ↵Chan 6 [0 –2.4] 4 1   FC6  CMA Hessian' ;


iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-27102020-003';
preMeta(iExp).expControlFN = '201027_Alfa_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-10-27-Alfa-01\2020-10-27-11-45-41' ) ;
preMeta(iExp).comments = '003   starts 11:45  ↵Chan 26  ↵FC6  CMAES  [-1.5 1.5] 3 1 ↵BigGAN  CMAES    [-1.5 1.5] 3 1 ↵SU 1.5 / 5 PSTH looks sparse and great!' ;


iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-27102020-004';
preMeta(iExp).expControlFN = '201027_Alfa_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-10-27-Alfa-02\2020-10-27-12-08-19' ) ;
preMeta(iExp).comments = 'Chan 9 ↵FC6  CMAES  [-1.3  -2.2] 4 1 ↵BigGAN  CMAES    [-1.3  -2.2] 4 1 ↵Loading up takes forever' ;


iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-29102020-002';
preMeta(iExp).expControlFN = '201029_Alfa_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-10-29-Alfa-01\2020-10-29-11-30-46' ) ;
preMeta(iExp).comments = '002 Generate BigGAN ↵So seems james'' code face some issue?  ↵27 [-0.5  0]  2  3  CMAES ↵27 [-0.5  0]  2  3  ' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-04112020-002';
preMeta(iExp).expControlFN = '201104_Alfa_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-11-04-Alfa-01\2020-11-04-11-50-56' ) ;
preMeta(iExp).comments = '↵28  [-2 –3.5]  5  1 CMAES FC6   ↵28  [-2 –3.5]  5  1 CMAES BigGAN' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-06112020-002';
preMeta(iExp).expControlFN = '201106_Alfa_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-11-06-Alfa-01\2020-11-06-11-29-50' ) ;
preMeta(iExp).comments = '002 G BG (1) ↵1 [-3 –3.5] 5 3 CMAES ↵1 [-3 –3.5] 5 3 CMAES' ;


iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-06112020-003';
preMeta(iExp).expControlFN = '201106_Alfa_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-11-06-Alfa-02\2020-11-06-11-53-09' ) ;
preMeta(iExp).comments = '19 [0 –2.4] 4 1 CMAES ↵19 [0 –2.4] 4 1 CMAES ↵BigGAN evolves a bit and it''s on top of other threads' ;


iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-11112020-003';
preMeta(iExp).expControlFN = '201111_Alfa_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-11-11-Alfa-01\2020-11-11-11-32-43' ) ;
preMeta(iExp).comments = '26 [-0.1 2.6 ] 4.5 2 FC6 , 26 [-0.1 2.6 ] 4.5 2 BigGAN' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-25112020-003';
preMeta(iExp).expControlFN = '201125_Alfa_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2020-11-25-Alfa-02\2020-11-25-10-21-04' ) ;
preMeta(iExp).comments = 'Ch 12 (0, -0.8) 3 1 hash 5/5. Fc6 vs BG with still images,' ;

% in lab notebook, these files appear in CMAES reduced dimension page
%  {'65'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   }
%     {'Alfa-06072020-002'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    }
%     {'200706_Alfa_generate_integrated(1)'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   }
%     {'N:\Stimuli\2020-BigGAN\2021-01-07-Alfa-01\2021-01-07-11-11-04'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        }
%     {'002 genrate BigGAN ↵Copied reference  ↵58  [-1 -2.4] 4 1 CMAES fc6  ↵58  [-1 -2.4] 4 1 CMAES BigGAN full  ↵Very fast experiment.  ↵Is there some common features for the units that are easy to evolve?  ↵Single units?  ↵This guy is definitely evolving~ Even the BigGAN thread.  ↵Wow…..The FC6 thread is evolving great!  ↵This unit is a treasure.  ↵FC6  ↵OK…. Waste 15 mins doing eye calibaration with no results and monkey is mad:( Sad  ↵He is not cooperating  ↵Starts from 5-6 -> 18-19 sp.s ↵FC6 thread super successful. BigGAN thread goes nowhere.  ↵45 mins finished Monkey has drank 350 +100ml.  ↵Yeah backup has ruined 10 mins…. Not good  ↵The think looks like a curved surface with shading on it.  ↵Completed '}

% {'67'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  }
%     {'Alfa-06072020-003'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   }
%     {'200706_Alfa_generate_integrated(2)'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  }
%     {'N:\Stimuli\2020-BigGAN\2021-01-07-Alfa-02\2021-01-07-12-11-50'                                                                                                                                                                                                                                                                                                                                                                                                                                                       }
%     {'003 ↵58  [-1 -2.4] 4 1 CMAES BigGAN class space ↵I feel BigGAN is not the style it likes ↵Stuck around 10 sp/s ↵Interestingly, the cell evolves a cat ball, which is like the surface in last exp.  ↵I feel one problem that BigGAN space have is that it’s not converging / attracting. ↵Still Stuck around 10 sp/s.   ↵450 mlNow ↵Cool interesting! This guy suddenly evolved at 30 blocks and stay relatively high.  ↵So this is actually successful.  ↵35 mins 37 blocks matched.  ↵500mL +  ↵Completed '}


iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-07012021-002';
preMeta(iExp).expControlFN = '210107_Alfa_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2021-01-07-Alfa-01\2021-01-07-11-11-04' ) ;
preMeta(iExp).comments = '58  [-1 -2.4] 4 1 CMAES fc6  ↵58  [-1 -2.4] 4 1 CMAES BigGAN full' ;


iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-12012021-002';
preMeta(iExp).expControlFN = '210112_Alfa_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2021-01-12-Alfa-01\2021-01-12-11-33-31' ) ;
preMeta(iExp).comments = '27 [0.1 0.1] 4 deg  unit 1 CMAES fc6 ↵27 [0.1 0.1] 4 deg  unit 1 CMAES BigGAN full' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-13012021-002';
preMeta(iExp).expControlFN = '210113_Alfa_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2021-01-13-Alfa-01\2021-01-13-11-40-24' ) ;
preMeta(iExp).comments = 'Chan 28 unit 1  [-2.3 -4.5]  5 deg FC6  ↵Chan 28 unit 1  [-2.3 -4.5]  5 deg BigGAN' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-20012021-003';
preMeta(iExp).expControlFN = '210120_Alfa_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2021-01-20-Alfa-01\2021-01-20-12-10-29' ) ;
preMeta(iExp).comments = 'Chan 32 ↵Use ref img from chan_32_exp_04' ;

% this one was aborted
% {'71'                                                                                                                                                                           }
%     {'Alfa-20012021-004'                                                                                                                                                            }
%     {'210120_Alfa_generate_BigGAN(2)'                                                                                                                                               }
%     {'N:\Stimuli\2020-BigGAN\2021-01-20-Alfa-02\2021-01-20-12-53-17'                                                                                                                }
%     {'004 Generate BigGAN ↵24  [-3.5  -5]  6  2  ↵24  [-3.5  -5]  6  2  ↵Seems this is evolving in the wrong way.....  ↵ANd also no visual response.  ↵Not Saving!!!! ↵Aborted '}


iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-20012021-005';
preMeta(iExp).expControlFN = '210120_Alfa_generate_BigGAN(3)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2021-01-20-Alfa-03\2021-01-20-13-06-42' ) ;
preMeta(iExp).comments = '55  [0  -1.5]  4  1  ↵55  [0  -1.5]  4  1' ;


iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-26012021-002';
preMeta(iExp).expControlFN = '210126_Alfa_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2021-01-26-Alfa-01\2021-01-26-12-03-50' ) ;
preMeta(iExp).comments = '58 [-0.6, -4.2] 5 deg 1  ↵FC6  ↵BigGAN full  ↵SU 2/5 ↵Pretty good, two threads are both quite successful and doing their job' ;


iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-26012021-003';
preMeta(iExp).expControlFN = '210126_Alfa_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2021-01-26-Alfa-02\2021-01-26-12-45-32' ) ;
preMeta(iExp).comments = '7 [-1 –2.5] 4 1 CMAES' ;

% monkey did not work on this one
%     {'75'                                                                                                                                                                                                                                                                                                                                                                                                           }
%     {'Alfa-28012021-002'                                                                                                                                                                                                                                                                                                                                                                                            }
%     {'210128_Alfa_generate_BigGAN'                                                                                                                                                                                                                                                                                                                                                                                  }
%     {'N:\Stimuli\2020-BigGAN\2021-01-28-Alfa-01\2021-01-28-10-38-03'                                                                                                                                                                                                                                                                                                                                                }
%     {'002 Generate BigGAN ↵24 [-4 –5]  6 deg  unit 2  ↵24 [-4 –5]  6 deg  unit 2  ↵MU 4.5/5 ↵Reference copied.  ↵Why the PSTH looks different from what we get from RF mapping exp.  ↵Visual response is bad. fixation is bad. Break fixation rate 98%!!!!  ↵20mins 9blocks evolve in the negative direction.....This won't work  ↵Terminates ↵01 folder  ↵Not completed only 10 blocks..... not successful '}
 
iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-28012021-003';
preMeta(iExp).expControlFN = '210128_Alfa_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2021-01-28-Alfa-02\2021-01-28-11-05-28' ) ;
preMeta(iExp).comments = '5[-2.5 –3.9]  4 deg  unit 1  ↵5[-2.5 –3.9]  4 deg  unit 1  ↵4.5/5 signal' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-18022021-003';
preMeta(iExp).expControlFN = '210218_Alfa_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2021-02-18-Alfa-01\2021-02-18-11-25-57' ) ;
preMeta(iExp).comments = '003 generate BigGAN↵Single unit! Good signal, weird PSTH, bad rf but strong selectivity↵The PSTH is pretty late' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-19022021-003';
preMeta(iExp).expControlFN = '210219_Alfa_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2021-02-19-Alfa-01\2021-02-19-13-03-12' ) ;
preMeta(iExp).comments = '003 evolving BigGAN FC6↵FC6 Hessian, BigGAN CMAES full space↵[-3 -3 ] 7 deg huge image' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-25022021-002';
preMeta(iExp).expControlFN = '210225_Alfa_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2021-02-25-Alfa-01\2021-02-25-11-23-18' ) ;
preMeta(iExp).comments = 'Generate BigGAN FC6↵3  [-1.5  -2.0]  6  1  CMAES_Hessian↵3  [-1.5  -2.0]  6  1  CMAES' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-26022021-003';
preMeta(iExp).expControlFN = '210226_Alfa_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2021-02-26-Alfa-01\2021-02-26-12-53-02' ) ;
preMeta(iExp).comments = '18 [0.5 -5] 1 7deg.↵BigGAN goes into face zone' ;

% this one had a huge threshold change.
% they tend to be very hard to analyze
% {'81'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         }
%     {'Alfa-12032021-002'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          }
%     {'210312_Alfa_generate_BigGAN'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                }
%     {'N:\Stimuli\2020-BigGAN\2021-03-12-Alfa-01\2021-03-12-11-39-05'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              }
%     {'002 BigGAN FC6 evolution  ↵Copy ref 8 ↵8 [-5 –2.5] 6 1 'CMAESHess'  fc6 ↵8 [-5 –2.5] 6 1 'CMAES' BigGAN ↵Promising sparse PSTH ↵Su 1.5/5 good signal.  ↵Caveat, in this exp at first the unit drifted a bit so I moved the sorting criterion to make it less sparse.  ↵This cell prefer FC6 > BigGAN to a large margin! FC6 ~ 45 BigGAN ~ 15 ↵Evolution for the FC6 is not stable... ↵This guy is doing up down up down  ↵May exclude this cell in analysius... ↵Seems not really successful.  ↵003 BigGAN FC6 evolution  ↵Copy ref 5 ↵5 [-0.5 -3 ]  5 1  'CMAESHess'  fc6 ↵5 [-0.5 -3 ]  5 1  'CMAES'  BigGAN ↵Network broke the connection...  '}

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-12032021-004';
preMeta(iExp).expControlFN = '210312_Alfa_generate_BigGAN(2)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2021-03-12-Alfa-02\2021-03-12-12-10-06' ) ;
preMeta(iExp).comments = '5 [-0.5 -3 ]  5 1  ''CMAESHess''  fc6 ↵5 [-0.5 -3 ]  5 1  ''CMAES''  BigGAN' ;


iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-01042021-002';
preMeta(iExp).expControlFN = '210401_Alfa_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2021-04-01-Alfa-01\2021-04-01-12-30-13' ) ;
preMeta(iExp).comments = '002 Generate BigGAN starts 12:30 ↵FC6 BigGAN ↵9 [-1.5 -1.5] 5 1  ↵CMAES_Hessian, FC6 ↵CMAES BigGAN_class' ;


iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-08042021-002';
preMeta(iExp).expControlFN = '210408_Alfa_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2021-04-08-Alfa-01\2021-04-08-11-48-20' ) ;
preMeta(iExp).comments = 'Chan 7 [-1.3 -1.4] 4deg unit 1↵Chan 7 [-1.3 -1.4] 4deg unit 1↵CMAES Hessian FC6↵CMAES BigGAN' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-08042021-004';
preMeta(iExp).expControlFN = '210408_Alfa_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2021-04-08-Alfa-02\2021-04-08-12-42-30' ) ;
preMeta(iExp).comments = 'Chan 2 [-1.4 -4] 5deg unit 1 | Chan 2 [-1.4 -4] 5deg unit 1 CMAES Hessian FC6  CMAES BigGAN' ;


iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-15042021-002';
preMeta(iExp).expControlFN = '210415_Alfa_generate_BigGAN';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2021-04-15-Alfa-01\2021-04-15-12-06-14' ) ;
preMeta(iExp).comments = 'Generate BigGAN↵26  [-0.9  -4.5]  5  1' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-16042021-003';
preMeta(iExp).expControlFN = '210416_Alfa_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2021-04-16-Alfa-01\2021-04-16-11-33-33' ) ;
preMeta(iExp).comments = '20  [-3  -4]  7  1  CMAES_Hessian 20  [-3  -4]  7  1  CMAES' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-12052021-004';
preMeta(iExp).expControlFN = '210512_Alfa_generate_BigGAN(4)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2021-05-12-Alfa-01\2021-05-12-12-48-19' ) ;
preMeta(iExp).comments = '30  [-5.3 -0.9]  5  1   ↵30  [-5.3 -0.9]  5  1   ↵CMAES Hessian↵BigGAN Class  CMAES' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-14052021-003';
preMeta(iExp).expControlFN = '210514_Alfa_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2021-05-14-Alfa-01\2021-05-14-13-11-03' ) ;
preMeta(iExp).comments = '003 Generate BigGAN↵29 [-0.3 -2.0] 5 1↵Hash 5/5↵FC6 CMAES Hessian↵BigGAN class CMAES↵Good BigGAN evolution' ;

iExp = iExp + 1; 
preMeta(iExp).ephysFN  =     'Alfa-02062021-004';
preMeta(iExp).expControlFN = '210602_Alfa_generate_BigGAN(1)';
preMeta(iExp).stimuli = fullfile('n:','Stimuli\2020-BigGAN\2021-06-02-Alfa-01\2021-06-02-12-48-00' ) ;
preMeta(iExp).comments = '004 Generate BigGAN FC6↵22 [     0    -3]  1  6 deg  2/5 SU↵FC6↵BigGAN class↵Wow this guy is really a black pointy ear detector!' ;

% iExp = iExp + 1; 
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
    
end

A.meta =        meta;
A.rasters =     rasters;
A.Trials =      Trials ;
A.lfps =        lfps ;
