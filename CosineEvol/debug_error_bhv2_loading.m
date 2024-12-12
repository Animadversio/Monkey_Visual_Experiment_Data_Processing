% debug error loading 
mlread("N:\Data-Behavior (BHV2)\220311_Beto_rfMapper_basic.bhv2")
mlread("N:\Data-Behavior (BHV2)\220311_Beto_selectivity_basic.bhv2")
mlread("N:\Data-Behavior (BHV2)\220311_Beto_generate_BigGAN_cosine(4).bhv2")

%%
[B,MLConfig,TrialRecord] = mlread("N:\Data-Behavior (BHV2)\210629_Beto_howSimilar(1).bhv2");% success 
[B,MLConfig,TrialRecord] = mlread("N:\Data-Behavior (BHV2)\210813_Beto_rfMapper_basic(1).bhv2");% FAIL!!!
[B,MLConfig,TrialRecord] = mlread("N:\Data-Behavior (BHV2)\211025_Beto_rfMapper_basic.bhv2");% FAIL!!!
[B,MLConfig,TrialRecord] = mlread("N:\Data-Behavior (BHV2)\220810_Beto_selectivity_basic(1).bhv2");%FAIL!!!
[B,MLConfig,TrialRecord] = mlread("N:\Data-Behavior (BHV2)\220810_Beto_selectivity_basic(1).bhv2");%FAIL!!!
%% Solution 
%  https://www.notion.so/binxus-mind-palace/Debug-notes-for-unable-to-load-bhv2-file-with-mlread-in-monkeylogic-d16a374f5a684fd3ad5c1206c6b0e409?pvs=4
%  Basically make sure the mlread is from a Monkeylogic installation with version higher than 2.2
%%
% https://monkeylogic.nimh.nih.gov/board/read.php?3,1328
% https://monkeylogic.nimh.nih.gov/board/read.php?2,1166,1166#msg-1166
setpref('NIMH_MonkeyLogic','DefaultEncoding','windows-1252')  % See fopen for all available encoding schemes
setpref('NIMH_MonkeyLogic','MachineFormat','ieee-le')