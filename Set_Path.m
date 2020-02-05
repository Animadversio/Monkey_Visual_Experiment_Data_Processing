%% Set Path 
system("subst S: D:\Network_Data_Sync") % set this alias! so that copying and syncing could work 
%result_dir = "C:\\Users\\ponce\\OneDrive - Washington University in St. Louis\\PC_space_tuning";
% it will load the newest version of ExpSpecTable and compute pref_chan_arr
% and norm_arr from it! 
ExpSpecTable_Aug = readtable("S:\ExpSpecTable_Augment.xlsx");
copyfile("S:\ExpSpecTable_Augment.xlsx", ".\ExpSpecTable_Augment.xlsx")
addpath(".\utils")
addpath(".\DNN")
addpath(".\NMF")
% Depends on 
% - brewermap
% - shadedErrorBar