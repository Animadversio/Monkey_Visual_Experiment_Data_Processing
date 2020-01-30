%% Set Path 
system("subst S: D:\Network_Data_Sync") % set this alias! so that copying and syncing could work 
result_dir = "C:\\Users\\ponce\\OneDrive - Washington University in St. Louis\\PC_space_tuning";
% it will load the newest version of ExpSpecTable and compute pref_chan_arr
% and norm_arr from it! 
ExpSpecTable_Aug = readtable("S:\ExpSpecTable_Augment.xls");
copyfile("S:\ExpSpecTable_Augment.xls", ".\ExpSpecTable_Augment.xls")
addpath(".\utils")
addpath(".\DNN")

% Depends on 
% - brewermap
% - shadedErrorBar