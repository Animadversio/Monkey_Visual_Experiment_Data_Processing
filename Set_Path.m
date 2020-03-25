%% Set Path 
system("subst S: E:\Network_Data_Sync") % set this alias! so that copying and syncing could work 
%result_dir = "C:\\Users\\ponce\\OneDrive - Washington University in St. Louis\\PC_space_tuning";
% it will load the newest version of ExpSpecTable and compute pref_chan_arr
% and norm_arr from it! 
ExpSpecTable_Aug = readtable("S:\ExpSpecTable_Augment.xlsx");
ExpSpecTable_Aug_alfa = readtable("S:\Exp_Record_Alfa.xlsx");
copyfile("S:\ExpSpecTable_Augment.xlsx", ".\ExpSpecTable_Augment.xlsx")
copyfile("S:\Exp_Record_Alfa.xlsx", ".\Exp_Record_Alfa.xlsx")
addpath(".\utils")
addpath(".\DNN")
addpath(".\NMF")
% Depends on 
% - brewermap
% - shadedErrorBar
switch Animal
    case "Alfa"
        ExpRecord = ExpSpecTable_Aug_alfa;
    case "Beto"
        ExpRecord = ExpSpecTable_Aug;
    case "Both"
        ExpRecord = [ExpSpecTable_Aug; ExpSpecTable_Aug_alfa];
end 