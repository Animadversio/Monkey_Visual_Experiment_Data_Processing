
function ExpRecord = loadBackupExpRecord(Animal)

%% Set Path  % set this alias! so that copying and syncing could work 
% if ~strcmp(getenv('COMPUTERNAME'), "DESKTOP-MENSD6S")
% system('subst N: \\research.files.med.harvard.edu\Neurobio\PonceLab')
% system("subst S: F:\Network_Data_Sync")
% % system('subst N: \\storage1.ris.wustl.edu\crponce\Active')
% end
% if strcmp(getenv('COMPUTERNAME'), 'LAPTOP-U8TSR4RE')
%     system('subst O: "D:\OneDrive - Washington University in St. Louis"')
%     system("subst S: D:\Network_Data_Sync") % set this alias! so that copying and syncing could work 
% else
%     system('subst O: "E:\OneDrive - Washington University in St. Louis"')
%     system("subst S: F:\Network_Data_Sync") % set this alias! so that copying and syncing could work 
% end
% it will load the newest version of ExpSpecTable and compute pref_chan_arr
% and norm_arr from it! 
if strcmp(getenv('COMPUTERNAME'), "DESKTOP-MENSD6S")  % At home
	fprintf("At home, check the date and version of your experiment Record table.(The one in folder may be the most recent. Loading from it)\n")
    keyboard;
    ExpSpecTable_Aug = readtable("ExpSpecTable_Augment.xlsx");
	ExpSpecTable_Aug_alfa = readtable("Exp_Record_Alfa.xlsx");
    ExpSpecTable_Aug_caos = readtable("Exp_Record_Caos.xlsx");
    ExpSpecTable_Aug_diab = readtable("Exp_Record_Diablito.xlsx");
    system("subst N: E:\Network_Data_Sync")
% 	copyfile("E:\Monkey_Data\ExpSpecTable_Augment.xlsx", ".\ExpSpecTable_Augment.xlsx")
% 	copyfile("E:\Monkey_Data\Exp_Record_Alfa.xlsx", ".\Exp_Record_Alfa.xlsx")
elseif strcmp(getenv('COMPUTERNAME'), 'LAPTOP-U8TSR4RE') % new BInxu laptop
%     keyboard;
    ExpSpecTable_Aug = readtable("ExpSpecTable_Augment.xlsx");
	ExpSpecTable_Aug_alfa = readtable("Exp_Record_Alfa.xlsx");
    ExpSpecTable_Aug_caos = readtable("Exp_Record_Caos.xlsx");
    ExpSpecTable_Aug_diab = readtable("Exp_Record_Diablito.xlsx");
    system("subst N: D:\Network_Data_Sync")
    copyfile(".\ExpSpecTable_Augment.xlsx", "S:\ExpSpecTable_Augment.xlsx")
	copyfile(".\Exp_Record_Alfa.xlsx", "S:\Exp_Record_Alfa.xlsx")
elseif exist("S:\",'dir') % Currently I set up S:\ at home as well, so everything should match
	ExpSpecTable_Aug = readtable("S:\ExpSpecTable_Augment.xlsx",'Format','auto');
	ExpSpecTable_Aug_alfa = readtable("S:\Exp_Record_Alfa.xlsx",'Format','auto');
    ExpSpecTable_Aug_caos = readtable("S:\Exp_Record_Caos.xlsx",'Format','auto');
    ExpSpecTable_Aug_diab = readtable("S:\Exp_Record_Diablito.xlsx",'Format','auto');
    % copy the files to local folder and save as csv for backup
	copyfile("S:\Exp_Record_Alfa.xlsx", ".\Exp_Record_Alfa.xlsx")
	copyfile("S:\ExpSpecTable_Augment.xlsx", ".\ExpSpecTable_Augment.xlsx")
	copyfile("S:\Exp_Record_Caos.xlsx", ".\Exp_Record_Caos.xlsx")
	copyfile("S:\Exp_Record_Diablito.xlsx", ".\Exp_Record_Diablito.xlsx")
    writetable(ExpSpecTable_Aug, ".\Exp_Record_Beto.csv")
    writetable(ExpSpecTable_Aug_alfa, ".\Exp_Record_Alfa.csv")
    writetable(ExpSpecTable_Aug_caos, ".\Exp_Record_Caos.csv")
    writetable(ExpSpecTable_Aug_diab, ".\Exp_Record_Diablito.csv")
else
    fprintf("load local exprecord in folder %s instead\n", pwd)
    keyboard;
	ExpSpecTable_Aug = readtable("ExpSpecTable_Augment.xlsx");
	ExpSpecTable_Aug_alfa = readtable("Exp_Record_Alfa.xlsx");
    ExpSpecTable_Aug_caos = readtable("Exp_Record_Caos.xlsx");
    ExpSpecTable_Aug_diab = readtable("Exp_Record_Diablito.xlsx");
end
if ~ (exist('Animal','var') == 1)
    Animal = "Both";
end
switch Animal
    case "Alfa"
        ExpRecord = ExpSpecTable_Aug_alfa;
    case "Beto"
        ExpRecord = ExpSpecTable_Aug;
    case "Caos"
        ExpRecord = ExpSpecTable_Aug_caos;
    case "Diablito"
        ExpRecord = ExpSpecTable_Aug_diab;
    case "Both"
        ExpRecord = [ExpSpecTable_Aug; ExpSpecTable_Aug_alfa];
    case "All"
        ExpRecord = [ExpSpecTable_Aug_alfa; ExpSpecTable_Aug; ExpSpecTable_Aug_caos; ExpSpecTable_Aug_diab];
    otherwise
        fprintf("Unknown animal %s, loading all\n", Animal);
        ExpRecord = [ExpSpecTable_Aug; ExpSpecTable_Aug_alfa; ExpSpecTable_Aug_caos; ExpSpecTable_Aug_diab];
end 
end
