% Spec file to code the norm of the sphere for each experiment and the id
% for the preferred channel

system("subst S: D:\Network_Data_Sync") % set this alias! so that copying and syncing could work 
% it will load the newest version of ExpSpecTable and compute pref_chan_arr
% and norm_arr from it! 
ExpSpecTable_Aug = readtable("S:\ExpSpecTable_Augment.xls");
pref_chan_arr = [];
for iExp = 1:max(ExpSpecTable_Aug.Expi)
    rowi = ExpSpecTable_Aug.Expi == iExp & contains(ExpSpecTable_Aug.expControlFN, 'selectivity');
    if iExp == 35
        rowi = find(rowi,1);
    end
    pref_chan_arr = [pref_chan_arr, ExpSpecTable_Aug.pref_chan(rowi)];
end
norm_arr = [];
for iExp = 1:max(ExpSpecTable_Aug.Expi)
    rowi = ExpSpecTable_Aug.Expi == iExp & contains(ExpSpecTable_Aug.expControlFN, 'selectivity');
    path = ExpSpecTable_Aug.stimuli(rowi);
    path = path{1};
    imglist = struct2table(dir(path));
    idx = find(contains(imglist.name, 'norm'));
    if isempty(idx)
        warning("No norm file found cannot get the norm value.")
    end
    norm = regexp(imglist.name(idx(1)), "norm_(\d*)_PC2", 'tokens');
    norm = str2num(norm{1}{1}{1});
    norm_arr = [norm_arr, norm];
end
result_dir = "C:\\Users\\ponce\\OneDrive - Washington University in St. Louis\\PC_space_tuning";
% norm_arr = [328, 326, 269, 329, 401, 386, 300, 332, 421, 316]; % from 429, it's PC23 + Pasupathy
% pref_chan_arr = [29, 6, 5, 20, 19, 13, 28, 15, 17, 63]; % 1st 10 experiment majorly IT, 1 V4 channel
% norm_arr = [328,326,269,329,401,386,300,332,421,316,429,350,296,354,288,298,309,342,359,292,406,348,308,393,424,319,408,297,302,300,297,310,296]; % Automatic generated
% pref_chan_arr = [29,6,5,20,19,13,28,15,17,63,26,20,40,42,35,36,34,37,49,60,56,54,58,50,57,55,26,39,39,45,45,48,48];

%%
Pasu_norm_arr = [429, 350, 296, 354, 288, 298, 309, 342, 359, 292, 406, 348, 308, 393, 424, 319]; % experiment involving Pasupathy patches
Pasu_pref_chan_arr = [26, 20, 40, 42, 35, 36, 34, 37, 49, 60, 56, 54, 58, 50, 57, 55]; % a few IT channel and lots of V1 channel 
 % a folder to store the files. 