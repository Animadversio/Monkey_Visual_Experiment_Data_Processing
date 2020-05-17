%% Evolution Select Axis
system("subst S: E:\Network_Data_Sync")
system("subst N: E:\Network_Data_Sync")
mat_dir = "C:\Users\binxu\OneDrive - Washington University in St. Louis\Mat_Statistics";
Animal = "Alfa";
load(fullfile(mat_dir,"Alfa_ManifPopDynamics.mat"));
load(fullfile(mat_dir,"Alfa_Manif_stats.mat"));
load(fullfile(mat_dir,"Alfa_Evol_stats.mat"));
%% 
global G 
G = FC6Generator("matlabGANfc6.mat");
evofig = figure(4);
evoldata = struct('mask', mask, 'Addmask', 0, "curcode", basis(1,:), "scoremap", scoremap);
% You want multi-views of your data. 
% Neural view (PSTH), Image View, Code Space View
Expi = 3;geni = 1;
EStats(Expi).evol.psth
guidata(evofig, evoldata);
GenSLD
RankSLD
LocalPCASLD

function playEvol_Callback(hObj, evt)
data = guidata(hObj);

end
function drawPSTH

end