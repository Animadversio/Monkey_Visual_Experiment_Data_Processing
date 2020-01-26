% As reading the imageNames in Trials in *.Mat file is very hard in Python. 
% We try to export the imageNames into csv files so that it's easy to read
% in python. 
%% Output Stimuli image names
locMATPath = "S:\Data-Ephys-MAT";
netMATPath = "\\storage1.ris.wustl.edu\crponce\Active\Data-Ephys-MAT";
EphsFN = "Beto64chan-30102019-001";
FN = EphsFN+"_formatted.mat";
if exist(fullfile(locMATPath,FN))==2
    data = load(fullfile(locMATPath,FN),'Trials');
elseif exist(fullfile(locMATPath,FN))==2
    data = load(fullfile(netMATPath,FN),'Trials');
else
    error("Formatted mat file doesn't exist")
end
%%
imgName = string(data.Trials.imageName);
fidOut=fopen(fullfile(locMATPath,EphsFN+"_imgName.csv"),'w');
for k=1:size(imgName,1)
fprintf(fidOut,"%s\n",imgName(k));
end
fidOut=fclose(fidOut);
