% Make final figures out of this
tabdir = "O:\Manif_Fitting\Kent_summary";
alfatab = readtable(fullfile(tabdir,"Alfa_Kstats.csv"));
betotab = readtable(fullfile(tabdir,"Beto_Kstats.csv"));
%%
alltab = [];
Animal_tab = array2table(repmat("Alfa",size(alfatab,1),1),'VariableNames',{'Animal'});
alltab = [alltab; Animal_tab, alfatab];
Animal_tab = array2table(repmat("Beto",size(betotab,1),1),'VariableNames',{'Animal'});
alltab = [alltab; Animal_tab, betotab];
%%
validmsk = ~((alltab.Animal=="Alfa")&(alltab.Expi==10));
figure;
histogram(alltab.R2(validmsk),20)

%% Show only the real driver units in the channel 
%  Compare the R2 histogram

