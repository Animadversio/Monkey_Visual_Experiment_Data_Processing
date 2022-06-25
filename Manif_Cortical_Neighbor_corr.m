
load("E:\OneDrive - Washington University in St. Louis\Mat_Statistics\Alfa_CortiDisCorr.mat")
%%\
Animal = "Alfa";
for Expi = 1:numel(CortiDisCorr)
M = struct();
areavec = area_map(CortiDisCorr(Expi).units.spikeID,"Alfa");
M.ITmsk = areavec=="IT"; M.V4msk = areavec=="V4"; M.V1msk = areavec=="V1";
Fmsk = [CortiDisCorr(Expi).FStats.F_P]'<0.01;

prefchan = CortiDisCorr(Expi).units.pref_chan;
Cortdist = CortiDisCorr(Expi).cortexDismat;
samemsk = abs(Cortdist - 0) < 1E-5;
neighbormsk = abs(Cortdist - 400) < 1E-5;
nonneighbormsk = (~neighbormsk)&(~samemsk)&~isnan(Cortdist);
% diagmsk = nan(size(Cortdist));
%%
nonneighcorrvec = CortiDisCorr(Expi).avgsph_corrmat(nonneighbormsk);
% corrvec  = CortiDisCorr(Expi).avgsph_corrmat(neighbormsk);
%%
fprintf("%s Exp%d chan %d\n",Animal,Expi, prefchan)
for area = ["V1","V4","IT"]
areamsk = M.(area+"msk");
unitmsk = areamsk & Fmsk;
ITcorrvec = CortiDisCorr(Expi).avgsph_corrmat((unitmsk*unitmsk') & neighbormsk);
ITnoncorrvec = CortiDisCorr(Expi).avgsph_corrmat((unitmsk*unitmsk') & nonneighbormsk);
ttest2_print(ITcorrvec, ITnoncorrvec,area+" neighbor",area+" nonneighbor");
end
end