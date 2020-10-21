
%% Movie Selectivity Analysis, Trial by Trial
Animal="Alfa";Set_Path;
ftr = contains(ExpRecord.ephysFN,"21102020");
[meta_new,rasters_new,~,Trials_new] = loadExperiments(find(ftr),Animal);

%%
[meta_,rasters_,lfps_,Trials_] = loadData('Alfa-21102020-006','expControlFN','201021_Alfa_selectivity_movie(2)', 'rasterWindow',[-250 2500]); %'sdf', 'raster') ;

[meta_,rasters_,lfps_,Trials_] = loadData('Alfa-21102020-005','expControlFN','201021_Alfa_selectivity_movie(1)', 'rasterWindow',[-250 2500]); %'sdf', 'raster') ;
