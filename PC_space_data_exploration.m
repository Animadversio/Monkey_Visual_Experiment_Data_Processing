%% Data Exploration 
%% See the fluctuation of baseline and how that connects with response
%Expi = 1;
%chani = 3;
figure,hold on  
plot(smooth(squeeze(mean(rasters{1}(3,1:40,:), 2))',5))
plot(smooth(squeeze(mean(rasters{1}(3,100:140,:), 2))',5))
plot(smooth(squeeze(mean(rasters{1}(3,100:140,:), 2))' - squeeze(mean(rasters{1}(3,1:40,:), 2))',5))
namelist = unique(Trials{1}.imageName)