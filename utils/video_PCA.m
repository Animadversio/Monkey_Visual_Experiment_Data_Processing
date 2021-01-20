webcam_dir = "N:\Data-WebCam";
figdir = "O:\MonkVidPCA";
vid = VideoReader(fullfile(webcam_dir,"Video 4.wmv"));
fps = vid.FrameRate;
%% Crop the relevant part of frame.
[I2 RECT] = imcrop(frameRGB);% RECT=[149.5    0.5  276  258]
%% Collect Video Frames RGB (cropped) and do PCA.
fr_pix_mat = [];
vid = VideoReader(fullfile(webcam_dir,"Video 4.wmv"));
while hasFrame(vid)
frameRGB = readFrame(vid);
frRGBCrp = imcrop(frameRGB,RECT);
fr_pix_mat = [fr_pix_mat;reshape(frRGBCrp,1,[])];
end
fr_size = size(frRGBCrp);
frN = size(fr_pix_mat,1);
%% PCA is expensive....
% [coeff,score,latent,tsquared,explained] = pca(single(fr_pix_mat/255.0),'NumComponents',50);
%% Compute Sparse SVD of the Frames. Takes long time
[U,S,V] = svds(double(fr_pix_mat) / 255.0, 50);
S = diag(S);
%%
save(fullfile(figdir,"video_PCA.mat"),"U","S","V",'fr_size')
%% SVD spectrum
expvar = S(2:end).^2 / sum(S(2:end).^2); % get rid of first SV as it's the mean
figure(11); set(11,'pos',[1000         192         560         786]);
T = tiledlayout(2,1,"Pad",'compact','TileSp','compact');
nexttile(1)
plot(expvar);
xlabel("Sing Vect # (Exclud. 1st SV)")
ylabel("Expl. Var.")
% box off
nexttile(2)
plot(cumsum(expvar));
xlabel("Sing Vect # (Exclud. 1st SV)")
ylabel("Cum. Expl. Var.")
title(T,compose("SVD Spectrum of Frame Image (RGB)\n%d frame, %dpix value",size(U,1),size(V,1)))
% box off
saveas(11,fullfile(figdir,"Vid_SVD_spectrum.png"))
savefig(11,fullfile(figdir,"Vid_SVD_spectrum.fig"))
%% Time Series of the Projection
figure(12);
T=tiledlayout(1,1,'Pad','compact');
nexttile(T,1)
timeline = (1:size(U,1))/fps;
plot(timeline,U(:,1:7),"LineWidth",1.2)
xlabel("Time in Exp(sec)")
ylabel("PC Projection")
title("PC Projection of Frames. Time Serie")
legend(compose("PC%d",1:7))
saveas(12,fullfile(figdir,"Frame_PCcoef_timeseries.png"))
savefig(12,fullfile(figdir,"Frame_PCcoef_timeseries.fig"))
%% Image Sequence of Moving along the top PCs
PC_imgs = reshape(V,[fr_size,size(V,2)]);
refframe = 100*PC_imgs(:,:,:,1)*sign(mean(U(:,1),1));%double(frRGBCrp)/255.0;
figure(6);
for PCi=1:30
amp = 10*std(U(:,PCi))*S(PCi);
imshow([refframe-amp*PC_imgs(:,:,:,PCi), refframe-amp/2*PC_imgs(:,:,:,PCi),refframe, refframe+amp/2*PC_imgs(:,:,:,PCi), refframe+amp*PC_imgs(:,:,:,PCi)]) %255*PC_imgs(:,:,:,5))
title(compose("PC%d Effects: PC1+-10 std PC%d",PCi,PCi))
saveas(6,fullfile(figdir,compose("PC%dpattern.png",PCi)))
pause
end
%% Movie of Moving along the top PCs
for PCi=1:30
amp = 10*std(U(:,PCi))*S(PCi);
imgcol = arrayfun(@(x) refframe - x*PC_imgs(:,:,:,PCi), linspace(-amp,amp,5),'uni',0);
imgs2gif(imgcol,fullfile(figdir,compose("PC%dmovie.gif",PCi)),0.1)
end
%% Montages of the Top PC Dev vectors
figure(10);
PCamps = std(U,1,1).*S(:)'.*sign(mean(U,1)) .* [10,150*ones(1,49)];
montage(PC_imgs.*reshape(PCamps,1,1,1,[])); 
saveas(10,fullfile(figdir,compose("TopPC_Montages.png")))
figure(10);
PCamps = std(U,1,1).*S(:)'.*sign(mean(U,1)) .* [10,-150*ones(1,49)];
montage(PC_imgs.*reshape(PCamps,1,1,1,[])); 
saveas(10,fullfile(figdir,compose("TopPC_neg_Montages.png")))
%%
%%
opticFlow = opticalFlowFarneback('NeighborhoodSize',7);
%('NoiseThreshold',0.0005);%opticalFlowHS;%opticalFlowFarneback;
% bwfimg = uint8(mean(fimg,3));
% opticFlow.estimateFlow(bwfimg)
h = figure(1);
movegui(h);
hViewPanel = uipanel(h,'Position',[0 0 1 1],'Title','Plot of Optical Flow Vectors');
hPlot = axes(hViewPanel);
while hasFrame(vid)
    frameRGB = readFrame(vid);
    frameGray = rgb2gray(frameRGB);
    flow = estimateFlow(opticFlow,frameGray);
    imshow(frameRGB)
    hold on
    plot(flow,'DecimationFactor',[5 5],'ScaleFactor',35,'Parent',hPlot);
    hold off
    pause(10^-2)
end
%%
%%
function imgs2gif(imgs,gifname,Delay)
% do a frame shot on a figure and append it to the gif collection.
if nargin == 3
    Delay = 0.05;
end
for fi = 1:numel(imgs)
[imind,cm] = rgb2ind(imgs{fi},256); % im = frame2im(imgs{fi}); 
% Write to the GIF File 
if fi == 1 
  imwrite(imind,cm,gifname,'gif', 'Loopcount',inf,'DelayTime',Delay); 
else 
  imwrite(imind,cm,gifname,'gif', 'WriteMode','append','DelayTime',Delay); 
end 
end
end