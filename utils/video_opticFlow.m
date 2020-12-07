webcam_dir = "N:\Data-WebCam";
figdir = "O:\MonkVidPCA";
vid = VideoReader(fullfile(webcam_dir,"Video 4.wmv"));
fps = vid.FrameRate;
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