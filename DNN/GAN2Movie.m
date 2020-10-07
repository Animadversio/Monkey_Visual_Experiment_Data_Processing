%% movie from GAN
%  There is a python version achieving the same thing! 
py.torch.cuda.empty_cache()
G = torchStyleGAN2('ffhq-512-avg-tpurun1.pt');
py.torch.cuda.empty_cache()
%%
vec1 = randn(1,512);
vec2 = randn(1,512);
interpV = SLERP(vec1, vec2, linspace(0,1,72));
imgs = G.visualize(interpV, 1);
% figure; montage(imgs)
%%
interpV = SLERP(vec1, vec2, linspace(0,1,72));
imgs = G.visualize(interpV, 1);
%%
interpV = TangSLERP(vec1, vec2, linspace(-0.2,0.2,72));
imgs = G.visualize(interpV, 1);
%%
workingDir = "E:\";
outputVideo = VideoWriter(fullfile(workingDir,'Face512.avi'));
outputVideo.FrameRate = 24; %shuttleVideo.FrameRate;
open(outputVideo)
% Loop through the image sequence, and then write it to the video.
for ii = 1:size(imgs,4)
   writeVideo(outputVideo,imgs(:,:,:,ii))
end
close(outputVideo) % Finalize the video file.
winopen(fullfile(workingDir,'Face512.avi'))
%%
outputVideo = VideoWriter(fullfile(workingDir,'Face512_rep.avi'));
outputVideo.FrameRate = 24; %shuttleVideo.FrameRate;
open(outputVideo)
imgn = size(imgs,4);
% Loop through the image sequence, and then write it to the video.
for ii = [1:imgn, imgn, imgn, imgn, imgn:-1:1, 1, 1, 1]
   writeVideo(outputVideo,imgs(:,:,:,ii))
end
close(outputVideo) % Finalize the video file.
winopen(fullfile(workingDir,'Face512_rep.avi'))
%% Sinusoidally modulated GAN image. 
sinticks = 0.2*sin(linspace(-pi/2,pi/2,72));
sin_interpV = TangSLERP(vec1, vec2, sinticks);
imgs = G.visualize(sin_interpV, 1);
%%
outputVideo = VideoWriter(fullfile(workingDir,'Face512_sin.mp4'),'MPEG-4');
% outputVideo = VideoWriter(fullfile(workingDir,'Face512_sin.avi'));
outputVideo.FrameRate = 24; %shuttleVideo.FrameRate;
open(outputVideo)
imgn = size(imgs,4);
% Loop through the image sequence, and then write it to the video.
for ii = [1:imgn, imgn:-1:1, 1:imgn]
   writeVideo(outputVideo,imgs(:,:,:,ii))
end
close(outputVideo) % Finalize the video file.
winopen(fullfile(workingDir,'Face512_sin.mp4'))
% winopen(fullfile(workingDir,'Face512_sin.avi'))
%%
py.torch.cuda.empty_cache()
G = torchStyleGAN2('ffhq-512-avg-tpurun1.pt');
py.torch.cuda.empty_cache()
%%
Hdata = py.numpy.load("E:\OneDrive - Washington University in St. Louis\Hessian_summary\StyleGAN2_Fix"+...
                        "\ffhq-512-avg-tpurun1_fix\H_avg_ffhq-512-avg-tpurun1_fix.npz");
eva_avg = Hdata.get('eva_avg').double;
evc_avg = Hdata.get('evc_avg').double;
evc_avg = evc_avg(:,end:-1:1);
eva_avg = eva_avg(end:-1:1);
%% Tune the exploration range here
eva_avg(1:20).^-0.4 / 3
%%
refvec = randn(1, 512);
for eigi = 1:20
maxang = eva_avg(eigi).^-0.4/3; 
imgs = GANVideo(G, refvec, evc_avg(:,eigi)', compose("Face512_eig%03d_sin",eigi), maxang, 72);
end
%%
refvec = randn(1, 512);
for eigi = 1:20
maxang = eva_avg(eigi).^-0.4 * 6; 
imgs = GANVideo(G, refvec, evc_avg(:,eigi)', compose("Face512_eig%03d_lin_sin",eigi), maxang, 72, false);
end
%% BigGAN version
py.torch.cuda.empty_cache()
G = torchBigGAN();
py.torch.cuda.empty_cache()
EmbedMat = G.get_embedding();
%%
Hdata = py.numpy.load("E:\OneDrive - Washington University in St. Louis\Hessian_summary\BigGAN\H_avg_1000cls.npz");
eva_all = Hdata.get('eigvals_avg').double;
evc_all = Hdata.get('eigvects_avg').double;
evc_all = evc_all(:,end:-1:1);
eva_all = eva_all(end:-1:1);
eva_cls = Hdata.get('eigvals_clas_avg').double;
evc_cls = Hdata.get('eigvects_clas_avg').double;
evc_cls = evc_cls(:,end:-1:1);
eva_cls = eva_cls(end:-1:1);
eva_nos = Hdata.get('eigvals_nois_avg').double;
evc_nos = Hdata.get('eigvects_nois_avg').double;
evc_nos = evc_nos(:,end:-1:1);
eva_nos = eva_nos(end:-1:1);
evc_nos_ag = [evc_nos;zeros(128)];
evc_cls_ag = [zeros(128);evc_cls];

%%
refvec = [0.7*randn(1,128), EmbedMat(:, 11)'];
G.space = "all";
img = G.visualize(refvec);
figure;imshow(img)
%%
clsid = 11;
refvec = [0.7*randn(1,128), EmbedMat(:, 11)'];
G.space = "all";
%
for eigi = 1:10
maxang = eva_nos(eigi).^-0.4 * 6; 
imgs = GANVideo(G, refvec, evc_nos_ag(:,eigi)', compose("BigGAN_cls%03d_eig%03d_lin_sin",clsid,eigi), maxang, 72, false);
end
%%
clsid = 374;
refvec = [0.7*randn(1,128), EmbedMat(:, clsid)'];
%%
G.space = "all";
for eigi = 1:10
maxang = eva_nos(eigi).^-0.4/3; 
imgs = GANVideo(G, refvec, evc_nos_ag(:,eigi)', compose("BigGAN_cls%03d_eig%03d_sin",clsid,eigi), maxang, 72, true);
end
%%
clsid = 374;
% refvec = [0.7*randn(1,128), EmbedMat(:, clsid)'];
% G.space = "all";
scalers = eva_cls.^-0.2*1;
for eigi = 1:10
maxang = scalers(eigi);% eva_cls(eigi).^-0.3*3.5; 
imgs = GANVideo(G, refvec, evc_cls_ag(:,eigi)', compose("BigGAN_cls%03d_clseig%03d_sin",clsid,eigi), maxang, 72, false);
end
%%
function imgs = GANVideo(G, vec1, vec2, namestr, maxang, ticks, sphere)
if nargin <= 5, sphere = true; end % spherical vs linear exploration
if nargin <= 4, ticks=72; end % n samples in a pass 
if nargin <= 3, maxang = 0.2; end % maximum exploration, angle for spherical
workingDir = "E:\tmp";
sinticks = maxang*sin(linspace(-pi/2,pi/2,ticks+1));
if sphere
sin_interpV = TangSLERP(vec1, vec2, sinticks);
else
sin_interpV = TangLERP(vec1, vec2, sinticks);
end 
imgs = G.visualize(sin_interpV);
outputVideo = VideoWriter(fullfile(workingDir, namestr+".avi"));%,'MPEG-4');
% outputVideo = VideoWriter(fullfile(workingDir,'Face512_sin.avi'));
outputVideo.FrameRate = 24; %shuttleVideo.FrameRate;
open(outputVideo)
imgn = size(imgs,4);
% Loop through the image sequence, and then write it to the video.
for ii = [[ticks/2:imgn, imgn:-1:1, 1:ticks/2+1], [ticks/2:imgn, imgn:-1:1, 1:ticks/2+1]]
   writeVideo(outputVideo,imgs(:,:,:,ii))
end
close(outputVideo) % Finalize the video file.
winopen(fullfile(workingDir, namestr+".avi"))
end