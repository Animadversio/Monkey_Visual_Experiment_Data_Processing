%% 
global ManifDyn Stats EStats
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
Animal = "Beto"; 
load(fullfile(mat_dir, Animal+'_Evol_stats.mat')) 
load(fullfile(mat_dir, Animal+'_Manif_stats.mat')) 
load(fullfile(mat_dir, Animal+'_ManifPopDynamics.mat'),'ManifDyn')
load(fullfile(mat_dir, Animal+'_ManifBasis.mat'), 'ManifBasisStats')
load(fullfile(mat_dir, Animal+'_Manif_SummaryStats.mat'),'SummaryStats')
savepath = "E:\OneDrive - Washington University in St. Louis\ManifDynCoding";
mkdir(savepath)
%%
G = FC6Generator();
%% Dynamics of `prototypes`
% What is the weight for center of mass. 
% - The activation. 
% - The activation - baseline
% - Activation - baseline - std
% - Activation - percentile of current frame.
Expi = 1; iCh = 1;
dyntsr = squeeze(ManifDyn(Expi).psth_tsr(iCh,:,:,:));
PSTHDynViewer(dyntsr, 0.01, 10)
% PSTHDynViewer(ManifDyn(Expi).psth_tsr(21,:,:,:), 0.02, 10)
%%

%%
PSTHDynViewer(dyntsr, 0.01, 10)
[PHI,THETA] = meshgrid(-90:18:90,-90:18:90);
coordtsr = cat(3,cosd(THETA).*cosd(PHI),... % PC1
                 sind(THETA).*cosd(PHI),... % PC2
                 sind(PHI));  % PC3
baseline = mean(dyntsr(1:50,:,:),'all');
basestd = std(mean(dyntsr(1:50,:,:),1),0,'all');
% Choose the weight tsr protocol: 
% weightsr = max(0,dyntsr - baseline - basestd);
weightsr = max(0, dyntsr - max(prctile(dyntsr,50,[2,3]), baseline + basestd)); 
comtrace = squeeze(sum(weightsr .* shiftdim(coordtsr,-1), [2, 3])) ./ sum(weightsr,[2,3]);
THEtrace = atan2d(comtrace(:,2),comtrace(:,1));
PHItrace = asind(comtrace(:,3));
% [azim,elev,radi] = cart2sph(comtrace(:,1),comtrace(:,2),comtrace(:,3));
figure(6);clf;hold on 
xlim([-90,90]);ylim([-90,90])
xlabel("PC3");ylabel("PC2");
comet(gca,PHItrace,THEtrace)
%% Visualize trajectory of images correspond to the peak of the field.
basis = ManifBasisStats(Expi).basis{1};
sphere_norm = ManifBasisStats(Expi).sphere_norm;
CoM_codes = sphere_norm * comtrace * basis;
img_traj = G.visualize(CoM_codes);
%%  
figure(7); sleep=0.03;
ims = imshow(img_traj(:,:,:,1));
for fi=1:200
ims.CData = img_traj(:,:,:,fi);
title(compose("%d ms Theta %.1f Phi %.1f",fi,THEtrace(fi),PHItrace(fi)))
drawnow;pause(sleep)
end
%% Correlating dynamics across image, or correlating selectivity across time
dynmat = reshape(dyntsr,200,[]);
imgsimmat = corr(dynmat);
timsimmat = corr(dynmat');
figure(8);
subplot(121);imagesc(imgsimmat);colorbar(); xlabel("image id"); axis image
line([0.5,121.5]'+zeros(1,12),[0:11:121;0:11:121]+0.5,'color','red')
subplot(122);imagesc(timsimmat);colorbar(); xlabel("time"); axis image
%% nnmf Factorization applied to spatial temporal data. 
%%
[Tfact,Sfact] = nnmf(dynmat,1);
residue = dynmat - Tfact*Sfact;
L1res = mean(abs(residue),'all');
L2res = mean(residue.^2,'all');
%%
[Tfact2,Sfact2] = nnmf(dynmat,2);
residue2 = dynmat - Tfact2*Sfact2;
L1res2 = mean(abs(residue2),'all');
L2res2 = mean(residue2.^2,'all');
%%
[Tfact3,Sfact3] = nnmf(dynmat,3);
residue3 = dynmat - Tfact3*Sfact3;
L1res3 = mean(abs(residue3),'all');
L2res3 = mean(residue3.^2,'all');
%% nnmf Factorization applied to spatial temporal data. 
FN = 3;
[TfactN,SfactN] = nnmf(dynmat,FN);
residueN = dynmat - TfactN*SfactN;
L1resN = mean(abs(residueN),'all');
L2resN = sqrt(mean(residueN.^2,'all'));
resVar = var(residueN,1,'all');
allVar = var(dynmat,1,'all');
figure(12);
T = tiledlayout(2,FN);
title(T,compose("Residue L1:%.1f L2 %.1f ExpVar %.1f/%.1f",L1resN,L2resN,resVar,allVar))
for Fi=1:FN
ax = nexttile(T,Fi);
plot(TfactN(:,Fi))
ax = nexttile(T,FN+Fi);
imagesc(reshape(SfactN(Fi,:),11,[]))
axis image
colorbar()
end
%%
PSTHDynViewer(dyntsr, 0.01, 10)
basestd = std(mean(dyntsr(1:50,:,:),1),0,'all');
weightsr = max(0,dyntsr - baseline - basestd);
% weightsr = max(0,dyntsr - prctile(dyntsr,40,[2,3]));
comtrace = squeeze(sum(weightsr .* shiftdim(coordtsr,-1), [2, 3])) ./ sum(weightsr,[2,3]);
THEtrace = atan2d(comtrace(:,2),comtrace(:,1));
PHItrace = asind(comtrace(:,3));
% [azim,elev,radi] = cart2sph(comtrace(:,1),comtrace(:,2),comtrace(:,3));
%%
figure(6);clf;hold on 
xlim([-90,90]);ylim([-90,90])
xlabel("PC3");ylabel("PC2");
comet(gca,PHItrace,THEtrace)
%%
figure(5);comet3(comtrace(:,1),comtrace(:,2),comtrace(:,3))
xlabel("PC1");ylabel("PC2");zlabel("PC3")
%%
basis = ManifBasisStats(Expi).basis{1};
sphere_norm = ManifBasisStats(Expi).sphere_norm;
CoM_codes = sphere_norm * comtrace * basis;
img_traj = G.visualize(CoM_codes);
%%
figure(7);sleep=0.02;
ims = imshow(img_traj(:,:,:,1));
for fi=1:200
ims.CData = img_traj(:,:,:,fi);
title(compose("%d ms Theta %.1f Phi %.1f",fi,THEtrace(fi),PHItrace(fi)))
drawnow;pause(sleep)
end
%% Instantaneous F stats (no time window averaging, only SDF spike convolution.)
%  Channels with very good quality can have a significant F stats without
%  any time averaging. (Still trial averaging.) 
Expi = 1;
repmap = shiftdim(ManifDyn(Expi).rep_num,-1);
totrep = sum(repmap,'all');
Fthrs = finv([0.99,0.995,0.999],120,totrep-121);
for iCh = 1:size(ManifDyn(Expi).psth_std_tsr,1)
stdmap = squeeze(ManifDyn(Expi).psth_std_tsr(iCh,:,:,:));
meamap = squeeze(ManifDyn(Expi).psth_tsr(iCh,:,:,:));
meanall = sum(meamap.*repmap,[2,3])/totrep;
SSwi = sum(stdmap.^2 .* (repmap-1),[2,3]);
SSbt = sum((meamap - meanall).^2 .* repmap,[2,3]);
Ftr = SSbt./SSwi / (121-1) * (totrep-121);
figure(4);plot(Ftr);ylim([0.6,2.5]);ylabel('F stats')
line([0,199]'+zeros(1,3),[Fthrs;Fthrs])
title(Stats(Expi).units.unit_name_arr(iCh))
legend(["F stat each ms","0.99 thresh","0.995 thresh","0.999 thresh"])
pause;
end
% for Expi = 20:numel(Stats)
% for si = 1:length(Stats(Expi).manif.psth) % space idx
% for ui = 1:length(Stats(Expi).units.pref_chan_id) % units idx in pref chan
%%

%%
Expi = 11;
ui = Stats(Expi).units.pref_chan_id(1);
PSTHDynViewer(ManifDyn(Expi).psth_tsr(ui,:,:,:))


function DynViewer(Expi, ui, pref)
global ManifDyn Stats
PSTHDynViewer(ManifDyn(Expi).psth_tsr(ui,:,:,:), sleep)
end
function PSTHDynViewer(psth_map, sleep, h)
% PSTH map is a 200 by 11 by 11 array
if nargin==1, sleep=0.01; h=[]; end
if nargin==2, h=[]; end
if ndims(psth_map)==4, psth_map=squeeze(psth_map); end
Wlen=20; shift_step=2;
% gifname = fullfile(savepath, compose('%s_Exp%d_manif%s_calign.gif',Animal,Expi,Unit_str));
start_arr =  0:shift_step:200-Wlen;
act_map_col = cell2mat(arrayfun(@(strt)mean(psth_map(strt+1:strt+Wlen,:,:), [1]),start_arr','Uni',false)); % length(start_arr), 11, 11
act_map_col = permute(act_map_col,[2,3,1]); % 11, 11, length(start_arr)
CMAX = prctile(act_map_col(:),98);
CMIN = prctile(act_map_col(:),2.5);
if isempty(h); h = figure(1); else, h = figure(h); end
set(0,"CurrentFigure",h); set(h,'position',[680   436   552   542]); 
imsc = imagesc(-90:18:90, -90:18:90, act_map_col(:,:,1));
caxis([CMIN CMAX]); colorbar();
axis image; ylabel("PC 2 degree");xlabel("PC 3 degree")
title(compose("ps: [%d,%d] ms",1,20))
fi = 1;
for start = 0:shift_step:200-Wlen
wdw = start+1:start+Wlen;
imsc.CData = act_map_col(:, :, fi); fi=fi+1;
title(compose("ps: [%d,%d] ms",wdw(1),wdw(end))) % ,Animal,Expi,Unit_str,
drawnow; pause(sleep);
end
end
