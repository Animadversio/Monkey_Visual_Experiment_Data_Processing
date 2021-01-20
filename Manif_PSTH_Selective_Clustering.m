%% 
% Selectively Cluster and Visualize the PSTH

global ManifDyn_P ManifDyn Stats EStats
global ExpTabPool
mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
Animal = "Beto"; 
load(fullfile(mat_dir, Animal+'_Evol_stats.mat')) 
load(fullfile(mat_dir, Animal+'_Manif_stats.mat')) 
load(fullfile(mat_dir, Animal+'_ManifPopDynamics.mat'),'ManifDyn')
load(fullfile(mat_dir, Animal+'_ManifBasis.mat'), 'ManifBasisStats')
load(fullfile(mat_dir, Animal+'_Manif_SummaryStats.mat'),'SummaryStats')
load(fullfile(mat_dir, Animal+'_Manif_PrecSummaryStats.mat'), 'PrecSummaryStats', 'wdw_vect')
savepath = "E:\OneDrive - Washington University in St. Louis\ManifDynCoding";
pooldir = "E:\OneDrive - Washington University in St. Louis\ManifExpPool";
mkdir(savepath)
%% Form a Pooled table of all the units in all the experiments, pool their stats for manifold exp
tic
ExpTabPool = table();
% a template for nan structure when there is nothing there.
nanstruct = struct('F', nan,'F_P', nan,'F_bsl', nan,'F_P_bsl', nan,'T', nan,'t_P', nan,'F_wdw', nan(1,size(wdw_vect,1)),'F_P_wdw', nan(1,size(wdw_vect,1)));
for Animal = ["Beto","Alfa"] 
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'), 'EStats') 
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'), 'Stats') 
load(fullfile(mat_dir, Animal+'_Manif_PrecSummaryStats.mat'), 'PrecSummaryStats', 'wdw_vect')
for Expi = 1:numel(Stats)
    if Animal=="Alfa" && Expi==10,continue; end
    nCh = numel(Stats(Expi).meta.spikeID);
    expstattab = table();
    expstattab.Animal = repmat(Animal, nCh, 1);
    expstattab.Expi = repmat(Expi, nCh, 1);
    expstattab.iCh =  [1:nCh]';
    expstattab.chan = Stats(Expi).meta.spikeID;
    expstattab.unit_num = Stats(Expi).units.unit_num_arr;
    expstattab.chan_str = Stats(Expi).units.unit_name_arr;
    expstattab.area = chan2area(Stats(Expi).meta.spikeID);
    expstattab.prefchan = repmat(Stats(Expi).units.pref_chan, nCh, 1);
    expstattab.prefchan_ui = repmat(EStats(Expi).evol.unit_in_pref_chan, nCh, 1); % unit in pref chan
    expstattab.imgpos = repmat(EStats(Expi).evol.imgpos, nCh, 1); % position of image
    expstattab.imgsize = repmat(EStats(Expi).evol.imgsize, nCh, 1); % position of image
    manifstattab = struct2table(PrecSummaryStats(Expi).manif{1});
    if isempty(PrecSummaryStats(Expi).pasu)
        pasustattab = struct2table(repmat(nanstruct,nCh,1));
    else
        pasustattab = struct2table(PrecSummaryStats(Expi).pasu);
    end
    if isempty(PrecSummaryStats(Expi).gabor)
        gaborstattab = struct2table(repmat(nanstruct,nCh,1));
    else
        gaborstattab = struct2table(PrecSummaryStats(Expi).gabor);
    end
    pasustattab.Properties.VariableNames = cellfun(@(V)"pasu"+V, pasustattab.Properties.VariableNames);
    gaborstattab.Properties.VariableNames = cellfun(@(V)"gab"+V, gaborstattab.Properties.VariableNames);
    ExpTabPool = [ExpTabPool; [expstattab, manifstattab, pasustattab, gaborstattab]];
end
end
writetable(ExpTabPool, fullfile(pooldir,"ManifExpPool.csv"))
writetable(ExpTabPool, fullfile(mat_dir,"ManifExpPool.csv"))
toc
%%
load(fullfile(mat_dir, "Alfa"+'_ManifPopDynamics.mat'),'ManifDyn')
ManifDyn_P.Alfa = ManifDyn;
load(fullfile(mat_dir, "Beto"+'_ManifPopDynamics.mat'),'ManifDyn')
ManifDyn_P.Beto = ManifDyn;
%%

%%
figure(3);clf;hold on;set(3,'pos',[240         432        1198         546])
histogram(ExpTabPool.gabF(ExpTabPool.area=="V1"),0:0.2:12)
histogram(ExpTabPool.gabF(ExpTabPool.area=="V4"),0:0.2:12)
histogram(ExpTabPool.gabF(ExpTabPool.area=="IT"),0:0.2:12)
[md,idx] = min(abs(ExpTabPool.gabF_P-0.01));
Fthresh = ExpTabPool.gabF(idx);
line([Fthresh,Fthresh],ylim())
legend(["V1","V4","IT","P=0.01 thresh"]);xlabel("F");title("Tuning F on gabor space")
saveas(3,fullfile(pooldir,Animal+"gabor_pooled_F.png"))
%%
figure(4);clf;hold on;set(4,'pos',[240         432        1198         546])
histogram(ExpTabPool.pasuF(ExpTabPool.area=="V1"),0:0.2:12)
histogram(ExpTabPool.pasuF(ExpTabPool.area=="V4"),0:0.2:12)
histogram(ExpTabPool.pasuF(ExpTabPool.area=="IT"),0:0.2:12)
[md,idx] = min(abs(ExpTabPool.pasuF_P-0.01));
Fthresh = ExpTabPool.pasuF(idx);
line([Fthresh,Fthresh],ylim())
legend(["V1","V4","IT","P=0.01 thresh"]);xlabel("F");title("Tuning F on Pasupathy space")
saveas(4,fullfile(pooldir,Animal+"pasu_pooled_F.png"))
%%
figure(5);clf;hold on;set(5,'pos',[240         432        1198         546])
histogram(ExpTabPool.F(ExpTabPool.area=="V1"),0:0.2:12)
histogram(ExpTabPool.F(ExpTabPool.area=="V4"),0:0.2:12)
histogram(ExpTabPool.F(ExpTabPool.area=="IT"),0:0.2:12)
[md,idx] = min(abs(ExpTabPool.F_P-0.01));
Fthresh = ExpTabPool.F(idx);
line([Fthresh,Fthresh],ylim())
legend(["V1","V4","IT","P=0.01 thresh"]);xlabel("F");title("Tuning F on Manifold space")
saveas(5,fullfile(pooldir,Animal+"manif_pooled_F.png"))

%%
global tunedrows
tunedrows = find(ExpTabPool.F_P<0.001);
rowi = datasample(tunedrows, 1);
Expi = ExpTabPool.Expi(rowi); iCh = ExpTabPool.iCh(rowi); Sbj = ExpTabPool.Animal(rowi);
titlestr = compose("Row %d\n%s Exp %d prefchan %d Ch %s\nF=%.1f T=%.1f",rowi,ExpTabPool.Animal(rowi),ExpTabPool.Expi(rowi),ExpTabPool.prefchan(rowi),ExpTabPool.chan_str(rowi),ExpTabPool.F(rowi),ExpTabPool.T(rowi));
% PSTHDynViewer(ManifDyn(Expi).psth_tsr(iCh,:,:,:), 0.01, 10, titlestr)
PSTHDynViewer(ManifDyn_P.(Sbj)(Expi).psth_tsr(iCh,:,:,:), 0.012, 10, titlestr)
%% 
fprintf("Tuned V1 unit %d/%d, Alfa %d Beto %d\n",sum(ExpTabPool.area(tunedrows)=="V1"),sum(ExpTabPool.area=="V1"),sum(ExpTabPool.area(tunedrows)=="V1"&ExpTabPool.Animal(tunedrows)=="Alfa"),sum(ExpTabPool.area(tunedrows)=="V1"&ExpTabPool.Animal(tunedrows)=="Beto"))
fprintf("Tuned V4 unit %d/%d, Alfa %d Beto %d\n",sum(ExpTabPool.area(tunedrows)=="V4"),sum(ExpTabPool.area=="V4"),sum(ExpTabPool.area(tunedrows)=="V4"&ExpTabPool.Animal(tunedrows)=="Alfa"),sum(ExpTabPool.area(tunedrows)=="V4"&ExpTabPool.Animal(tunedrows)=="Beto"))
fprintf("Tuned IT unit %d/%d, Alfa %d Beto %d\n",sum(ExpTabPool.area(tunedrows)=="IT"),sum(ExpTabPool.area=="IT"),sum(ExpTabPool.area(tunedrows)=="IT"&ExpTabPool.Animal(tunedrows)=="Alfa"),sum(ExpTabPool.area(tunedrows)=="IT"&ExpTabPool.Animal(tunedrows)=="Beto"))
%%
gabrows = find(ExpTabPool.gabF_P<0.001);
pasurows = find(ExpTabPool.pasuF_P<0.001);
fprintf("Tuning to Gabor Patches\n")
fprintf("Tuned V1 unit %d/%d, Alfa %d Beto %d\n",sum(ExpTabPool.area(gabrows)=="V1"),sum(ExpTabPool.area=="V1"),sum(ExpTabPool.area(gabrows)=="V1"&ExpTabPool.Animal(gabrows)=="Alfa"),sum(ExpTabPool.area(gabrows)=="V1"&ExpTabPool.Animal(gabrows)=="Beto"))
fprintf("Tuned V4 unit %d/%d, Alfa %d Beto %d\n",sum(ExpTabPool.area(gabrows)=="V4"),sum(ExpTabPool.area=="V4"),sum(ExpTabPool.area(gabrows)=="V4"&ExpTabPool.Animal(gabrows)=="Alfa"),sum(ExpTabPool.area(gabrows)=="V4"&ExpTabPool.Animal(gabrows)=="Beto"))
fprintf("Tuned IT unit %d/%d, Alfa %d Beto %d\n",sum(ExpTabPool.area(gabrows)=="IT"),sum(ExpTabPool.area=="IT"),sum(ExpTabPool.area(gabrows)=="IT"&ExpTabPool.Animal(gabrows)=="Alfa"),sum(ExpTabPool.area(gabrows)=="IT"&ExpTabPool.Animal(gabrows)=="Beto"))
fprintf("Tuning to Pasupathy Patches\n")
fprintf("Tuned V1 unit %d/%d, Alfa %d Beto %d\n",sum(ExpTabPool.area(pasurows)=="V1"),sum(ExpTabPool.area=="V1"),sum(ExpTabPool.area(pasurows)=="V1"&ExpTabPool.Animal(pasurows)=="Alfa"),sum(ExpTabPool.area(pasurows)=="V1"&ExpTabPool.Animal(pasurows)=="Beto"))
fprintf("Tuned V4 unit %d/%d, Alfa %d Beto %d\n",sum(ExpTabPool.area(pasurows)=="V4"),sum(ExpTabPool.area=="V4"),sum(ExpTabPool.area(pasurows)=="V4"&ExpTabPool.Animal(pasurows)=="Alfa"),sum(ExpTabPool.area(pasurows)=="V4"&ExpTabPool.Animal(pasurows)=="Beto"))
fprintf("Tuned IT unit %d/%d, Alfa %d Beto %d\n",sum(ExpTabPool.area(pasurows)=="IT"),sum(ExpTabPool.area=="IT"),sum(ExpTabPool.area(pasurows)=="IT"&ExpTabPool.Animal(pasurows)=="Alfa"),sum(ExpTabPool.area(pasurows)=="IT"&ExpTabPool.Animal(pasurows)=="Beto"))


%% Collect the tuned rows from 2 monkeys
% selDynMaps = arrayfun(@(ri) ManifDyn(ExpTabPool.Expi(ri)).psth_tsr(ExpTabPool.iCh(ri),:,:,:),tunedrows,'Uni',0);
selDynMaps = arrayfun(@(ri) ManifDyn_P.(ExpTabPool.Animal(ri))(ExpTabPool.Expi(ri)).psth_tsr(ExpTabPool.iCh(ri),:,:,:),tunedrows,'Uni',0);
selDynMaps = cell2mat(selDynMaps);
%%
smSelDynMaps = smoothdata(selDynMaps,2,'gaussian',7);
smSelDynMaps = smoothdata(smSelDynMaps,3,'gaussian',3);
smSelDynMaps = smoothdata(smSelDynMaps,4,'gaussian',3);
smSelDynMat = reshape(smSelDynMaps,numel(tunedrows),[]);
%%
[reduction,umap,clusterIdentifiers,extras] = run_umap(smSelDynMat, 'n_neighbors',25, "metric","correlation",'n_epochs',1000);
%%
areacode = containers.Map(["IT","V4","V1"],[49,25,16]);
figure(8);
S = scatter(reduction(:,1),reduction(:,2),arrayfun(@(ar)areacode(ar),ExpTabPool.area(tunedrows)),ExpTabPool.Expi(tunedrows));
S.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Subj", ExpTabPool.Animal(tunedrows));
S.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Exp", ExpTabPool.Expi(tunedrows));
S.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Area", ExpTabPool.area(tunedrows));
S.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Chan", ExpTabPool.chan_str(tunedrows));
S.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Driver", ExpTabPool.prefchan(tunedrows));
S.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Pos", compose("(%.1f %.1f)",ExpTabPool.imgpos(tunedrows,:)));
S.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Size", ExpTabPool.imgsize(tunedrows));
S.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("F", ExpTabPool.F(tunedrows));
S.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("T", ExpTabPool.T(tunedrows));
S.DataTipTemplate.DataTipRows(1:2)=[];
S.ButtonDownFcn = @VisDynamics;
colorbar()
%%
[reduction_L,umap_L,clusterIdentifiers_L,extras_L] = run_umap(smSelDynMat, 'n_neighbors',75, "metric","correlation",'n_epochs',1000);
areacode = containers.Map(["IT","V4","V1"],[49,25,16]);
%%
figure(10);
S = scatter(reduction_L(:,1),reduction_L(:,2),arrayfun(@(ar)areacode(ar),ExpTabPool.area(tunedrows)),ExpTabPool.Expi(tunedrows));
S.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Exp", ExpTabPool.Expi(tunedrows));
S.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Area", ExpTabPool.area(tunedrows));
S.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Chan", ExpTabPool.chan_str(tunedrows));
S.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Driver", ExpTabPool.prefchan(tunedrows));
S.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Pos", compose("(%.1f %.1f)",ExpTabPool.imgpos(tunedrows,:)));
S.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Size", ExpTabPool.imgsize(tunedrows));
S.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("F", ExpTabPool.F(tunedrows));
S.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("T", ExpTabPool.T(tunedrows));
S.DataTipTemplate.DataTipRows(1:2)=[];
colorbar()
%%
meanSelDynMaps = mean(selDynMaps(:,50:200,:,:),2);
meanSelDynMat = reshape(meanSelDynMaps,numel(tunedrows),[]);
%%
[reduc_m,umap_m,clusterId_m,extras_m] = run_umap(meanSelDynMat, 'n_neighbors',25, "metric","correlation",'n_epochs',1000);
%%
areacode = containers.Map(["IT","V4","V1"],[49,25,16]);
figure(12);
S = scatter(reduc_m(:,1),reduc_m(:,2),arrayfun(@(ar)areacode(ar),ExpTabPool.area(tunedrows)),ExpTabPool.imgsize(tunedrows));
S.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Subj", ExpTabPool.Animal(tunedrows));
S.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Exp", ExpTabPool.Expi(tunedrows));
S.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Area", ExpTabPool.area(tunedrows));
S.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Chan", ExpTabPool.chan_str(tunedrows));
S.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Driver", ExpTabPool.prefchan(tunedrows));
S.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Pos", compose("(%.1f %.1f)",ExpTabPool.imgpos(tunedrows,:)));
S.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("Size", ExpTabPool.imgsize(tunedrows));
S.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("F", ExpTabPool.F(tunedrows));
S.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("T", ExpTabPool.T(tunedrows));
S.DataTipTemplate.DataTipRows(1:2)=[];
S.ButtonDownFcn = @VisDynamics;
colorbar()
%%
% array2table(nan(67,8),'VariableNames',{'F','F_P','F_bsl','F_P_bsl','T','t_P','F_wdw','F_P_wdw'})
function area_arr = chan2area(chans)
area_arr = strings(size(chans));
area_arr(chans<=32)="IT";
area_arr(33<=chans & chans<=48)="V1";
area_arr(49<=chans)="V4";
end
function [coordinateSelected, rowi] = VisDynamics(hObj, event)
global ManifDyn_P ExpTabPool tunedrows
x = hObj.XData; 
y = hObj.YData; 
pt = event.IntersectionPoint(1:2);       % The (x0,y0) coordinate you just selected
coordinates = [x(:),y(:)];     % matrix of your input coordinates
dist = pdist2(pt,coordinates); 
[~, rel_rowi] = min(dist);
rowi = tunedrows(rel_rowi);
Expi = ExpTabPool.Expi(rowi); iCh = ExpTabPool.iCh(rowi); Subj = ExpTabPool.Animal(rowi);
titlestr = compose("row %d %s Exp %d prefCh %d Ch %s\n pos %s size %d F=%.1f T=%.1f",rowi,...
    ExpTabPool.Animal(rowi),ExpTabPool.Expi(rowi),ExpTabPool.prefchan(rowi),ExpTabPool.chan_str(rowi),compose("(%.1f,%.1f)",ExpTabPool.imgpos(rowi,:)),ExpTabPool.imgsize(rowi),ExpTabPool.F(rowi),ExpTabPool.T(rowi));
fprintf("Row %d in table, %s\n",rowi,titlestr)
PSTHDynViewer(ManifDyn_P.(Subj)(Expi).psth_tsr(iCh,:,:,:), 0.01, 20, titlestr)
end
function PSTHDynViewer(psth_map, sleep, h, titlestr)
% PSTH map is a 200 by 11 by 11 array
if nargin==1, sleep=0.01; h=[]; titlestr=""; end
if nargin==2, h=[]; titlestr=""; end
if nargin==3, titlestr=""; end
if ndims(psth_map)==4, psth_map=squeeze(psth_map); end
Wlen=20; shift_step=2;
% gifname = fullfile(savepath, compose('%s_Exp%d_manif%s_calign.gif',Animal,Expi,Unit_str));
start_arr =  0:shift_step:200-Wlen;
act_map_col = cell2mat(arrayfun(@(strt)mean(psth_map(strt+1:strt+Wlen,:,:), [1]),start_arr','Uni',false)); % length(start_arr), 11, 11
act_map_col = permute(act_map_col,[2,3,1]); % 11, 11, length(start_arr)
CMAX = prctile(act_map_col(:),98);
CMIN = prctile(act_map_col(:),2.5);
if isempty(h); h = figure(1); else, h = figure(h); end
set(0,"CurrentFigure",h); set(h,'position',[50   350   552   542]); 
imsc = imagesc(-90:18:90, -90:18:90, act_map_col(:,:,1));
caxis([CMIN CMAX]); colorbar();
axis image; ylabel("PC 2 degree");xlabel("PC 3 degree")
title(compose("ps: [%d,%d] ms",1,20))
fi = 1;
for start = 0:shift_step:200-Wlen
wdw = start+1:start+Wlen;
imsc.CData = act_map_col(:, :, fi); fi=fi+1;
title(compose("%s\nps: [%d,%d] ms",titlestr,wdw(1),wdw(end))) % ,Animal,Expi,Unit_str,
drawnow; pause(sleep);
end
end