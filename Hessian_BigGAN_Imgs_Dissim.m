%%
% The figure correlate the response to BigGAN Hessian Manifold Images
% and the image distance matrix 
%%
Animal="Beto";Set_Path;
ftridx = find(contains(ExpRecord.Exp_collection,"BigGAN_Hessian") & contains(ExpRecord.expControlFN,"select") );
[meta_new,rasters_new,~,Trials_new] = loadExperiments(ftridx, Animal);
%%
MatStatsDir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
Animal="Beto";Set_Path;
load(fullfile(MatStatsDir, Animal+"HessBigGANStats.mat"),"HessBGStats")
%%
CorrImDist = repmat(struct(),1,numel(HessBGStats));
%% Collect all the correlation into a structure! 
for Expi = 1:numel(HessBGStats)
CorrImDist(Expi).units = HessBGStats(Expi).units;
for iCh = 1:size(HessBGStats(Expi).units.spikeID,1)
% Correlation in the Class Space
rsp_vec = reshape(HessBGStats(Expi).class.resp_mat(iCh,:,:),1,[]);
valmsk = ~isnan(rsp_vec); % add this valid mask to get rid of nan
distmat = HessBGStats(Expi).class.distmat;
[cc_vec,p_vec] = corrcoef_vec(rsp_vec(valmsk)',distmat(valmsk,:));
[cc_min, minidx] = min(cc_vec);
fprintf("Chan%s %.3f %.1e\n", HessBGStats(Expi).units.unit_name_arr(iCh),cc_min,p_vec(minidx))
CorrImDist(Expi).class.squ.corr(iCh) = cc_min;
CorrImDist(Expi).class.squ.p_cc(iCh) = p_vec(minidx);

distmat_nan = distmat + diag(nan(1,size(distmat,1)));
[cc_vec,p_vec] = corrcoef_vecnan(rsp_vec(valmsk)',distmat_nan(valmsk,:));
[cc_min, minidx] = min(cc_vec);
fprintf("Chan%s %.3f %.1e\n", HessBGStats(Expi).units.unit_name_arr(iCh),cc_min,p_vec(minidx))
CorrImDist(Expi).class.squ_noc.corr(iCh) = cc_min;
CorrImDist(Expi).class.squ_noc.p_cc(iCh) = p_vec(minidx);

% Correlation in the Noise Space
rsp_vec = reshape(HessBGStats(Expi).noise.resp_mat(iCh,:,:),1,[]);
valmsk = ~isnan(rsp_vec); % add this valid mask to get rid of nan
distmat = HessBGStats(Expi).noise.distmat;
[cc_vec,p_vec] = corrcoef_vec(rsp_vec(valmsk)',distmat(valmsk,:));
[cc_min, minidx] = min(cc_vec);
fprintf("Chan%s %.3f %.1e\n", HessBGStats(Expi).units.unit_name_arr(iCh),cc_min,p_vec(minidx))
CorrImDist(Expi).noise.squ.corr(iCh) = cc_min;
CorrImDist(Expi).noise.squ.p_cc(iCh) = p_vec(minidx);

distmat_nan = distmat + diag(nan(1,size(distmat,1)));
[cc_vec,p_vec] = corrcoef_vecnan(rsp_vec(valmsk)',distmat_nan(valmsk,:));
[cc_min, minidx] = min(cc_vec);
fprintf("Chan%s %.3f %.1e\n", HessBGStats(Expi).units.unit_name_arr(iCh),cc_min,p_vec(minidx))
CorrImDist(Expi).noise.squ_noc.corr(iCh) = cc_min;
CorrImDist(Expi).noise.squ_noc.p_cc(iCh) = p_vec(minidx);
end
end
%%
MatStatsDir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
savefast(fullfile(MatStatsDir, Animal+"_HessCorrImDist.mat"), 'CorrImDist')
%%
Expi = 11;
for iCh = 11%1:size(HessBGStats(Expi).class.resp_mat,1)
% Correlation in the Class Space
rsp_vec = reshape(HessBGStats(Expi).class.resp_mat(iCh,:,:),1,[]);
valmsk = ~isnan(rsp_vec);
distmat = HessBGStats(Expi).class.distmat;
[cc_vec,p_vec] = corrcoef_vec(rsp_vec(valmsk)',distmat(valmsk,:));
[cc_min, minidx] = min(cc_vec);
fprintf("Chan%s %.3f %.1e\n", HessBGStats(Expi).units.unit_name_arr(iCh),cc_min,p_vec(minidx))
CorrImDist(Expi).class.squ.corr(iCh) = cc_min;
CorrImDist(Expi).class.squ.p_cc(iCh) = p_vec(minidx);
% Correlation in the Noise Space
rsp_vec = reshape(HessBGStats(Expi).noise.resp_mat(iCh,:,:),1,[]);
valmsk = ~isnan(rsp_vec);
distmat = HessBGStats(Expi).noise.distmat;
[cc_vec,p_vec] = corrcoef_vec(rsp_vec(valmsk)',distmat(valmsk,:));
[cc_min, minidx] = min(cc_vec);
fprintf("Chan%s %.3f %.1e\n", HessBGStats(Expi).units.unit_name_arr(iCh),cc_min,p_vec(minidx))
CorrImDist(Expi).noise.squ.corr(iCh) = cc_min;
CorrImDist(Expi).noise.squ.p_cc(iCh) = p_vec(minidx);
end


%% Print out the significantly correlated channels
Thr = 1E-4;
for Expi = 1:numel(CorrImDist)
    fprintf("Exp %d class prefCh %d Significant Channels %s\n", Expi, CorrImDist(Expi).units.pref_chan, ...
        join([" ",CorrImDist(Expi).units.unit_name_arr(CorrImDist(Expi).class.squ.p_cc<Thr)']))
    fprintf("Exp %d noise prefCh %d Significant Channels %s\n", Expi, CorrImDist(Expi).units.pref_chan, ...
        join([" ",CorrImDist(Expi).units.unit_name_arr(CorrImDist(Expi).noise.squ.p_cc<Thr)']))
end
%% Print out the significantly correlated channels
Thr = -0.4;
for Expi = 1:numel(CorrImDist)
    fprintf("Exp %d class prefCh %d Significant Channels %s\n", Expi, CorrImDist(Expi).units.pref_chan, ...
        join([" ",CorrImDist(Expi).units.unit_name_arr(CorrImDist(Expi).class.squ.corr<Thr)']))
    fprintf("Exp %d noise prefCh %d Significant Channels %s\n", Expi, CorrImDist(Expi).units.pref_chan, ...
        join([" ",CorrImDist(Expi).units.unit_name_arr(CorrImDist(Expi).noise.squ.corr<Thr)']))
end
%% Visualize the correlation of each channel?
Expi = 2;
for iCh = 1:size(HessBGStats(Expi).class.resp_mat,1)
% Correlation in the Class Space
rsp_vec = reshape(HessBGStats(Expi).class.resp_mat(iCh,:,:),1,[]);
distmat = sqrt(HessBGStats(Expi).class.distmat);
[cc_vec,p_vec] = corrcoef_vec(rsp_vec',distmat);
[cc_min, minidx] = min(cc_vec);
fprintf("Chan%s %.3f %.1e\n", HessBGStats(Expi).units.unit_name_arr(iCh),cc_min,p_vec(minidx))
h=figure(1); clf;
scatter(distmat(:,minidx), rsp_vec')
XLIM = xlim();xlinsp = [0:0.025:XLIM(2)];hold on
gprMdl = fitrgp(distmat(:,minidx), rsp_vec');
gpr_fit = gprMdl.predict(xlinsp');
plot(xlinsp, gpr_fit)
xlabel("LPIPS distance");ylabel("Firing Rate")
title([compose("Exp %d (pref Ch%d) Chan%s %.3f %.1e", Expi, HessBGStats(Expi).units.pref_chan, ...
    HessBGStats(Expi).units.unit_name_arr(iCh),cc_min,p_vec(minidx)),"Class space"])

% Correlation in the Noise Space
rsp_vec = reshape(HessBGStats(Expi).noise.resp_mat(iCh,:,:),1,[]);
distmat = sqrt(HessBGStats(Expi).noise.distmat);
[cc_vec,p_vec] = corrcoef_vec(rsp_vec',distmat);
[cc_min, minidx] = min(cc_vec);
fprintf("Chan%s %.3f %.1e\n", HessBGStats(Expi).units.unit_name_arr(iCh),cc_min,p_vec(minidx))
h=figure(2); clf;
scatter(distmat(:,minidx), rsp_vec')
XLIM = xlim();xlinsp = [0:0.025:XLIM(2)];hold on
gprMdl = fitrgp(distmat(:,minidx), rsp_vec');
gpr_fit = gprMdl.predict(xlinsp');
plot(xlinsp, gpr_fit)
xlabel("LPIPS distance");ylabel("Firing Rate")
title([compose("Exp %d (pref Ch%d) Chan%s %.3f %.1e", Expi, HessBGStats(Expi).units.pref_chan, ...
    HessBGStats(Expi).units.unit_name_arr(iCh),cc_min,p_vec(minidx)),"Noise space"])
pause;
end
%% Plot the Hessian's correlation plot with Evolution's plot
figroot = "E:\OneDrive - Washington University in St. Louis\ImMetricTuning_BigGAN";
set(1,'pos',[184         404        1578         574])
for Expi = 12:numel(CorrImDist)
    % Response vector to evolved images.
    figdir = fullfile(figroot, compose("%s_Exp%02d", Animal, Expi));
    mkdir(figdir)
    evol_rspmat = cell2mat(HEStats(Expi).evol.rspmat(iThr,:)); % Evolution's response matrix
    
    for iCh = 1:numel(CorrImDist(Expi).units.spikeID)%CorrImDist(Expi).units.pref_chan_id' %
    h=figure(1); clf; T = tiledlayout(1,3, 'Padding', 'compact');
    nexttile(1)%subplot(131)
    % Hessian Class space data
    rsp_vec = reshape(HessBGStats(Expi).class.resp_mat(iCh,:,:),1,[]);
    valmsk = ~isnan(rsp_vec);
    distmat = HessBGStats(Expi).class.distmat;
    [cc_vec,p_vec] = corrcoef_vec(rsp_vec(valmsk)',distmat(valmsk,:));
    [cc_min, minidx] = min(cc_vec);
    fprintf("Chan%s class space %.3f %.1e\n", HessBGStats(Expi).units.unit_name_arr(iCh),cc_min,p_vec(minidx))
    scatter(distmat(:,minidx), rsp_vec')
    XLIM = xlim();xlinsp = [0:0.025:XLIM(2)];hold on
    gprMdl = fitrgp(distmat(:,minidx), rsp_vec');
    gpr_fit = gprMdl.predict(xlinsp');
    plot(xlinsp, gpr_fit)
    xlabel("LPIPS distance");ylabel("Firing Rate")
    title([compose("Exp %d (pref Ch%d) Chan%s %.3f %.1e", Expi, HessBGStats(Expi).units.pref_chan, ...
        HessBGStats(Expi).units.unit_name_arr(iCh),cc_min,p_vec(minidx)),"Class space Hessian"])
    % Hessian Noise space data
    rsp_vec = reshape(HessBGStats(Expi).noise.resp_mat(iCh,:,:),1,[]);
    valmsk = ~isnan(rsp_vec); % add this valid mask to get rid of nan
    distmat = HessBGStats(Expi).noise.distmat;
    [cc_vec,p_vec] = corrcoef_vec(rsp_vec(valmsk)',distmat(valmsk,:));
    [cc_min, minidx] = min(cc_vec);
    fprintf("Chan%s noise space %.3f %.1e\n", HessBGStats(Expi).units.unit_name_arr(iCh),cc_min,p_vec(minidx))
    nexttile(2)%subplot(132)
    scatter(distmat(:,minidx), rsp_vec')
    XLIM = xlim();xlinsp = [0:0.025:XLIM(2)];hold on
    gprMdl = fitrgp(distmat(:,minidx), rsp_vec');
    gpr_fit = gprMdl.predict(xlinsp');
    plot(xlinsp, gpr_fit)
    xlabel("LPIPS distance");ylabel("Firing Rate")
    title([compose("Exp %d (pref Ch%d) Chan%s %.3f %.1e", Expi, HessBGStats(Expi).units.pref_chan, ...
        HessBGStats(Expi).units.unit_name_arr(iCh),cc_min,p_vec(minidx)),"Noise space Hessian"])
    % Evolution data
    iCh_Evo = find( (HECorrImDist(Expi).units.unit_num_arr == CorrImDist(Expi).units.unit_num_arr(iCh)) & ...
            (HECorrImDist(Expi).units.spikeID == CorrImDist(Expi).units.spikeID(iCh))); 
    if isempty(iCh_Evo), continue;end
    rsp_vec = evol_rspmat(iCh_Evo,:);
    distmat = HEStats(Expi).evol.distmat_BG;
    [cc_vec_all, p_vec_all] = corrcoef_vec(rsp_vec', distmat); % Common Pearson 
    [cc_all_min, minidx] = min(cc_vec_all);
    fprintf("Chan%s evolution %.3f %.1e\n", HEStats(Expi).units.unit_name_arr(iCh_Evo),cc_all_min,p_vec_all(minidx))
    nexttile(3)%subplot(133)
    scatter(distmat(:,minidx), rsp_vec')
    XLIM = xlim();xlinsp = [0:0.025:XLIM(2)];hold on
    gprMdl = fitrgp(distmat(:,minidx), rsp_vec');
    gpr_fit = gprMdl.predict(xlinsp');
    plot(xlinsp, gpr_fit)
    xlabel("LPIPS distance");ylabel("Firing Rate")
    title([compose("Exp %d (pref Ch%d) Chan%s %.3f %.1e", Expi, HEStats(Expi).units.pref_chan(iThr), ...
         HEStats(Expi).units.unit_name_arr(iCh_Evo),cc_all_min, p_vec_all(minidx)), "Evolution"])
    title(T, "Comparing Correlation of Image Dissimilarity and Firing")
    saveas(h, fullfile(figdir, compose("%s_Exp%02d_pref%d_Ch%s_best.png", Animal, Expi, ...
        HessBGStats(Expi).units.pref_chan, HessBGStats(Expi).units.unit_name_arr(iCh))))
%     saveas(h, fullfile(figroot, compose("%s_Exp%02d_pref%d_Ch%s_best.png", Animal, Expi, ...
%         HessBGStats(Expi).units.pref_chan, HessBGStats(Expi).units.unit_name_arr(iCh))))
%     pause
    end
end

%%
Expi=3; CorrImDist(Expi).units.unit_name_arr(CorrImDist(Expi).class.squ.p_cc<1E-5)

%%
h=figure(1); clf;
scatter(distmat(:,minidx), rsp_vec')
XLIM = xlim();xlinsp = [0:0.025:XLIM(2)];hold on
gprMdl = fitrgp(distmat(:,minidx), rsp_vec');
gpr_fit = gprMdl.predict(xlinsp');
plot(xlinsp, gpr_fit)
%% Experimental: EigenVector of the distance matrix 
distmat = HessBGStats(Expi).class.distmat;
N = size(distmat, 1);
centH = eye(N) - ones(N) / N;
Bmat = centH * -distmat * centH;
Bmat = (Bmat + Bmat')/2;
[Bevc, Beva] = eig(Bmat, 'vector');
%% Test Triangle Inequality
for Expi = 5:numel(HessBGStats)
test_triangle_inequality(HessBGStats(Expi).class.distmat)
test_triangle_inequality(HessBGStats(Expi).noise.distmat)
end
% As a result, the distance matrices we have don't have triangle inequality
%%
for Expi = 1:numel(HessBGStats)
test_triangle_inequality(sqrt(HessBGStats(Expi).class.distmat))
test_triangle_inequality(sqrt(HessBGStats(Expi).noise.distmat))
end
%%

function test_triangle_inequality(distmat)
EPS = 1E-3;
N = size(distmat,1);
for i = 1:N
    for j = i+1:N
        for k = j+1:N
            assert(EPS+distmat(i,j)+distmat(i,k)>=distmat(j,k), "d(%d,%d)+d(%d,%d)>=d(%d,%d) Failed! %.4f + %.4f = %.4f < %.4f", i,j,i,k,j,k,distmat(i,j),distmat(i,k),distmat(i,j)+distmat(i,k),distmat(j,k))
            assert(EPS+distmat(i,j)+distmat(j,k)>=distmat(i,k), "d(%d,%d)+d(%d,%d)>=d(%d,%d) Failed! %.4f + %.4f = %.4f < %.4f", i,j,j,k,i,k,distmat(i,j),distmat(j,k),distmat(i,j)+distmat(j,k),distmat(i,k))
            assert(EPS+distmat(j,k)+distmat(i,k)>=distmat(i,j), "d(%d,%d)+d(%d,%d)>=d(%d,%d) Failed! %.4f + %.4f = %.4f < %.4f", j,k,i,k,i,j,distmat(j,k),distmat(i,k),distmat(j,k)+distmat(i,k),distmat(i,j))
        end
    end
end
end

function [cc_vec, p_vec] = corrcoef_vecnan(X, Y, varargin)
cc_vec = []; p_vec = [];
for i = 1:length(Y)
    [cc, P] = corr([X, Y(:,i)],'Rows','complete', varargin{:}); %'Type', 'Pearson'
    cc_vec(i) = cc(1,2);
    p_vec(i) = P(1,2);
end
end