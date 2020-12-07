%% ImDistCorr

%%
MatStatsDir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(MatStatsDir, Animal+"_HessBGEvolStats.mat"), 'HEStats')
%%
HECorrImDist = repmat(struct(),1,14);
%%
for Expi = 1:14
iThr = 2;
fprintf("Analyze Exp %d",Expi)
HECorrImDist(Expi).units = HEStats(Expi).units;
% Load the data
gen_rspmat = cell2mat(HEStats(Expi).evol.rspmat(iThr,:));
% gen_rspmat = squeeze(cell2mat(reshape(HEStats(Expi).evol.psth(iThr,:),1,1,[]))); % cell by trial
% gen_rspmat_col = cellfun(@(idx) mean(rasters(:, 51:200, idx),[2]), gen_idx_seq, 'Uni', 0);
distmat = HEStats(Expi).evol.distmat_BG;
distmat_nan = distmat + diag(nan(1,size(distmat,1)));
for iCh = 1:size(gen_rspmat,1)
    rsp_vec = gen_rspmat(iCh,:);
    [cc_vec_all, p_vec_all] = corrcoef_vec(rsp_vec', distmat); % Common Pearson 
    [cc_all_min, minidx] = min(cc_vec_all);
    HECorrImDist(Expi).evol.squ.corr(iCh) = cc_all_min;
    HECorrImDist(Expi).evol.squ.p_cc(iCh) = p_vec_all(minidx);
    [cc_vec, p_vec] = corrcoef_vec(rsp_vec', distmat, 'Type', 'Pearson');
    [cc_min, minidx] = min(cc_vec);
    fprintf("Chan%s %.3f %.1e\n", HEStats(Expi).units.unit_name_arr(iCh),cc_min,p_vec(minidx))
    HECorrImDist(Expi).evol.squ_noc.corr(iCh) = cc_min;
    HECorrImDist(Expi).evol.squ_noc.p_cc(iCh) = p_vec(minidx);
%     h=figure(3); clf;
%     scatter(distmat(:,minidx), rsp_vec')
%     XLIM = xlim();xlinsp = [0:0.025:XLIM(2)];hold on
%     gprMdl = fitrgp(distmat(:,minidx), rsp_vec');
%     gpr_fit = gprMdl.predict(xlinsp');
%     plot(xlinsp, gpr_fit)
%     xlabel("LPIPS distance");ylabel("Firing Rate")
%     title([compose("Evolution Exp %d (pref Ch%d) Chan%s %.3f %.1e", Expi, HEStats(Expi).units.pref_chan(2), ...
%          HEStats(Expi).units.unit_name_arr(iCh),cc_min,p_vec(minidx))])
%     pause
end
end
%%
MatStatsDir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
savefast(fullfile(MatStatsDir, Animal+"_HECorrImDist.mat"), 'HECorrImDist')

function [cc_vec, p_vec] = corrcoef_vecnan(X, Y, varargin)
cc_vec = []; p_vec = [];
for i = 1:length(Y)
    [cc, P] = corr([X, Y(:,i)],'Rows','complete', varargin{:}); %'Type', 'Pearson'
    cc_vec(i) = cc(1,2);
    p_vec(i) = P(1,2);
end
end