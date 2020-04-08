mat_dir = "C:\Users\binxu\OneDrive - Washington University in St. Louis\Mat_Statistics";
Animal = "Beto";
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'))
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'))
%%
expconv = fittype( @(A, B, tau, gen) A - (A-B) .* exp(- (gen) ./ tau), ...
    'independent', {'gen'},'dependent',{'score'});
sigmoid = fittype( @(A, B, t0, tau, gen) B + A ./ (1 + exp(- (gen - t0) ./ tau)), ...
    'independent', {'gen'},'dependent',{'score'});
%%
result_dir = "C:\Users\binxu\OneDrive - Washington University in St. Louis\Evol_Fitting";
Traj_stats = repmat(struct(), 1, length(EStats));
for i = 1:length(EStats)
%    savepath = fullfile(result_dir, sprintf("%s_Exp%02d", Animal, EStats(i).Expi));
%    mkdir(savepath)
   channel_j = EStats(i).units.pref_chan_id(1);
   fprintf("%s Exp%d Channel %s\n", EStats(i).Animal, EStats(i).Expi, ...
            EStats(i).units.unit_name_arr{channel_j});
   Traj_stats(i).Expi = EStats(i).Expi;
   Traj_stats(i).chan = EStats(i).units.pref_chan;
   Traj_stats(i).chan_id = EStats(i).units.pref_chan_id;
   Traj_stats(i).unit = EStats(i).units.unit_num_arr(EStats(i).units.pref_chan_id);
   
   act_traj = cellfun(@(psth) squeeze(mean(psth(1,51:200,:),[1,2])), EStats(i).evol.psth,'UniformOutput',false);
   meanact_traj = cellfun(@mean, act_traj,'UniformOutput',true);
   [H,P,CI,STATS] = ttest2(cat(1,act_traj{end-2:end-1}), cat(1,act_traj{2:3}));
   Traj_stats(i).P_act = P;
   Traj_stats(i).t_act = STATS.tstat;
   score_traj = cellfun(@(psth) squeeze(mean(psth(1,51:200,:),[1,2])) - squeeze(mean(psth(1,1:40,:),[1,2])), EStats(i).evol.psth,'UniformOutput',false);
   meanscore_traj = cellfun(@mean, score_traj,'UniformOutput',true);
   [H,P,CI,STATS] = ttest2(cat(1,score_traj{end-2:end-1}), cat(1,score_traj{2:3}));
   Traj_stats(i).P_scr = P;
   Traj_stats(i).t_scr = STATS.tstat;
   
   [expfit,expgof] = fit([1:length(meanact_traj)]', meanact_traj',expconv,...
                'StartPoint', [max(meanact_traj), min(meanact_traj), length(meanact_traj)/2], ...
                'Lower', [min(meanact_traj)-10, min(meanact_traj)-10,  0], ...
                'Upper', [ max(meanact_traj)+50,  max(meanact_traj)+50,  2*length(meanact_traj)]);
    expgof.coefname = string(coeffnames(expfit)');
    expgof.coef = coeffvalues(expfit);
    expgof.confint = confint(expfit);
    [sigmfit,sigmgof] = fit([1:length(meanact_traj)]', meanact_traj',sigmoid,...
                'StartPoint', [max(meanact_traj), min(meanact_traj), 0, length(meanact_traj)/2], ...
                'Lower', [min(meanact_traj)-10, min(meanact_traj)-10,  0, 0], ...
                'Upper', [ max(meanact_traj)+50,  max(meanact_traj)+50,  2*length(meanact_traj), 2*length(meanact_traj)]);
    sigmgof.coefname = string(coeffnames(sigmfit)');
    sigmgof.coef = coeffvalues(sigmfit);
    sigmgof.confint = confint(sigmfit);
    Traj_stats(i).expfit = expgof;
    Traj_stats(i).sigmfit = sigmgof;
end
%%
Traj_Tab = repmat(struct(), 1, length(EStats));
for i = 1:length(EStats)
   Traj_Tab(i).Expi = EStats(i).Expi;
   Traj_Tab(i).chan = EStats(i).units.pref_chan;
   Traj_Tab(i).chan_id = EStats(i).units.pref_chan_id;
   Traj_Tab(i).unit = EStats(i).units.unit_num_arr(EStats(i).units.pref_chan_id);
   Traj_Tab(i).P_act= Traj_stats(i).P_act;
   Traj_Tab(i).t_act= Traj_stats(i).t_act;
   Traj_Tab(i).P_scr= Traj_stats(i).P_scr;
   Traj_Tab(i).t_scr= Traj_stats(i).t_scr;
   Traj_Tab(i).expR2 = Traj_stats(i).expfit.rsquare;
   Traj_Tab(i).expA = Traj_stats(i).expfit.coef(1);
   Traj_Tab(i).expB = Traj_stats(i).expfit.coef(2);
   Traj_Tab(i).exptau = Traj_stats(i).expfit.coef(3);
   Traj_Tab(i).sigmR2 = Traj_stats(i).sigmfit.rsquare;
   Traj_Tab(i).sigmA = Traj_stats(i).sigmfit.coef(1);
   Traj_Tab(i).sigmB = Traj_stats(i).sigmfit.coef(2);
   Traj_Tab(i).sigmt0 = Traj_stats(i).sigmfit.coef(3);
   Traj_Tab(i).sigmtau = Traj_stats(i).sigmfit.coef(4);
end
%
Traj_Tab = struct2table(Traj_Tab)
%%
figure;clf;hold on;%set(7,'position',[712   362   493   616])
jitter = sort(randn(1,size(Traj_Tab,1)) * 0.05);
scatter(1 * ones(1,size(Traj_Tab,1)) + jitter, Traj_Tab.t_act)

%%
% Table is really good for doing statistics (filtering and ordering)
save(fullfile(result_dir, compose("%s_Trajstats.mat", Animal)), 'Traj_Tab', 'Traj_stats')
writetable(Traj_Tab, fullfile(result_dir, compose("%s_Trajstats.csv", Animal)))
%%
expconv = fittype( @(A, B, tau, t) A - (A-B) .* exp(- (t) ./ tau), ...
    'independent', {'t'},'dependent',{'y'});
[expfit,expgof] = fit([1:length(meanact_traj)]', meanact_traj',expconv,...
                'StartPoint', [max(meanact_traj), min(meanact_traj), length(meanact_traj)/2], ...
                'Lower', [min(meanact_traj)-10, min(meanact_traj)-10,  0], ...
                'Upper', [ max(meanact_traj)+50,  max(meanact_traj)+50,  2*length(meanact_traj)]);
expgof.coefname = string(coeffnames(expfit)');
expgof.coef = coeffvalues(expfit);
expgof.confint = confint(expfit);
figure;hold on;
plot(expfit);
plot(meanact_traj);
%%
sigmoid = fittype( @(A, B, t0, tau, gen) B + A ./ (1 + exp(- (gen - t0) ./ tau)), ...
    'independent', {'gen'},'dependent',{'score'});
[sigmfit,sigmgof] = fit([1:length(meanact_traj)]', meanact_traj',sigmoid,...
                'StartPoint', [max(meanact_traj), min(meanact_traj), 0, length(meanact_traj)/2], ...
                'Lower', [min(meanact_traj)-10, min(meanact_traj)-10,  0, 0], ...
                'Upper', [ max(meanact_traj)+50,  max(meanact_traj)+50,  2*length(meanact_traj), 2*length(meanact_traj)]);
sigmgof.coefname = string(coeffnames(sigmfit)');
sigmgof.coef = coeffvalues(sigmfit);
sigmgof.confint = confint(sigmfit);
figure;plot(sigmfit);hold on;plot(meanact_traj);
            %%
