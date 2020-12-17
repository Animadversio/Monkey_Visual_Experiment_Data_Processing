figure
autocorr(double(spikeChans(70,:)))
%%
figroot = "O:\SpkTimeScale";
Animal = "Alfa"; Set_Path;
Bigmatdir = "S:\Data-Ephys-MAT";

figdir = fullfile(figroot,meta.ephysFN);
mkdir(figdir);
%%
intvl = 20;
spkCh_movsum = movsum(spikeChans, intvl, 2); % Moving sum / mean
spkCh_movsum = spkCh_movsum(:,1:intvl:end); % subsampling to make it non-overlapping
tic
fitStat = [];
H = figure;T=tiledlayout('flow');
for iCh = 1:size(spikeChans,1)
y = spkCh_movsum(iCh,:);
[normedACF, lags, bounds] = autocorr(double(y),'NumLags',20);
unnormACF = normedACF*var(y,1);
lags = lags*intvl;
fitType = @(A,B,tau,x) A.*(B + exp(-x/tau));
[fitparam,gof] = fit(lags(2:end)',normedACF(2:end)',fitType,...
        'StartPoint', [1, 0.01, 15], ...
        'Lower', [0, 0, 0], ...
        'Robust', 'LAR' );
gof.iCh = iCh;
gof.chan = meta.spikeID(iCh);
if isfield(meta, "unitID"), gof.unit = meta.unitID(iCh); end
gof.A = fitparam.A;
gof.B = fitparam.B;
gof.tau = fitparam.tau;
fitStat = [fitStat;gof];
nexttile();%'flow'
hold on
plot(lags(2:end),normedACF(2:end))
plot(lags(2:end),feval(fitparam, lags(2:end)))
title(compose("Ch%d tau=%.1f\nA=%.1f B=%.1f (r2=%.3f)",gof.chan,gof.tau,gof.A,gof.B,gof.rsquare))
hold off
end
title(T,meta.ephysFN)
toc
saveas(H,fullfile(figdir, compose("%s.png",meta.ephysFN)))
%%
fitTab = struct2table(fitStat);
validmsk = fitTab.rsquare>0.6;
ITmsk = fitTab.chan<33;
V4msk = fitTab.chan>48;
V1msk = fitTab.chan>32 & fitTab.chan<49;
%%
[~,P,~,TSTAT] = ttest2(fitTab.tau(ITmsk&validmsk),fitTab.tau(V1msk&validmsk))
[~,P,~,TSTAT] = ttest2(fitTab.tau(ITmsk&validmsk),fitTab.tau(V4msk&validmsk))
[~,P,~,TSTAT] = ttest2(fitTab.tau(V4msk&validmsk),fitTab.tau(V1msk&validmsk))
figure;
scatter(fitTab.chan,fitTab.tau,36)%,fitTab.unit
vline(32.5);vline(48.5)
figure;
scatter(fitTab.chan,fitTab.rsquare,36)%,fitTab.unit
ylim([-1,1]);vline(32.5);vline(48.5);
%%
