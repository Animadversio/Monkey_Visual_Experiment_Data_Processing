%% 
Animal = "Beto";Set_Path;
mat_dir = "O:\Mat_Statistics";
load(fullfile(mat_dir, Animal+'_Manif_stats.mat'),'Stats') % Manifold Stats 
load(fullfile(mat_dir, Animal+'_Evol_stats.mat'),'EStats') % Evolution Stats
load(fullfile(mat_dir,Animal+"_ManifRadiTuneStats.mat"),'RadTuneStats') % Collection of Radial Tuning Curves in all Spaces
load(fullfile(mat_dir,Animal+"_EvoRef_ImDist.mat"),"EvoRefImDistStat") % Image Distance between reference images. 
%%
metric_list = ["squ"];%,"SSIM","L2","FC6","FC6_corr"]; % all the metrics we want to use to measure distance 
label_list = ["LPIPS (SqueezeNet)"];%, "SSIM", "L2", "FC6 (L2)", "FC6 (1 - corr)"];
tic
for Expi=11 %1:numel(Stats)
score_vec = cellfun(@(psth)mean(psth(1,51:200,:),'all'),reshape(EStats(Expi).ref.psth_arr,[],1));
score_sem_vec = cellfun(@(psth)mean(psth(1,51:200,:),'all'),reshape(EStats(Expi).ref.psth_arr,[],1));
[sortScore,sortId] = sort(score_vec,'Descend');
[maxScore,maxId] = max(score_vec);
% RadTuneStats(Expi).evoref.maxScore = maxScore;
for mi = 1:numel(metric_list)
metname = metric_list(mi);
dist_vec = EvoRefImDistStat(Expi).(metname)(:,maxId);
figure(1);clf;hold on
scatter(dist_vec, score_vec, 'DisplayName', 'Data Original')
[gpr_fit, ~, AUC_gpr, gprMdl] = GPRfitting(dist_vec, score_vec); % Gaussian Process Smoothing or Fitting
[smth_fit, xlinsp, AUC_smth, smthpnt] = BinSmoothing(dist_vec, score_vec);
% [smth_fit, xlinsp, AUC_smth, smthpnt] = BinSmoothing(dist_vec, score_vec, 'spline');
% [smth_fit, xlinsp, AUC_smth, smthpnt] = BinSmoothing(dist_vec, score_vec, 'pchip');
% [smth_fit, xlinsp, AUC_smth, smthpnt] = BinSmoothing(dist_vec, score_vec, 'makima');
% [smth_fit, xlinsp, AUC_smth, smthpnt] = BinSmoothing(dist_vec, score_vec, 'cubic');

[R2, slope, intercept] = Linearfitting(dist_vec, score_vec);
legend()
% titstr{mi} = { compose("Manif: pear %.3f spear %.3f",
end
end
toc
%%
[sorted_dist, sortid] = sort(dist_vec);
trapz(sorted_dist, score_vec(sortid))
%%
[smth_fit, xlinsp, AUC_smth, smthpnt] = BinSmoothing(dist_vec, score_vec);

function [AUC_raw] = trapz_raw(X,Y)
[sorted_X, sortid] = sort(X);
AUC_raw = trapz(sorted_X, Y(sortid));
end

function [gpr_fit, xlinsp, AUC, gprMdl] = GPRfitting(X, Y)
if nargin==2, varargin={};end
gprMdl = fitrgp(X, Y, 'SigmaLowerBound', 5E-3); % Edit this lower bound if encoutering fitting problem! 
xlinsp = linspace(0,max(X),100);
gpr_fit = gprMdl.predict(xlinsp');
AUC = trapz(xlinsp, gpr_fit);  % integrate under the gaussian process fitting curve. 
% Plot curve on a fig
plot(xlinsp, gpr_fit, 'DisplayName','Gaussian Process Regression', varargin{:})% , 'HandleVisibility','off'% adding these arguments will remove it from legend
end

function [smth_fit, xlinsp, AOC, smthcrv] = BinSmoothing(X, Y, interpmethod)
if nargin==2, varargin={}; interpmethod = 'linear'; 
elseif nargin==3, varargin={}; end
smthcrv = smooth(X,Y);
xlinsp = linspace(0,max(X),100);
[uniqX, iX, iUniq] = unique(X); % get rid of redundancy in X if there is any.
smth_fit = interp1(X(iX),smthcrv(iX),xlinsp,interpmethod,'extrap'); %'spline'
AOC = trapz(xlinsp, smth_fit);  % integrate under the gaussian process fitting curve. 
% Plot curve on a fig
scatter(X, smthcrv, 'DisplayName', 'Data Smooth')
plot(xlinsp, smth_fit,'Disp',strcat('smooth+interp ',interpmethod), varargin{:})%,'HandleVisibility','off' % adding these arguments will remove it from legend
end

function [r,m,b] = Linearfitting(X, Y, varargin) % Util function to add a linear reg line to a scatter
if nargin==2, varargin={};end
[r,m,b] = regression(reshape(X,1,[]), reshape(Y,1,[]));
% Plot line on a fig
xmin = min(X); xmax = max(X);
p = plot([xmin,xmax],[xmin,xmax].*m+b, 'DisplayName', 'Linear Fit',varargin{:});
% set(get(get(p(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end