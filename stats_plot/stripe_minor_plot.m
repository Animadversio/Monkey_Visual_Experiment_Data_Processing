function h = stripe_minor_plot(tab, statname, masks, labels, minormsks, minorlabels, titstr, savestr, Tpairs, varargin)
% Parameter:
%   
% Example:
%   h = stripe_minor_plot(NPtab, "normAUS_bsl", {Fmsk&drvmsk,Fmsk&~drvmsk}, ["Driver","NonDriver"], {Alfamsk, Betomsk}, ["Alfa", "Beto"],...
%                    "all chan ANOVA P<1E-3", "drv_cmp_anim_sep", {[1,2]},'MarkerEdgeAlpha',0.3);

if nargin<7, Tpairs = {}; end
marker = 'o*x^v';
global figdir
h = figure;clf;hold on; set(h,'pos',[1686         323         369         531])
labelcol = []; Ns = [];
for i = 1:numel(masks)
    for j = 1:numel(minormsks)
    M = masks{i}&minormsks{j};
    statcol = tab.(statname)(M);
    mean_i = mean(statcol);
    sem_i = sem(statcol);
    N_i = sum(M);
    xjitter = 0.15*randn(N_i,1);
    legstr = compose("%s %.3f(%.3f) (%d)",labels(i)+minorlabels(j),mean_i,sem_i,N_i);
    scatter(i+xjitter,statcol,'Marker',marker(j),'DisplayName',legstr,varargin{:})
    labelcol = [labelcol, labels(i)+minorlabels(j)];
    Ns = [Ns, sum(M)];
    end
end
xticks(1:numel(masks));xticklabels(labels)
ylabel(statname,'interpreter','none')

statscol = cellfun(@(M)tab.(statname)(M), masks,'uni',0);
FStat = anova_cells(statscol);
title_str = compose("Comparison of %s for\n %s channels %s\nANOVA F:%.3f(p=%.1e df=%d)",...
    statname,join(labels),titstr,FStat.F,FStat.F_P,FStat.STATS.df);
for pi = 1:numel(Tpairs)
pair = Tpairs{pi};
[~,P,~,TSTAT]=ttest2(statscol{pair(1)},statscol{pair(2)});
title_str = title_str+compose("\n%s - %s: t=%.2f (p=%.1e (%d))",labels(pair(1)),labels(pair(2)),TSTAT.tstat,P,TSTAT.df);
end
title(title_str,'interpreter','none')
% Ns = cellfun(@sum,masks);
% legend(compose("%s(%d)",labelcol',Ns'),'Location','best') % ["Driver", "Non-Drivers"]
legend('Location','best')%
saveas(h,fullfile(figdir,compose("%s_%s_stripcmpmarker.png", statname, savestr)))
saveas(h,fullfile(figdir,compose("%s_%s_stripcmpmarker.pdf", statname, savestr)))
savefig(h,fullfile(figdir,compose("%s_%s_stripcmpmarker.fig", statname, savestr)))
end