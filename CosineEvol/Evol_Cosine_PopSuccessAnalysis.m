% Analyze the Cosine objective increase / improve 
% 
matdir = "E:\OneDrive - Harvard University\Mat_Statistics";
load(fullfile(matdir, "Alfa_CosStats.mat"))
%%
figdir1 = "E:\OneDrive - Harvard University\PhDDefense_Talk\Figures\Cosine_Figure";
figdir2 = "E:\OneDrive - Harvard University\CosinePoster_GermanConf\Figures";
%%
figure('Position',[366          97        1260         900]);
T=tiledlayout("flow",'TileSpacing','tight','Padding','tight');
for Expi = 1:numel(CosStats)
nexttile
nblocks = numel(CosStats(Expi).scores.ol_gen_avg);
shadedErrorBar(1:nblocks,CosStats(Expi).scores.ol_gen_avg,CosStats(Expi).scores.ol_gen_sem,'lineProps','-r');
shadedErrorBar(1:nblocks,CosStats(Expi).scores.ol_nat_avg,CosStats(Expi).scores.ol_nat_sem,'lineProps','-g');
title(compose("Exp%02d %s",Expi,CosStats(Expi).evol.score_mode),'Interpreter','none')
if Expi == 1
   legend(["gen","nat"])
end
end
title(T, "Cosine Evolution, Online objective comparison [Monkey A]")
saveallform([figdir1,figdir2],"online_score_traj_montage",gcf)
%%
score_mode_col = arrayfun(@(Expi) CosStats(Expi).evol.score_mode,1:numel(CosStats),'Uni',0);
online_obj_gen_init = arrayfun(@(Expi) CosStats(Expi).scores.ol_gen_avg(1),1:numel(CosStats));
online_obj_gen_end = arrayfun(@(Expi) CosStats(Expi).scores.ol_gen_avg(end),1:numel(CosStats));
online_obj_gen_init_sem = arrayfun(@(Expi) CosStats(Expi).scores.ol_gen_sem(1),1:numel(CosStats));
online_obj_gen_end_sem = arrayfun(@(Expi) CosStats(Expi).scores.ol_gen_sem(end),1:numel(CosStats));
online_obj_nat_init = arrayfun(@(Expi) CosStats(Expi).scores.ol_nat_avg(1),1:numel(CosStats));
online_obj_nat_end = arrayfun(@(Expi) CosStats(Expi).scores.ol_nat_avg(end),1:numel(CosStats));
online_obj_nat_init_sem = arrayfun(@(Expi) CosStats(Expi).scores.ol_nat_sem(1),1:numel(CosStats));
online_obj_nat_end_sem = arrayfun(@(Expi) CosStats(Expi).scores.ol_nat_sem(end),1:numel(CosStats));
%%
dotmsk = contains(score_mode_col,"dot");
corrmsk = contains(score_mode_col,"corr");
MSEmsk = contains(score_mode_col,"MSE");
figure;set(gcf,'Position',[680   414   840   560]);
T=tiledlayout(1,3,"TileSpacing","tight","Padding","compact");
% paired_stripe_error_plot(online_obj_gen_init,online_obj_gen_end,online_obj_gen_init_sem,online_obj_gen_end_sem)
ax=nexttile(1);msk=dotmsk;
paired_strip_plot(online_obj_gen_init(msk),online_obj_gen_end(msk),online_obj_gen_init_sem(msk),online_obj_gen_end_sem(msk))
xticklabels(["init","end"]);title("Dot Product",'fontsize',16);
ax.XAxis.FontSize = 14;ax.YAxis.FontSize = 14;
ax=nexttile(2);msk=corrmsk;
paired_strip_plot(online_obj_gen_init(msk),online_obj_gen_end(msk),online_obj_gen_init_sem(msk),online_obj_gen_end_sem(msk))
xticklabels(["init","end"]);title("Pearson Correlation",'fontsize',16);
ax.XAxis.FontSize = 14;ax.YAxis.FontSize = 14;
ax=nexttile(3);msk=MSEmsk;
paired_strip_plot(online_obj_gen_init(msk),online_obj_gen_end(msk),online_obj_gen_init_sem(msk),online_obj_gen_end_sem(msk))
xticklabels(["init","end"]);title("MSE",'fontsize',16);
ax.XAxis.FontSize = 14;ax.YAxis.FontSize = 14;
title(T,"Cosine Evolution, Objective Comparison [Monkey A]")
saveallform([figdir1,figdir2],"cosine_evol_obj_comparison",gcf)

%%
figure;set(gcf,'Position',[680   414   840   560]);
T=tiledlayout(1,3,"TileSpacing","tight","Padding","compact");
% paired_stripe_error_plot(online_obj_gen_init,online_obj_gen_end,online_obj_gen_init_sem,online_obj_gen_end_sem)
ax=nexttile(1);msk=dotmsk;
paired_strip_plot(online_obj_nat_init(msk),online_obj_nat_end(msk),online_obj_nat_init_sem(msk),online_obj_nat_end_sem(msk),{'k','k'});hold on;
paired_strip_plot(online_obj_gen_init(msk),online_obj_gen_end(msk),online_obj_gen_init_sem(msk),online_obj_gen_end_sem(msk))
xticklabels(["init","end"]);title("Dot Product",'fontsize',16);
ax.XAxis.FontSize = 14;ax.YAxis.FontSize = 14;
ax=nexttile(2);msk=corrmsk;
paired_strip_plot(online_obj_nat_init(msk),online_obj_nat_end(msk),online_obj_nat_init_sem(msk),online_obj_nat_end_sem(msk),{'k','k'});hold on;
paired_strip_plot(online_obj_gen_init(msk),online_obj_gen_end(msk),online_obj_gen_init_sem(msk),online_obj_gen_end_sem(msk))
xticklabels(["init","end"]);title("Pearson Correlation",'fontsize',16);
ax.XAxis.FontSize = 14;ax.YAxis.FontSize = 14;
ax=nexttile(3);msk=MSEmsk;
paired_strip_plot(online_obj_nat_init(msk),online_obj_nat_end(msk),online_obj_nat_init_sem(msk),online_obj_nat_end_sem(msk),{'k','k'});hold on;
paired_strip_plot(online_obj_gen_init(msk),online_obj_gen_end(msk),online_obj_gen_init_sem(msk),online_obj_gen_end_sem(msk))
xticklabels(["init","end"]);title("MSE",'fontsize',16);
ax.XAxis.FontSize = 14;ax.YAxis.FontSize = 14;
title(T,"Cosine Evolution, Objective Comparison [Monkey A]")
saveallform([figdir1,figdir2],"cosine_evol_obj_ctrl_comparison",gcf)


%%
function paired_strip_plot(list1, list2, sem1, sem2, clrs)
    if nargin == 4
        clrs = {'b','r'};
    end
    % Check if the input lists are of the same length
    if length(list1) ~= length(list2) || length(list1) ~= length(sem1) || length(list2) ~= length(sem2)
        error('All input lists must have the same length.');
    end
    
    % Add a jitter to the x-values
    jitter_amount = 0.35; % You can adjust the jitter amount as needed
    jitter1 = jitter_amount * rand(size(list1)) - jitter_amount/2;
    % jitter2 = jitter_amount * rand(size(list2)) - jitter_amount/2;
    
    % Plot the first list
    scatter(ones(size(list1)) + jitter1, list1, clrs{1}, 'filled');
    hold on;

    % Plot the error bars for the first list
    errorbar(ones(size(list1)) + jitter1, list1, sem1, clrs{1}, 'LineStyle', 'none');
    
    % Plot the second list
    scatter(2 * ones(size(list2)) + jitter1, list2, clrs{2}, 'filled');
    
    % Plot the error bars for the second list
    errorbar(2 * ones(size(list2)) + jitter1, list2, sem2, clrs{2}, 'LineStyle', 'none');
    
    % Connect the paired points with gray lines
    for i = 1:length(list1)
        plot([1 + jitter1(i), 2 + jitter1(i)], [list1(i), list2(i)], 'Color', [0.5, 0.5, 0.5]);
    end
    
    % Add labels and title
    set(gca, 'XTick', [1, 2], 'XTickLabel', {'List 1', 'List 2'});
    ylabel('Objective Values');
    % title('Paired Strip Plot');
    hold off;
end