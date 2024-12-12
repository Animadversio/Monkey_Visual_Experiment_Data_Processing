function h = xscatter_error_plot(tab, statname1, statname2, stat_errname1, stat_errname2, ...
    masks, labels, titstr, savestr, sumdir, varargin)
% Scatter 2 variables in a table. 
% with special treatments for theta, phi and others. 
% marker = 'o*x^v';
% Corder = colororder;
% clrs = Corder([1,3,5,7],:);
if isempty(masks), masks={ones(size(tab,1),1,'logical')}; labels=["all"]; 
end
h = figure;clf;hold on;set(h,'pos',[1106         327         560         526]);
% statscol = cellfun(@(M)StatsTab_sum.(statname)(M), masks,'uni',0);
for i = 1:numel(masks)
    M = masks{i};
    statcol1 = tab.(statname1)(M);
    statcol2 = tab.(statname2)(M);
    mean_1 = mean(statcol1);
    sem_2 = sem(statcol1);
    mean_2 = mean(statcol2);
    sem_2 = sem(statcol2);
    N_i = sum(M);
    legstr = compose("%s %.3f(%.3f) %.3f(%.3f) (%d)",labels(i),mean_1,sem_2,mean_2,sem_2,N_i);
    if size(tab.(stat_errname1),2) == 2
        CI1 = tab.(stat_errname1)(M,:);
        xpos = CI1(:,2) - statcol1;
        xneg = CI1(:,1) - statcol1;
    elseif size(tab.(stat_errname1),2) == 1
        CI1 = tab.(stat_errname1)(M,:);
        xpos = CI1;
        xneg = CI1;
    end
    if size(tab.(stat_errname2),2) == 2
        CI2 = tab.(stat_errname2)(M,:);
        ypos = CI2(:,2) - statcol2;
        yneg = CI2(:,1) - statcol2;
    elseif size(tab.(stat_errname2),2) == 1
        CI2 = tab.(stat_errname2)(M,:);
        ypos = CI2;
        yneg = CI2;
    end
    errorbar(statcol1, statcol2,...
             yneg, ypos, xneg, xpos,'o',...
            'DisplayName',legstr,varargin{:});
    % scatter(StatsTab_sum.(statname1)(masks{i}),StatsTab_sum.(statname2)(masks{i}),...
    %     'DisplayName',legstr,varargin{:})
end
title_str = compose("Scatter of %s - %s for\n %s channels %s",...
    statname1,statname2,join(labels),titstr);
Ns = cellfun(@sum,masks); % Num of points within each musk
xlabel(statname1,'interpreter','none');ylabel(statname2,'interpreter','none'); %
title(title_str,'interpreter','none')
legend('Location','best')
% legend(compose("%s (%d)",labels',Ns')) % Num of points marked on label
if islogical(savestr) && (~savestr), 
return % if savestr==false, then no saving. 
else
saveallform(sumdir, compose("%s_%s_xscat_%s", statname1, statname2, savestr), h)
end
end