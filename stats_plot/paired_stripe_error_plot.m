function paired_stripe_error_plot(var_cols, var_ci_cols, varnms, msks, labels,varargin)
Nrow = numel(var_cols{1});
Ncol = numel(var_cols);
assert(numel(var_cols) == numel(var_ci_cols))
if size(var_ci_cols{1},2)==2, has_CI = true; % has 2 side CI for each value 
elseif size(var_ci_cols{1},2)==1, has_CI = false; % symmetric, assume it's SEM 
end
varmat = [];
for i = 1:Ncol
    assert(numel(var_cols{i})==Nrow)
    varmat(:,i) = var_cols{i};
end
if isempty(msks), msks={ones(Nrow,1,'logical')}; labels="all";
else
for mi=1:numel(msks) % substitute empty mask
    if isempty(msks{mi}), msks{mi}=ones(Nrow,1,'logical');end
end
end
xjit = 0.15 * randn(Nrow,1);
mskofs = 0.8 / numel(msks);
offsets = ([1:numel(msks)] - 1 - (numel(msks)-1)/2) * mskofs;
figure;hold on;set(gcf,'pos',[1000         357         430         625])
for i=1:Ncol
for mi=1:numel(msks)
    offset = offsets(mi);
    yvalue = varmat(msks{mi},i);
    varmean = mean(yvalue);
    varsem = sem(yvalue);
    varN = numel(yvalue);
    displabel = compose("%s-%s %.1f+-%.1f(N=%d)",varnms(i),labels(mi),varmean,varsem,varN);
    if has_CI
        CI_col = var_ci_cols{i}(msks{mi},:);
        ypos = CI_col(:,2) - yvalue;
        yneg = CI_col(:,1) - yvalue;
    else
        CI_col = var_ci_cols{i}(msks{mi},:);
        ypos = CI_col;
        yneg = CI_col;
    end
    errorbar(i+offset+xjit(msks{mi}), yvalue, yneg, ypos, 'o',...
                'DisplayName',displabel,varargin{:});
end
end
for mi=1:numel(msks)
    plot(offsets(mi)+[1:Ncol]'+xjit(msks{mi})',varmat(msks{mi},:)','color',[0,0,0,0.1],'HandleVisibility','off')
end
newcolors = [[1, 0, 0];
             [0, 0, 1]];
colororder(newcolors)
xticks(1:Ncol); xticklabels(varnms)
legend('Location','best')
% title(T,compose("%s %s Pos Correlated Voxel Percent Across Layers",Animal,ExpType))
end