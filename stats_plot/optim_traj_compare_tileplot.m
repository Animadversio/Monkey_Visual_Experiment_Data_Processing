function figh = optim_traj_compare_tileplot(score_m_traj_extrap_col, block_traj_extrap_col, extrap_mask_col, ...
	mskcol, labcol, mskrow, labrow, thread_labels, plot_indiv, Xlim)
if nargin<10, Xlim = true;end
expN = numel(block_traj_extrap_col);
if isempty(mskcol), mskcol={ones(expN,1,'logical')}; labcol="";
else, mskcol = cellfun(@(msk)reshape(msk,[],1),mskcol,'uni',0); end
if isempty(mskrow), mskrow={ones(expN,1,'logical')}; labrow="";
else, mskrow = cellfun(@(msk)reshape(msk,[],1),mskrow,'uni',0); end

Corder = colororder();
blockN_arr = cellfun(@(M)sum(~M), extrap_mask_col); 
ncol = numel(mskcol);
nrow = numel(mskrow);
figh = figure('position', [200,100,300*ncol,300*nrow]);
T = tiledlayout(nrow,ncol,'TileSp','compact','Pad','compact');
for ci = 1:ncol
	for ri = 1:nrow
	ax = nexttile(T,ci+ncol*(ri-1));hold on
	msk = mskcol{ci} & mskrow{ri};
	[score_C_col_m,score_C_col_s,block_colvec] = sort_scoreblock(block_traj_extrap_col(msk),score_m_traj_extrap_col(msk,1));
	[score_G_col_m,score_G_col_s,block_colvec] = sort_scoreblock(block_traj_extrap_col(msk),score_m_traj_extrap_col(msk,2));
	% [score_C_col_m,score_C_col_s,block_colvec] = sort_scoreblock(cat(2,block_traj_extrap_col{msk}),cat(2,score_m_traj_extrap_col{msk,1}));
	% [score_G_col_m,score_G_col_s,block_colvec] = sort_scoreblock(cat(2,block_traj_extrap_col{msk}),cat(2,score_m_traj_extrap_col{msk,1}));
	shadedErrorBar(block_colvec,score_C_col_m,score_C_col_s,'lineProps',{'Color',Corder(2,:)})
	shadedErrorBar(block_colvec,score_G_col_m,score_G_col_s,'lineProps',{'Color',Corder(1,:)})
	title(ax,compose("%s %s (N=%d)",labrow(ri),labcol(ci),sum(msk)))
    if plot_indiv
        for iTr=1:expN
		    if ~msk(iTr),continue; end
		    blkmsk = extrap_mask_col{iTr};
		    plot(block_traj_extrap_col{iTr}(~blkmsk),movmean(score_m_traj_extrap_col{iTr,1}(~blkmsk),3),'Color',[Corder(2,:),0.3],'LineWidth',0.5)
		    plot(block_traj_extrap_col{iTr}(~blkmsk),movmean(score_m_traj_extrap_col{iTr,2}(~blkmsk),3),'Color',[Corder(1,:),0.3],'LineWidth',0.5)
		    % show the extrapolation part with blkmsk
			plot(block_traj_extrap_col{iTr}(blkmsk),movmean(score_m_traj_extrap_col{iTr,1}(blkmsk),3),'Color',[Corder(2,:),0.3],'LineWidth',0.5,'LineStyle','-.')
		    plot(block_traj_extrap_col{iTr}(blkmsk),movmean(score_m_traj_extrap_col{iTr,2}(blkmsk),3),'Color',[Corder(1,:),0.3],'LineWidth',0.5,'LineStyle','-.')
		end
    end
    if Xlim
%         disp(sort(blockN_arr(msk))')
		Xmax = prctile(blockN_arr(msk),85);
		xlim([0,Xmax])
    end
	end
end
nexttile(T,1);legend(thread_labels,'location','best')
title(T,"Trajectory Compare")
end