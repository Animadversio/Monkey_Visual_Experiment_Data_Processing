load('D:\Monkey_Data\monkey64chan-10062019-001_formatted.mat')

figure(1);imagesc(rasters(:,:,100));colorbar()
%% Create Color Table! 
fun = @(m)srgb_to_Lab(m);
color_seq = maxdistcolor(30,fun); % Use this to generate maximally distinguishable color sequence
color_seq = brewermap(101, 'spectral'); % Use this to generate gradual changing sequence
%% Natural Stimuli Response
cd 'NS_response'
for channel_j = 1:65 %1, 6, 30, 31 are bad silent channel 
figure('Position',[0,0,1000,600]);clf;hold on;
for i = 1:max(natural_stim_i)
    shadedErrorBar([],nat_stim_fr(i,:,channel_j),nat_stim_fr_sem(i,:,channel_j),...
    'lineprops',{'Color',[color_seq(i, :),0.6]},'transparent',1,'patchSaturation',0.075)
    % plot(nat_stim_fr(i,:,channel_j))
end
xlabel("time (ms)")
title(['Trial averaged PSTH of Natural Stimuli channel ', num2str(channel_j)])
saveas(gcf,sprintf("NS_rsp_channel%d.png",channel_j))
end
cd ..\
%% Evolved Image Response
mkdir 'Evolv_response'
cd 'Evolv_response'
for channel_j = 1:65 %1, 6, 30, 31 are bad silent channel 
h = figure(1);clf;hold on;
h.Position = [ 1          41        2560         963];
h.Visible = 'off'; 
gen_list = 0:100;
for i = 1:length(gen_list)
    shadedErrorBar([],evol_stim_fr(i, :, channel_j),evol_stim_sem(i, :, channel_j),...
    'lineprops',{'Color',[color_seq(i, :),0.85]},'transparent',1,'patchSaturation',0.075)
    % plot(evol_stim_fr(i,:,channel_j))
end
YL=ylim;YL(1)=0;ylim(YL);
XL=xlim;XL(1)=0;xlim(XL);
xlabel("time (ms)")
title(['Generation averaged PSTH of Evolved Stimuli channel ', num2str(channel_j)])
saveas(gcf,sprintf("Evolv_rsp_channel%d.png",channel_j))
hold off
end
cd ..\
%%
cbh = figure('Position', [500,400,1800,50]);
% for y = 1:ymx
%     num = raw(y).num;
%     typ = raw(y).typ;
%     map = raw(y).rgb(bmIndex(num,num,typ),:)/255; % downsample
y = 1; ymx=1.5; xmx = size(color_seq,1);
axh = axes('Parent',cbh, 'Color','none',...
	'XTick',0.5:10:xmx-0.5, 'YTick',0.5:ymx,'XTickLabel',{gen_list(1:10:end)},...
        'YTickLabel',{"Generation Color Code"},'YDir','reverse');
for x = 1:size(color_seq,1)
    patch([x-1,x-1,x,x],[y-1,y,y,y-1],1, 'FaceColor',color_seq(x,:), 'Parent',axh)
end
xlim([0, xmx])
axh.YAxis.FontSize = 14;
axh.XAxis.FontSize = 10;
set(gcf,'Visible','on')
% text(xmx+0.1,y-0.5,typ, 'Parent',axh, 'FontName',axf)
% end
saveas(gcf, ".\Evolv_response\Generation_Color_code.png")


%%
% `gen_num_i` maps from sorted id to generation number, map natural images
% as -1
trial_id_mask = sort_idx(gen_num_i~=-1);
part_gen_num = gen_num_i(gen_num_i~=-1);
for channel_j = 16:65
cluster_input = rasters(trial_id_mask, :, channel_j);
if sum(cluster_input,'all')==0
   continue 
end
Z = linkage(cluster_input, 'average', 'euclidean');%''correlation
UL = prctile(cluster_input(:), 98)+1;
LL =  prctile(cluster_input(:), 2)-1;
%%
figure(17);clf;
set(gcf, "Position",[0,40,2560,960])
suptitle(sprintf("Evolving Image PSTH Sorted by Hierachical clustering Channel %d",channel_j))
ax1 = subplot('Position', [0.06, 0.73, 0.88, 0.20]);
[~, T , perm_indx] = dendrogram(Z,0);
ax1.XTick=[];
ax2 = subplot('Position', [0.06, 0.15, 0.88, 0.55]);
imagesc(cluster_input(perm_indx, :)',[0,UL])
colormap(ax2,'parula')
ch1 = colorbar();
ch1.Position = [0.95, 0.15, 0.008, 0.55]; 
set(get(ch1, 'Label'), 'string','Firing rate','Fontsize',12);
ylabel("Time (ms)   ",'Fontsize',14)
ax2.YLabel.Rotation = 0;
ax2.XTick={};
ax3 = subplot('Position', [0.06, 0.08, 0.88, 0.05]);
sorted_gen_num = part_gen_num(perm_indx); 
imagesc(sorted_gen_num')
xlabel("Sorted Image id",'Fontsize',14)
ax3.YTick = [1];
ax3.YTickLabel = {"Generation #"};
ax3.YAxis.FontSize = 14;
colormap(ax3,color_seq)
caxis([min(gen_list),max(gen_list)]);
ch2 = colorbar();
ch2.Position = [0.95, 0.08, 0.008, 0.05]; 
set(get(ch2, 'Label'), 'string','Generation #','Fontsize',12);

saveas(gcf,sprintf(".\\Evolv_Sort_Rsp\\Evolv_rsp_sort_channel%d.png",channel_j))
end
%%
trial_id_mask = sort_idx(gen_num_i~=-1);
part_gen_num = gen_num_i(gen_num_i~=-1);
for channel_j = 1:65
cluster_input = rasters(trial_id_mask, :, channel_j);
if sum(cluster_input,'all')==0
   continue 
end
Z = linkage(cluster_input, 'average', 'correlation');%''
UL = prctile(cluster_input(:), 98)+1;
%
figure(18);clf;
set(gcf, "Position",[0,40,2560,960])
suptitle(sprintf("Evolving Image PSTH Sorted by Hierachical clustering (correlation) Channel %d",channel_j))
ax1 = subplot('Position', [0.06, 0.73, 0.88, 0.20]);
[~, T , perm_indx] = dendrogram(Z,0);
ax1.XTick=[];
ax2 = subplot('Position', [0.06, 0.15, 0.88, 0.55]);
imagesc(cluster_input(perm_indx, :)',[0,UL])
colormap(ax2,'parula')
ch1 = colorbar();
ch1.Position = [0.95, 0.15, 0.008, 0.55]; 
set(get(ch1, 'Label'), 'string','Firing rate','Fontsize',12);
ylabel("Time (ms)   ",'Fontsize',14)
ax2.YLabel.Rotation = 0;
ax2.XTick={};
ax3 = subplot('Position', [0.06, 0.08, 0.88, 0.05]);
sorted_gen_num = part_gen_num(perm_indx); 
imagesc(sorted_gen_num')
xlabel("Sorted Image id",'Fontsize',14)
ax3.YTick = [1];
ax3.YTickLabel = {"Generation #"};
ax3.YAxis.FontSize = 14;
colormap(ax3,color_seq)
caxis([min(gen_list),max(gen_list)]);
ch2 = colorbar();
ch2.Position = [0.95, 0.08, 0.008, 0.05]; 
set(get(ch2, 'Label'), 'string','Generation #','Fontsize',12);

saveas(gcf,sprintf(".\\Evolv_Corr_Sort_Rsp\\Evolv_rsp_sort_channel%d.png",channel_j))
end

%%
figure("Position",[0,0,1500,500])
imagesc(cluster_input(:, :)',[0,500])
colorbar()
title(sprintf("Evolving Image PSTH Channel %d",channel_j))
ylabel("Time")
xlabel("Image id")
