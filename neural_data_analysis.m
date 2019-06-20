load('D:\Monkey_Data\monkey64chan-10062019-001_formatted.mat')

figure(1);imagesc(rasters(:,:,100));colorbar()
%%
fun = @(m)srgb_to_Lab(m);
color_seq = maxdistcolor(30,fun); % Use this to generate maximally distinguishable color sequence

color_seq = brewermap(101, 'spectral'); % Use this to generate gradual changing sequence
%%
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
%%
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