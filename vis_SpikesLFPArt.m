Animal = "Alfa";Set_Path;
ftr = find(ExpRecord.Exp_collection=="Manifold"&ExpRecord.Expi==3);
no_return = false; no_lfp = false;no_conv= true;
[meta_new,rasters_new,lfps_new,Trials_new] = loadExperiments(ftr, Animal, no_return, no_lfp, no_conv);
%%
Triali = 1;
rasters = rasters_new{Triali};
lfps = lfps_new{Triali};
meta = meta_new{Triali};
Trials = Trials_new{Triali};

pref_chan = Trials.TrialRecord.User.prefChan;
unit_in_pref_chan = cell2mat(Trials.TrialRecord.User.evoConfiguration(:,4))';
thread_num = size(Trials.TrialRecord.User.evoConfiguration, 1);
if thread_num == 2, assert(pref_chan(1) == pref_chan(2)); end
% Exp_label_str = sprintf("Exp%d pref chan %d", Expi, pref_chan(1));
if contains(meta.ephysFN, "Beto"), Animal = "Beto"; elseif contains(meta.ephysFN, "Alfa"), Animal = "Alfa"; end
% savepath = fullfile(result_dir, compose("Manifold_%s_Evol%02d_chan%02d", Animal, Expi, pref_chan(1)));
% mkdir(savepath);
%%
imgnm = Trials.imageName;
% seperate the thread natural images and generated images 
row_gen = contains(imgnm, "gen") & ... % contains gen
        contains(imgnm, "block") & ... % contains block 
        cellfun(@(c) isempty(regexp(c(1:2), "\d\d")), imgnm) & ...% first 2 characters are not digits
        cellfun(@(c) ~contains(c(end-4:end), "_nat"), imgnm) ; % last 4 characters are not `_nat`
row_nat = ~row_gen;%contains(imgnm, "nat") & cellfun(@(c) ~isempty(regexp(c(1:2), "\d\d")), imgnm);
% seperate the thread 1 and 2 (and maybe thread 3 4)
thread_msks = cell(1, thread_num);
for threadi = 1:thread_num
    msk = contains(imgnm, compose("thread%03d", threadi - 1));
    thread_msks{threadi} = msk; % store masks in a structure for the ease to iterate
end
assert(sum(cellfun(@sum, thread_msks)) == length(imgnm)) % images comes from all these threads
block_arr = cell2mat(Trials.block);
block_num = max(block_arr);
%%
gen_idx_seq = cell(thread_num, block_num);
nat_idx_seq = cell(thread_num, block_num);
for threadi = 1:thread_num 
    for blocki = min(block_arr):max(block_arr)
        gen_msk = row_gen & block_arr == blocki & thread_msks{threadi}; 
        nat_msk = row_nat & block_arr == blocki & thread_msks{threadi};
        gen_idx_seq{threadi, blocki} = find(gen_msk);
        nat_idx_seq{threadi, blocki} = find(nat_msk);
    end
end
%%
iCh = find(meta.spikeID==pref_chan);
%%
evo_raster = squeeze(rasters(iCh(1),:,:));
%%
rasterplot(find(evo_raster(:,mask)),3525,200)
%%
mask = contains(Trials.imageName,"gen");
figure;plot(sum(evo_raster(51:200,mask)/1E3,1))
%%
figure;
for blocki = min(block_arr):max(block_arr)
    plot(squeeze(lfps(pref_chan,:,gen_idx_seq{1,blocki})))
    pause
end
%%
figure;
for iTr=2000:2100
clf;hold on
plot(squeeze(lfps(pref_chan,:,iTr)))
stem(squeeze(rasters(iCh(1),:,iTr)),'Marker','none')
stem(squeeze(rasters(iCh(2),:,iTr)),'Marker','none')
drawnow
pause
end
%%
bpFilt = designfilt('bandpassfir','FilterOrder',30, ...
         'CutoffFrequency1',3,'CutoffFrequency2',200, ...
         'SampleRate',1000);
fvtool(bpFilt)
     %%
% hpFilt = designfilt('highpassfir', 'StopbandFrequency',1,...
%          'PassbandFrequency',2,'PassbandRipple',0.5, ...
%          'StopbandAttenuation',20,'DesignMethod','kaiserwin','SampleRate',1000);

hpFilt = designfilt('highpassiir','FilterOrder',8, ...
         'PassbandFrequency',3,'PassbandRipple',0.2, ...
         'SampleRate',1000);
fvtool(hpFilt)
%%
lpFilt = designfilt('lowpassiir','FilterOrder',8, ...
         'PassbandFrequency',150,'PassbandRipple',0.2, ...
         'SampleRate',1000);
fvtool(lpFilt)
% lpFilt = designfilt('lowpassfir','PassbandFrequency',200, ...
%          'StopbandFrequency',300,'PassbandRipple',0.5, ...
%          'StopbandAttenuation',30,'DesignMethod','kaiserwin','SampleRate',1000);
% fvtool(lpFilt)
%%
% trace = lfps(pref_chan,:,iTr);
traces = squeeze(lfps(pref_chan,:,:));
lptraces = filtfilt(lpFilt,traces);
bptraces = filtfilt(hpFilt,lptraces);
%%
figure;hold on;
plot(trace)
plot(lptrace)
plot(bptrace)
%%
figure(9);
for blocki = min(block_arr):max(block_arr)
    clf;hold on 
    plot(bptraces(:,gen_idx_seq{1,blocki}),'color',[0.4,0.4,0.4,0.2])
    semtrace = sem(bptraces(:,gen_idx_seq{1,blocki}),2);
    meantrace = mean(bptraces(:,gen_idx_seq{1,blocki}),2);
%     plot(meantrace,'Color','red','LineWidth',2)
    shadedErrorBar([],meantrace,semtrace)
    ylim([-1000,1000])
    pause
end
%%
figure(13);clf;
for blocki = min(block_arr):max(block_arr)
    hold on 
    plot3(1:200,blocki*ones(1,200),bptraces(:,gen_idx_seq{1,blocki}),'color',[0.4,0.4,0.4,0.05],'LineWidth',0.4)
    sem(bptraces(:,gen_idx_seq{1,blocki}),2);
    meantrace = mean(bptraces(:,gen_idx_seq{1,blocki}),2);
    plot3(1:200,blocki*ones(1,200),meantrace,'LineWidth',3)
    pause
end
%%
figure(14);clf;
for blocki = min(block_arr):max(block_arr)
    hold on 
%     plot3(1:200,blocki*ones(1,200),bptraces(:,gen_idx_seq{1,blocki}),'color',[0.4,0.4,0.4,0.2],'LineWidth',0.4)
    sem(bptraces(:,gen_idx_seq{1,blocki}),2);
    meantrace = mean(bptraces(:,gen_idx_seq{1,blocki}),2);
    plot3(1:200,blocki*ones(1,200),meantrace,'LineWidth',3)
    pause
end
%%
power_per_gen = [];
var_per_gen = [];
for blocki = min(block_arr):max(block_arr)
    hold on 
    meantrace = mean(bptraces(:,gen_idx_seq{1,blocki}),2);
    semtrace = sem(bptraces(:,gen_idx_seq{1,blocki}),2);
    power_per_gen(blocki) = var(meantrace);
    var_per_gen(blocki) = mean(semtrace);
end
figure;plot(var_per_gen)
%%
figure(16);clf;
for blocki = min(block_arr):max(block_arr)
    hold on 
%     plot3(1:200,blocki*ones(1,200),bptraces(:,gen_idx_seq{1,blocki}),'color',[0.4,0.4,0.4,0.2],'LineWidth',0.4)
    sem(bptraces(:,gen_idx_seq{1,blocki}),2);
    meantrace = squeeze(mean(rasters(iCh(1),:,gen_idx_seq{1,blocki}),3));
    plot3(1:200,blocki*ones(1,200),meantrace,'LineWidth',3)
    pause
end
%%
figure(17);clf;
for blocki = min(block_arr):max(block_arr)
    hold on 
%     plot3(1:200,blocki*ones(1,200),bptraces(:,gen_idx_seq{1,blocki}),'color',[0.4,0.4,0.4,0.2],'LineWidth',0.4)
%     sem(bptraces(:,gen_idx_seq{1,blocki}),2);
    meantrace = squeeze(mean(rasters(iCh(1),:,gen_idx_seq{1,blocki}),3));
    sh(blocki) = stem3(1:200,blocki*ones(1,200),meantrace,'Marker','none','LineWidth',1,'Color',[0.3,0.3,0.3,0.05]);
%     cm_idx = max( conv(k1)-min(conv), 1 );
%     set(sh(blocki), 'Color',cm(cm_idx,:)
end
surf(1:200,min(block_arr):max(block_arr),mean_sdf_pergen,'FaceAlpha',0.9,'EdgeColor','none')
%%
ker = normpdf(-6:6,0,2);
raw_spikes = squeeze(rasters(iCh(1),:,:))';
raw_sdf = imfilter(raw_spikes,ker,'conv','replicate');
figure(18);clf;
for blocki = min(block_arr):max(block_arr)
    hold on 
    plot3(1:200,blocki*ones(1,200),raw_sdf(gen_idx_seq{1,blocki},:),'color',[0.4,0.4,0.4,0.02],'LineWidth',0.4)
    sem(bptraces(:,gen_idx_seq{1,blocki}),2);
    meantrace = squeeze(mean(raw_sdf(gen_idx_seq{1,blocki},:),1));
    plot3(1:200,blocki*ones(1,200),meantrace,'Marker','none','LineWidth',1.5)
    pause
end
%%
mean_sdf_pergen = arrayfun(@(blocki)squeeze(mean(raw_sdf(gen_idx_seq{1,blocki},:),1)),min(block_arr):max(block_arr),'Un',0);
mean_sdf_pergen = cell2mat(mean_sdf_pergen');
%%
figure(21);
surf(1:200,min(block_arr):max(block_arr),mean_sdf_pergen,'FaceAlpha',1,'EdgeColor','none')
