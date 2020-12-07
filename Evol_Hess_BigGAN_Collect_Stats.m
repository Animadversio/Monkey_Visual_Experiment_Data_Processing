%% Find the Paired or unPaired Evolution that Accompany Hessian Experiments. 


Animal="Beto"; Set_Path; 
ftridx = find(contains(ExpRecord.Exp_collection,"BigGAN_Hessian") & contains(ExpRecord.expControlFN,"select") );
% Code to procedurely find the corresponding Evolution code. 
%   Basically find the BigGAN evolution that preceeds the hessian with same size
%   and channel. 
ftridx_evol = [];
for idx = ftridx'
    for dif=1:5
    try
        assert(contains(ExpRecord.expControlFN(idx-dif),"generate_BigGAN"))
        assert(ExpRecord.pref_chan(idx-dif) == ExpRecord.pref_chan(idx))
        assert(ExpRecord.stim_size(idx-dif) == ExpRecord.stim_size(idx))
        break
    catch
        continue
    end
    end
    fprintf("BigGAN Evolution %s pref %d size %.1f\n", ExpRecord.expControlFN{idx-dif}, ExpRecord.pref_chan(idx-dif), ExpRecord.stim_size(idx-dif))
    fprintf("%s\n",ExpRecord.comments{idx-dif})
    fprintf("with Hessian Exp %s pref %d size %.1f\n", ExpRecord.expControlFN{idx}, ExpRecord.pref_chan(idx), ExpRecord.stim_size(idx))
    fprintf("%s\n\n",ExpRecord.comments{idx})
    ftridx_evol(end+1)=idx-dif;
end
%% Load the Spikes and Trials
[meta_new, rasters_new, lfps_new, Trials_new] = loadExperiments(ftridx_evol, Animal);
%%
HEStats = repmat(struct(),1,numel(ftridx));
D = torchImDist();
%%
%% Load the Images 
for Triali = 1:14
meta = meta_new{Triali};
rasters = rasters_new{Triali};
% lfps = lfps_new{Triali};
Trials = Trials_new{Triali};
if Triali == 5, meta.stimuli = "N:\Stimuli\2020-BigGAN\2020-08-13-Beto-01\2020-08-13-10-11-18";end
%% Record Basic Infos
Expi = Triali;% FIXME! Expi should be decoded from the ExpRecord.
HEStats(Expi).Animal = Animal;
HEStats(Expi).Expi = Expi; 
HEStats(Expi).imageName = Trials.imageName;
HEStats(Expi).meta = meta;

pref_chan = Trials.TrialRecord.User.prefChan;
unit_in_pref_chan = cell2mat(Trials.TrialRecord.User.evoConfiguration(:,4))';
imgpos = cell2mat(Trials.TrialRecord.User.evoConfiguration(:,2));
imgsize = cell2mat(Trials.TrialRecord.User.evoConfiguration(:,3));
thread_num = size(Trials.TrialRecord.User.evoConfiguration, 1);
unit_name_arr = generate_unit_labels(meta.spikeID);
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);
% note if the evolved channel is marked as null 0 then this can fail! 
pref_chan_id = zeros(1, thread_num);
for i = 1:thread_num % note here we use the unit_num_arr which exclude the null channels.
    pref_chan_id(i) = find(meta.spikeID==pref_chan(i) & ... % the id in the raster and lfps matrix 
                    unit_num_arr==unit_in_pref_chan(i)); % match for unit number
end
HEStats(Expi).units.pref_chan = pref_chan;
HEStats(Expi).units.unit_name_arr = unit_name_arr;
HEStats(Expi).units.unit_num_arr = unit_num_arr;
HEStats(Expi).units.activ_msk = activ_msk;
HEStats(Expi).units.spikeID = meta.spikeID;
HEStats(Expi).units.pref_chan_id = pref_chan_id;
HEStats(Expi).units.unit_in_pref_chan = unit_in_pref_chan;

if size(Trials.TrialRecord.User.evoConfiguration,2) == 5
Optim_names = [];
for i = 1:thread_num
    Optim_names = [Optim_names, strrep(string(Trials.TrialRecord.User.evoConfiguration{i,end}),"_"," ")]; % use this substitution for better printing effect
end
elseif size(Trials.TrialRecord.User.evoConfiguration,2) == 4 % old fashioned config
    Optim_names = repmat("CMAES", 1,thread_num);
end
HEStats(Expi).evol.optim_names = Optim_names;
HEStats(Expi).evol.thread_num = thread_num;
HEStats(Expi).evol.imgpos = imgpos;
HEStats(Expi).evol.imgsize = imgsize;
HEStats(Expi).evol.unit_in_pref_chan = unit_in_pref_chan; 
HEStats(Expi).evol.pref_chan = cell2mat(Trials.TrialRecord.User.evoConfiguration(:,1));
% Sort the image names
% seperate the thread natural images and generated images 
imgnm = Trials.imageName;
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
HEStats(Expi).stim.gen_msk = row_gen;
HEStats(Expi).stim.nat_msk = row_nat;
HEStats(Expi).stim.thread_msks = thread_msks;
%% 
block_arr = cell2mat(Trials.block);
block_list = min(block_arr):max(block_arr);
block_num = length(block_list);
color_seq = brewermap(block_num, 'spectral');
HEStats(Expi).color_seq = color_seq;
HEStats(Expi).evol.block_arr = block_arr;
HEStats(Expi).evol.block_n = block_num;
gen_idx_seq = cell(thread_num, block_num); % generated image idx cell as a horizontal array. 
nat_idx_seq = cell(thread_num, block_num);
for threadi = 1:thread_num 
    for blocki = min(block_arr):max(block_arr)
        gen_msk = row_gen & block_arr == blocki & thread_msks{threadi}; 
        nat_msk = row_nat & block_arr == blocki & thread_msks{threadi};
        gen_idx_seq{threadi, blocki} = find(gen_msk);
        nat_idx_seq{threadi, blocki} = find(nat_msk);
    end
end
gen_psth_col = cellfun(@(idx) rasters(pref_chan_id, :, idx), gen_idx_seq, 'Uni', 0);
nat_psth_col = cellfun(@(idx) rasters(pref_chan_id, :, idx), nat_idx_seq, 'Uni', 0);
gen_rspmat_col = cellfun(@(idx) squeeze(mean(rasters(:, 51:200, idx),[2])), gen_idx_seq, 'Uni', 0);
nat_rspmat_col = cellfun(@(idx) squeeze(mean(rasters(:, 51:200, idx),[2])), nat_idx_seq, 'Uni', 0);
HEStats(Expi).evol.idx_seq = gen_idx_seq;
HEStats(Expi).evol.psth = gen_psth_col;
HEStats(Expi).evol.rspmat = gen_rspmat_col;
HEStats(Expi).ref.idx_seq = nat_idx_seq;
HEStats(Expi).ref.psth = nat_psth_col;
HEStats(Expi).ref.rspmat = nat_rspmat_col;
% %% BigGAN images loading
% BGmsk = row_gen & thread_msks{2};
% BGimgnms = imgnm(BGmsk);
% tic
% BGimgs = arrayfun(@(nm)imread(fullfile(meta.stimuli, nm+".bmp")),BGimgnms,"Uni",0);
% %% imall = cat(4,BGimgs{:,:});
% tic
% distmat_evBG = D.distmat_B(cat(4,BGimgs{:,:})); % 240 sec for 950 images
% toc
% %
% HEStats(Expi).evol.distmat_BG = distmat_evBG;
end
%%
MatStatsDir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics";
savefast(fullfile(MatStatsDir, Animal+"_HessBGEvolStats.mat"), 'HEStats')
%%
distmat = distmat_evBG;
distmat = distmat + diag(nan(1,size(distmat,1)));
gen_rspmat = squeeze(cell2mat(reshape(gen_rspmat_col(2,:),1,1,[])));
for iCh = 1:size(gen_rspmat,1)
    rsp_vec = gen_rspmat(iCh,:);
%     [cc_vec, p_vec] = corrcoef_vec(rsp_vec', distmat);
    cc_vec = []; p_vec = [];
    for i = 1:length(distmat)
        [cc, P] = corr([rsp_vec',distmat(:,i)],'Rows','complete', 'Type', 'Spearman');
        cc_vec(i) = cc(1,2);
        p_vec(i) = P(1,2);
    end
    [cc_min, minidx] = min(cc_vec);
    fprintf("Chan%s %.3f %.1e\n", HEStats(Expi).units.unit_name_arr(iCh),cc_min,p_vec(minidx))
    h=figure(3); clf;
    scatter(distmat(:,minidx), rsp_vec')
    XLIM = xlim();xlinsp = [0:0.025:XLIM(2)];hold on
    gprMdl = fitrgp(distmat(:,minidx), rsp_vec');
    gpr_fit = gprMdl.predict(xlinsp');
    plot(xlinsp, gpr_fit)
    xlabel("LPIPS distance");ylabel("Firing Rate")
    title([compose("Evolution Exp %d (pref Ch%d) Chan%s %.3f %.1e", Expi, HEStats(Expi).units.pref_chan(2), ...
         HEStats(Expi).units.unit_name_arr(iCh),cc_min,p_vec(minidx))])
    pause
end
%% Load the codes and the find the match! 
[codes_BG, img_ids_BG, code_BG_geni] = load_codes_all(meta.stimuli, 2);
%  Find Matching Code for each trial. Leave nan there if no code is found (init gen)
code_idx = nan(numel(BGimgnms),1);
code4trial = nan(numel(BGimgnms), 128);
for i=1:numel(BGimgnms)
    foundidx = find(contains(img_ids_BG, BGimgnms(i)));
    if isempty(foundidx)
    code_idx(i) = nan;
    else
    code_idx(i) = foundidx;
    code4trial(i,:) = codes_BG(foundidx,:);
    end
end
%%
BGrspmat = squeeze(mean(rasters(:, 51:200, BGmsk),2));
BGpsthmat = arrayfun(@(beg) mean(rasters(:, beg:beg+10-1, BGmsk),2), 51:10:191, 'Uni', 0);
BGpsthmat = cat(2, BGpsthmat{:});
%% Regress the codes onto the PSTH mat
n2c_beta = mvregress(BGrspmat(:,31:end)',code4trial(31:end,:));
c2n_beta = mvregress(code4trial(31:end,:),BGrspmat(:,31:end)');
%%
nsamp = size(BGrspmat,2);
n2c_beta = pinv([ones(nsamp,1),BGrspmat']) * code4trial;
c2n_beta = pinv([ones(nsamp,1),code4trial]) * BGrspmat';
%%
[evc_all, eva_all, evc_cls, eva_cls, evc_nos, eva_nos] = loadHessian();
%%
fixnoise = Trials.TrialRecord.User.space_cfg{2}{2};
G = torchBigGAN();
G = G.select_space("class", fixnoise);
%%
imgidx = 20:40; 
figure(2);
subplot(121)
montage(G.visualize(code4trial(imgidx,:)))
subplot(122)
montage(G.visualize([ones(numel(imgidx), 1), BGrspmat(:,imgidx)'] * n2c_beta))