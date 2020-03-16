clearvars -except meta_new rasters_new lfps_new Trials_new ExpSpecTable_Aug ExpRecord
%% Visualizing Code Evolution Traj

Animal = "Beto"; Set_Path;
result_dir = "C:\Users\ponce\OneDrive - Washington University in St. Louis\Optimizer_Tuning";
% load Generator
G = FC6Generator("matlabGANfc6");
%% Loading the Exp data
Triali = 4;
meta = meta_new{Triali};
rasters = rasters_new{Triali};
Trials = Trials_new{Triali};
exp_rowi = find(contains(ExpRecord.ephysFN, meta.ephysFN));
% Check the Expi match number
Expi = ExpRecord.Expi(exp_rowi);
fprintf("Processing  Exp %d:\n",Expi)
disp(ExpRecord.comments(exp_rowi))
% Fetch basic info of the experimebnt
pref_chan = Trials.TrialRecord.User.prefChan;
unit_in_pref_chan = cell2mat(Trials.TrialRecord.User.evoConfiguration(:,4))';
thread_num = size(Trials.TrialRecord.User.evoConfiguration, 1);
Exp_label_str = sprintf("Exp%d pref chan %d", Expi, pref_chan(1));
savepath = fullfile(result_dir, compose("%s_Evol%02d_chan%02d", Animal, Expi, pref_chan(1)));
Optim_names = [];
for i = 1:thread_num
    Optim_names = [Optim_names, strrep(string(Trials.TrialRecord.User.evoConfiguration{i,end}),"_"," ")];
end
unit_name_arr = generate_unit_labels(meta.spikeID);
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);
for i = 1:thread_num
    pref_chan_id(i) = find(meta.spikeID==pref_chan(i) & ... % the id in the raster and lfps matrix 
                    unit_num_arr==unit_in_pref_chan(i)); % match for unit number
end

imgnm = Trials.imageName;
% seperate the thread natural images and generated images 
row_gen = contains(imgnm, "gen") & ... % contains gen
        contains(imgnm, "block") & ... % contains block 
        cellfun(@(c) isempty(regexp(c(1:2), "\d\d")), imgnm) & ...% first 2 characters are not digits
        cellfun(@(c) ~contains(c(end-4:end), "_nat"), imgnm) ; % last 4 characters are not `_nat`
row_nat = ~row_gen;%contains(imgnm, "nat") & cellfun(@(c) ~isempty(regexp(c(1:2), "\d\d")), imgnm);
% seperate the thread 1 and 2 
row_thread0 = contains(imgnm, compose("thread%03d", 0));
row_thread1 = contains(imgnm, compose("thread%03d", 1));
assert(sum(row_thread0)+sum(row_thread1) == length(imgnm))
thread_msks = {row_thread0, row_thread1}; % store masks in a structure for the ease to iterate
% get the generation number 
block_arr = cell2mat(Trials.block);

% Sort image name and get scores in the given window
block_list = min(block_arr):max(block_arr);
scores_tsr = squeeze(mean(rasters(:, 51:200, :), 2) - mean(rasters(:, 1:40, :), 2)); % [unit_num, img_nums]
meanscore_syn = nan(size(rasters, 1), length(block_list), 2); % [unit_num, gen_nums, threads]
stdscore_syn = nan(size(rasters, 1), length(block_list), 2); 
meanscore_nat = nan(size(rasters, 1), length(block_list), 2);
stdscore_nat = nan(size(rasters, 1), length(block_list), 2);
for threadi = 1:thread_num 
    for blocki = min(block_arr):max(block_arr)
        gen_msk = row_gen & block_arr == blocki & thread_msks{threadi}; 
        nat_msk = row_nat & block_arr == blocki & thread_msks{threadi};
        meanscore_syn(:, blocki, threadi) = mean(scores_tsr(:, gen_msk), 2);
        meanscore_nat(:, blocki, threadi) = mean(scores_tsr(:, nat_msk), 2);
        stdscore_syn(:, blocki, threadi)  = std(scores_tsr(:, gen_msk), 1, 2) / sqrt(sum(gen_msk));
        stdscore_nat(:, blocki, threadi)  = std(scores_tsr(:, nat_msk), 1, 2) / sqrt(sum(nat_msk));
    end
end
%%
for threadi=1:2
channel_j = pref_chan_id(threadi);
scores_pref_ch = scores_tsr(pref_chan_id(threadi), :);
bsl_pref_ch = squeeze(mean(rasters(pref_chan_id(threadi), 1:40, :), 2));
rsp_pref_ch = squeeze(mean(rasters(pref_chan_id(threadi), 50:200, :), 2));
corr(bsl_pref_ch,max(0,scores_pref_ch') )
corr(bsl_pref_ch,scores_pref_ch') 
corr(bsl_pref_ch,rsp_pref_ch) 
% note the correlation between the baseline and score is significantly negative -0.8938! 
% Correlation bettwen rsponse and baseline rate is -0.1112, not
% significant. 
% The change of score is majorly (negatively) driven by fluctuation in baseline firing
% rate
%% 
% Find and load all codes 
[codes_all, img_ids, code_geni] = load_codes_all(meta.stimuli, threadi);
code_scores = nan(1, length(img_ids));
code_rsps = nan(1, length(img_ids));
% img_ids_wsfx = cellfun(@(c) c(1:end-4), img_ids, 'UniformOutput', false);
% For each generation in the experiment, sort the scores onto codes
for blocki = 2:max(block_arr)
    gen_msk = row_gen & block_arr == blocki & thread_msks{threadi}; 
    %nat_msk = row_nat & block_arr == blocki & thread_msks{threadi};
    gen_imgnms = imgnm(gen_msk);
    code_idx = zeros(1, length(gen_imgnms));
    for imgi = 1:numel(gen_imgnms)
        code_idx(imgi) = find(contains(img_ids, gen_imgnms{imgi}));
    end
    code_scores(code_idx) = scores_pref_ch(gen_msk);
    code_rsps(code_idx) = rsp_pref_ch(gen_msk);
end
%% get all codes in this gen, basis 
% do local PCA on codes? Or norm ? 
% Plot the image at the given location 
figure(4);clf;
colorvar = "gen";
for geni = 1:max(code_geni) - 1
    codes_basis = codes_all(code_geni <= geni+3 & code_geni >= geni-2, :);
    scores_basis = code_scores(code_geni <= geni+3 & code_geni >= geni-2); 
    codes_cur = codes_all(code_geni <= 45, :);
    scores_cur = code_scores(code_geni <= 45); 
    [~, RDmap] = pca(codes_basis, 3);
    codeRD = codes_cur * RDmap.M;
    sz_arr = 50*ones(1, length(scores_cur)); sz_arr(1:41:end) = 180;
    if colorvar == "score", clr_arr = code_scores; elseif colorvar == "gen", clr_arr = code_geni; end
    scatter3(codeRD(:,1),codeRD(:,2), codeRD(:,3), sz_arr, clr_arr, "filled")%,"MarkerFaceColor",0.5)
    xlabel("PC1");ylabel("PC2");zlabel("PC3")
    title([sprintf("Dimension reducing all code into PC123 of gen %d-%d", ...
        max(1, geni-2), min(max(code_geni), geni+3)),num2str(geni)])
    colorbar()
    pause
end

%% Visualize the gradient by interpolation / the next sample 
addpath ..\CMAES_optimizer_matlab\
figure(5);clf;set(5,'Position',[532         215        1028         763]);
figure(6);clf;set(6,'Position',[108         408        1894         284])
for geni = 1:max(code_geni) - 1
%     codes_basis = codes_all(code_geni <= geni+1 & code_geni >= geni, :);
%     scores_basis = code_scores(code_geni <= geni+1 & code_geni >= geni); 
%     [~, RDmap] = pca(codes_basis, 2);
    codes_cur = codes_all(code_geni == geni, :);
    scores_cur = code_rsps(code_geni == geni); % code_scores
    codes_next = codes_all(code_geni == geni + 1, :);
    scores_next = code_rsps(code_geni == geni + 1); % code_scores
    if ~contains(Optim_names{threadi},"CMAES")
    cur_basis = codes_cur(1,:);
    next_basis = codes_next(1,:); 
    [~, cur_idx] = sort(scores_cur(2:end),'descend');
    [~, next_idx] = sort(scores_next(2:end),'descend');
    ang_cur_samp = ang_dist(codes_cur(2:end, :),cur_basis);
    ang_next_samp = ang_dist(codes_next(2:end, :),next_basis);
    else
    weights = rankweight(40,20);
    [~, cur_idx] = sort(scores_cur,'descend');
    [~, next_idx] = sort(scores_next,'descend');
    cur_basis = weights * codes_cur(cur_idx,:);
    next_basis = weights * codes_next(next_idx,:);
    ang_cur_samp = ang_dist(codes_cur,cur_basis);
    ang_next_samp = ang_dist(codes_next,next_basis);
    end
    ang_among_cur = real(acos(1-pdist(codes_cur,"cosine")));%pdist(codes_cur)
    ang_among_next = real(acos(1-pdist(codes_next,"cosine")));
    set(0,"CurrentFigure",5)
    subplot(211)
    if geni == 1, imgs_cur = G.visualize(codes_cur(end-39:end,:));else, imgs_cur = imgs_next; end
    montage(imgs_cur(:,:,:,cur_idx),'Size',[4,10]) % sort the images from the most exciting ones to least
    title(sprintf("Gen %02d: Mean %.1f (%.2f) Max %.1f\n code dist %.3f (%.1f deg), ang dist to basis %.3f (%.1f deg)", ...
        geni+1, mean(scores_cur), std(scores_cur), max(scores_cur), mean(ang_among_cur), mean(ang_among_cur) / pi *180,...
        mean(ang_cur_samp), mean(ang_cur_samp) / pi * 180))
%     subplot(224)
    subplot(212)
    imgs_next = G.visualize(codes_next(end-39:end,:));
    montage(imgs_next(:,:,:,next_idx),'Size',[4,10])
    title(sprintf("Gen %02d: Mean %.1f (%.2f) Max %.1f\n code dist %.3f (%.1f deg), ang dist to basis %.3f (%.1f deg)", ...
        geni+2, mean(scores_next), std(scores_next), max(scores_next), mean(ang_among_next), mean(ang_among_next) / pi *180, ...
        mean(ang_next_samp), mean(ang_next_samp) / pi *180))
    
    interp_codes = sphere_interp(cur_basis, next_basis, linspace(-0.5,2, 6));
    interp_imgs = G.visualize(interp_codes);
    set(0,"CurrentFigure",6);clf%figure(6);
    montage(interp_imgs,'Size',[1,size(interp_imgs,4)])
    set(6,'Position', [1352         640        1107         236])
    ang_basis = ang_dist(cur_basis,next_basis);
    title(sprintf("Gen %02d: From %.1f to %.1f; ang dist %.3f (%.1f deg)", ...
        geni+2, scores_cur(1), scores_next(1), ang_basis, ang_basis / pi * 180))
%     subplot(221)
%     imshow(G.visualize(cur_basis))
%     title(num2str(scores_cur(1),"%.2f"))
%     subplot(223)
%     imshow(G.visualize(next_basis))
%     title(num2str(scores_next(1),"%.2f"))
%     
%     codeRD = codes_cur * RDmap.M;
%     sz_arr = 50*ones(1, length(scores_cur)); sz_arr(1:41:end) = 150;
%     scatter3(codeRD(:,1),codeRD(:,2), codeRD(:,3), sz_arr, scores_cur,"filled")%,"MarkerFaceColor",0.5)
%     xlabel("PC1");ylabel("PC2");zlabel("PC3")
%     colorbar()
     keyboard
end
%% Is the new samples closer to the mean than 
% Local geometric analysis 
figure(4);clf
geom_char = repmat(struct("ang_basis_step",[],"ang_among", [], "ang_to_basis", [], "ang_to_next", []), max(code_geni), 1);
for geni = 2:max(code_geni)-1
    codes_last = codes_all(code_geni == geni - 1, :);
    scores_last = code_scores(code_geni == geni - 1); % code_scores
    codes_cur = codes_all(code_geni == geni, :);
    scores_cur = code_scores(code_geni == geni); % code_scores
    codes_next = codes_all(code_geni == geni + 1, :);
    scores_next = code_scores(code_geni == geni + 1); % code_scores
    if contains(Optim_names{threadi},"ZOHA")
    % For ZOHA, the first code is the basis
    cur_basis = codes_cur(1,:); % basis for current generation of codes
    next_basis = codes_next(1,:); 
    [~, cur_idx] = sort(scores_cur(2:end),'descend');
    [~, next_idx] = sort(scores_next(2:end),'descend');
    ang_basis_cur_samp = ang_dist(codes_cur(2:end, :), cur_basis);
    ang_cur_samp_next_basis = ang_dist(codes_cur(2:end, :), next_basis);
    ang_among_cur = real(acos(1-pdist(codes_cur(2:end, :), "cosine")));
    ang_cur_samp_next_basis = ang_cur_samp_next_basis(cur_idx);
    elseif contains(Optim_names{threadi},"CMAES")
    % For ZOHA, the basis is implicit, has to be computed from linearly
    % recombine the codes
    weights = rankweight(40,20);
    [~, last_idx] = sort(scores_last,'descend');
    [~, cur_idx] = sort(scores_cur,'descend');
    [~, next_idx] = sort(scores_next,'descend');
    cur_basis = weights * codes_last(last_idx, :); % basis for current generation of codes
    next_basis = weights * codes_cur(cur_idx, :); 
    ang_basis_cur_samp = ang_dist(codes_cur, cur_basis);
    ang_cur_samp_next_basis = ang_dist(codes_cur, next_basis);
    ang_among_cur = real(acos(1-pdist(codes_cur,"cosine")));
    ang_cur_samp_next_basis = ang_cur_samp_next_basis(cur_idx);
    end
    % Compute some geometric charaters
    geom_char(geni).ang_basis_step = ang_dist(cur_basis, next_basis);
    geom_char(geni).ang_among = ang_among_cur;
    geom_char(geni).ang_to_basis = ang_basis_cur_samp;
    geom_char(geni).ang_to_next = ang_cur_samp_next_basis;
    set(0,"CurrentFigure",4);hold on %clf;
    plot(ang_cur_samp_next_basis,'linew',1.5,"color",[0,0,1,0.3])%"blue",'edgealpha',0.4)
    ylim([0.2, 1.5])
    xlabel("current samples number (sorted by score)") % from next basis to 
    ylabel("ang dist")
    title(sprintf("Gen %d %s\n red hline, ang to basis %.3f\n black hline, ang among samps %.3f\n green hline, step ang between basis %.3f\n blue ang from next basis to cur samps",...
        geni, Optim_names{threadi},mean(geom_char(geni).ang_to_basis), mean(geom_char(geni).ang_among), geom_char(geni).ang_basis_step))
%     hL = fig_horzline(mean(geom_char(geni).ang_among), axis);set(hL,'color','k','linew',1.5,'lines','-')
%     hL = fig_horzline(mean(geom_char(geni).ang_among) + std(geom_char(geni).ang_among), axis);set(hL,'color','k','linew',1.5,'lines',':')
%     hL = fig_horzline(mean(geom_char(geni).ang_among) - std(geom_char(geni).ang_among), axis);set(hL,'color','k','linew',1.5,'lines',':')
%     hL = fig_horzline(mean(geom_char(geni).ang_to_basis), axis);set(hL,'color','r','linew',1.5,'lines','-')
%     hL = fig_horzline(mean(geom_char(geni).ang_to_basis) + std(geom_char(geni).ang_to_basis), axis);set(hL,'color','r','linew',1.5,'lines',':')
%     hL = fig_horzline(mean(geom_char(geni).ang_to_basis) - std(geom_char(geni).ang_to_basis), axis);set(hL,'color','r','linew',1.5,'lines',':')
    hL = fig_horzline(geom_char(geni).ang_basis_step, axis);set(hL,'color',[0,1,0,0.4],'lines','-')
    % pause
end
%%
figure(2);clf;

%%
geom_char_ZOHA =geom_char;
%%
figure(2);clf;
geom_cell = struct2cell(geom_char);
shadedErrorBar(1:max(code_geni), cellfun(@mean, geom_cell(1,:)), cellfun(@std,  geom_cell(1,:)),... % step between basis
    'lineprops',{'Color',[0,1,0,0.7]},'transparent',1,'patchSaturation',0.075)
shadedErrorBar(1:max(code_geni), cellfun(@mean, geom_cell(2,:)), cellfun(@std,  geom_cell(2,:)),... % ang among samples
    'lineprops',{'Color',[0,0,0,0.7]},'transparent',1,'patchSaturation',0.075) 
shadedErrorBar(1:max(code_geni), cellfun(@mean, geom_cell(3,:)), cellfun(@std,  geom_cell(3,:)),... % cur samp to cur basis
    'lineprops',{'Color',[1,0,0,0.7]},'transparent',1,'patchSaturation',0.075)
hold on
plot(1:max(code_geni), cellfun(@mean, geom_cell(4,:)), "-.", 'color', [0,0,1,0.5]);hold on
cellfun(@nanmin, geom_cell(4,2:end-1))
plot(2:max(code_geni)-1, cellfun(@nanmin, geom_cell(4,2:end-1)), 'color', [0,0,1,0.8])
% shadedErrorBar(1:max(code_geni), cellfun(@mean, geom_cell(4,:)), cellfun(@std,  geom_cell(4,:)),... % cur samp to next basis
%     'lineprops',{'Color',[0,0,1,0.7]},'transparent',1,'patchSaturation',0.075)
title([sprintf("%s", Optim_names{threadi}),sprintf("red line, ang to basis \n black line, ang among samps \n green line, step ang between basis \n blue line, ang from next basis to cur samps (solid is min, -. is mean)")])
xlabel("generations")
ylabel("angular distance")

%%
cellfun(@mean, geom_cell(2,:))
cellfun(@std,  geom_cell(2,:))
cellfun(@mean, geom_cell(3,:))
cellfun(@std,  geom_cell(3,:))
cellfun(@mean, geom_cell(4,:))
cellfun(@std,  geom_cell(4,:))
%%
geom_char_CMAES = geom_char;
%%
figure(2);clf
plot(ang_cur_samp_next_basis,'linew',1.5)
xlabel("current samples number (sorted by score)") % from next basis to 
ylabel("ang dist")
title(sprintf("Gen %d\n red hline, ang to basis %.3f\n black hline, ang among samps %.3f\n green hline, step ang between basis %.3f\n blue ang from next basis to cur samps",...
    geni, mean(geom_char(geni).ang_to_basis), mean(geom_char(geni).ang_among), geom_char(geni).ang_basis_step))
hL = fig_horzline(mean(geom_char(geni).ang_among), axis);set(hL,'color','k','linew',1.5,'lines','-')
hL = fig_horzline(mean(geom_char(geni).ang_among) + std(geom_char(geni).ang_among), axis);set(hL,'color','k','linew',1.5,'lines',':')
hL = fig_horzline(mean(geom_char(geni).ang_among) - std(geom_char(geni).ang_among), axis);set(hL,'color','k','linew',1.5,'lines',':')
hL = fig_horzline(mean(geom_char(geni).ang_to_basis), axis);set(hL,'color','r','linew',1.5,'lines','-')
hL = fig_horzline(mean(geom_char(geni).ang_to_basis) + std(geom_char(geni).ang_to_basis), axis);set(hL,'color','r','linew',1.5,'lines',':')
hL = fig_horzline(mean(geom_char(geni).ang_to_basis) - std(geom_char(geni).ang_to_basis), axis);set(hL,'color','r','linew',1.5,'lines',':')
hL = fig_horzline(geom_char(geni).ang_basis_step, axis);set(hL,'color','green','lines','-')
%%
ang_dist(codes_next, cur_basis)
%%
interp_codes = sphere_interp(cur_basis, next_basis, 9);
interp_imgs = G.visualize(interp_codes);
set(0,"CurrentFigure",6);clf%figure(6);
montage(interp_imgs,'Size',[1,9])
title(sprintf("Gen %02d: From %.1f to %.1f; ang dist %.f", ...
        geni+2, scores_cur(1), scores_next(1), ang_dist(cur_basis,next_basis)))
%%
figure,plot(code_scores);hold on ;plot(code_rsps);xlabel("code num");ylabel("firing rate")
title(["Optim Tuning",Exp_label_str,sprintf("Chan %s", unit_name_arr(channel_j))])
%%
figure,scatter(code_geni, code_scores,"filled");hold on ;scatter(code_geni, code_rsps,"filled");xlabel("generations");ylabel("firing rate")
title(["Optim Tuning",Exp_label_str,sprintf("Chan %s", unit_name_arr(channel_j))])
%%
function ang = ang_dist(V1,V2)
nV1 = V1 ./ norm_axis(V1,2);
nV2 = V2 ./ norm_axis(V2,2);
cosang = nV1 * nV2';
ang = real(acos(cosang));
end
function norms = norm_axis(Vecs, dim)
norms = sqrt(sum(Vecs.^2,dim));
end
function Vecs = sphere_interp(V1, V2, samp_n)% So called SLERP
nV1 = V1 / norm(V1);
nV2 = V2 / norm(V2);
cosang = nV1 * nV2';
ang = acos(cosang);
if size(samp_n) == 1
    lingrid = linspace(0,1,samp_n)';
else
    lingrid = reshape(samp_n, [], 1);
end
Vecs = (sin((1-lingrid) * ang) * V1 + sin(lingrid * ang) * V2) / sin(ang);
end