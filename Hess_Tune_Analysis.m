% Hessian Analaysis
% open myPaths
Animal = "Alfa";Set_Path;
% ftr = contains(ExpRecord.ephysFN, "25062020") ;%| contains(ExpRecord.ephysFN, "25062020");
ftr = contains(ExpRecord.Exp_collection, "Hessian");
rowlist = find(ftr);
[meta_new,rasters_new,lfps_new,Trials_new] = Project_Manifold_Beto_loadRaw(rowlist, Animal, false, true);
%%
Triali = 2;
meta = meta_new{Triali};
rasters = rasters_new{Triali};
% lfps = lfps_new{Triali};
Trials = Trials_new{Triali};

unit_name_arr = generate_unit_labels(meta.spikeID);
[activ_msk, unit_name_arr, unit_num_arr] = check_channel_active_label(unit_name_arr, meta.spikeID, rasters);

%%
imgname_uniq = unique(Trials.imageName); 
namepart_uniq = regexp(imgname_uniq,"norm(?<norm>\d*)_PC(?<pc_id>\d*)_ang(?<angle>[-\d]*)",'names');
namepart_uniq = namepart_uniq(~cellfun(@isempty,namepart_uniq));
pc_id_arr = unique(cellfun(@(U)str2num(U.pc_id),namepart_uniq));
angle_arr = unique(cellfun(@(U)str2num(U.angle),namepart_uniq));
sphere_norm = unique(cellfun(@(U)str2num(U.norm),namepart_uniq));

gen_msk = contains(Trials.imageName, "norm");
didGabor = false;
if sum(contains(imgname_uniq,"gab")) > 6
didGabor = true;
end
%% get the idx_arr
imgnm_arr = strings(length(pc_id_arr), length(angle_arr));
idx_arr = cell(length(pc_id_arr), length(angle_arr));
for i = 1:length(pc_id_arr)
    pc_id = pc_id_arr(i);
    for j = 1:length(angle_arr)
    ang = angle_arr(j);
    imgnm = compose("norm%d_PC%d_ang%d",sphere_norm,pc_id,ang);
    idx_arr{i,j} = find(contains(Trials.imageName, imgnm));
    imgnm_arr(i,j) = imgnm;
    end
end
%%
iCh = 5;
psth_arr = cellfun(@(idx)rasters(iCh, :, idx),idx_arr,'UniformOutput',false);
mean_psth_arr = cellfun(@(psth)mean(psth,3), psth_arr,'UniformOutput',false);
score_arr = cellfun(@(psth)mean(psth(:,51:200,:),[2,3]), psth_arr);
%
figure(1);
imagesc(score_arr);
axis image
title(compose("Alfa Channel %s Tuning",unit_name_arr(iCh)))
ylabel("Numbering of PC vector");xlabel("angle");
yticks(1:length(pc_id_arr));yticklabels(pc_id_arr)
xticks(1:2:11);xticklabels(angle_arr(1:2:11))
colorbar()
%%
wdw_vect = [1, 20] + 10 * [0:18]';% subsample to decrease redunancy
wdw_vect = [wdw_vect; [1,50]+[0:50:150]'; [51,200]];
for fi = 1:size(wdw_vect,1)
score_arr(:,:,fi) = cellfun(@(psth)mean(psth(:,wdw_vect(fi,1):wdw_vect(fi,2),:),[2,3]), psth_arr);
end
CLIM_arr = [ones(19,1)*prctile(score_arr(:,:,1:19),[1,99],'all')';...
            ones(4, 1)*prctile(score_arr(:,:,20:23),[1,99],'all')';...
            ones(1, 1)*prctile(score_arr(:,:,24),[1,99],'all')'];
    
%%
Save = true;
figure(3);set(3,'position',[1000         270         560         700])
IMS = imagesc(score_arr(:,:,1));
axis image
TIT = title(compose("Alfa Hessian Exp Channel %s Tuning",unit_name_arr(iCh)));
ylabel("Numbering of PC vector");xlabel("angle");
yticks(1:length(pc_id_arr));yticklabels(pc_id_arr)
xticks(1:2:11);xticklabels(angle_arr(1:2:11))
colorbar()
for fi = 1:size(wdw_vect,1)
    IMS.CData = score_arr(:,:,fi);
    TIT.String = compose("Alfa Hessian Exp Channel %s Tuning\nwindow [%d,%d] ms",unit_name_arr(iCh),wdw_vect(fi,1),wdw_vect(fi,2));
    caxis(CLIM_arr(fi,:));
    pause(0.5)
end

%
% figure(2); 
% plot(cell2mat(mean_psth_arr(1,:)')')
