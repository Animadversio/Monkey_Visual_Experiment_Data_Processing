%% Evol Ref Images Analysis
Animal = "Beto";
MatStats_path = "C:\Users\ponce\OneDrive - Washington University in St. Louis\Mat_Statistics";
load(fullfile(MatStats_path, compose("%s_Evol_stats.mat", Animal)), 'EStats')
%%
figure(1);set(1,'position',[ 469   430   771   548])
result_dir = "E:\Evolution_Exp";
for Expi = 3:45
savedir = fullfile(result_dir,compose("Manifold_Evol%02d_chan%02d", ...
    Expi, EStats(Expi).units.pref_chan));
mkdir(savedir)
fprintf("Processing Exp %d\n",Expi)
rasters = rasters_new{Expi};
imgnm = EStats(Expi).imageName;
row_gen = EStats(Expi).stim.gen_msk;
row_nat = ~row_gen;%contains(imgnm, "nat") & cellfun(@(c) ~isempty(regexp(c(1:2), "\d\d")), imgnm);
all_gabor_ref = all(contains(EStats(Expi).ref.imgnm,"gab"));

nat_imgidx = EStats(Expi).ref.idx_arr;
nat_psth_avg = cellfun(@(idx) mean(rasters(:, :, idx), 3), ...
        nat_imgidx, 'UniformOutput', false);
nat_psth_sem = cellfun(@(idx) std(rasters(:, :, idx), 1, 3)/sqrt(length(idx)), ...
        nat_imgidx, 'UniformOutput', false);
for ch_j = 1:size(rasters,1)
set(0, "CurrentFigure", 1);clf;hold on
for i = 1:length(EStats(Expi).ref.imgnm)
    %plot(nat_psth_avg{i})
    shadedErrorBar([],nat_psth_avg{i}(ch_j,:),nat_psth_sem{i}(ch_j,:),'lineprops',...
        {'markerfacecolor',EStats(Expi).color_seq(i,:)}, 'transparent',1);
end
%set(gca, 'position', [0.05,0.05,0.9,0.9])
xlabel("time (ms)")
if all_gabor_ref
    title(compose("%s Manifold Exp %d\n PSTH of unit %s to ref images (gabors)",Animal,Expi,EStats(Expi).units.unit_name_arr(ch_j)))
else
    title(compose("%s Manifold Exp %d\n PSTH of unit %s to ref images",Animal,Expi,EStats(Expi).units.unit_name_arr(ch_j)))
end
saveas(1, fullfile(savedir, compose("ref_psth_%s.png",EStats(Expi).units.unit_name_arr(ch_j))))
%break
end
%break
end
%%
