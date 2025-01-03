matdir = "E:\OneDrive - Harvard University\Mat_Statistics";
load(fullfile(matdir, "CD_BigGAN_Hessian_Manifold_Stats.mat"),"HEStats");
%%
for Expi = 1:length(HEStats)
    fprintf("Experiment %d\n", Expi);
    disp(HEStats(Expi).meta)
    HEStats(Expi).noise.resp_sgtr_col;
    % HEStats(Expi).class.resp_sgtr_col
end
