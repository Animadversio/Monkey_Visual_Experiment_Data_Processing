mat_dir = "E:\OneDrive - Washington University in St. Louis\Mat_Statistics"; 
for Animal = ["Both"] %, "Beto""Beto"
load(fullfile(mat_dir, Animal + "_BigGAN_FC6_Evol_Stats.mat"), 'BFEStats'); 
end
%%
[codes_all, img_ids, code_geni] = load_codes_all(BFEStats(3).meta.stimuli, 2);
%%
savedir = "E:\Network_Data_Sync\BigGAN_latent_codes";
%%
for Expi = 1:numel(BFEStats)
for threadi = [1, 2]
    try
        [codes_all, img_ids, code_geni] = load_codes_all(BFEStats(Expi).meta.stimuli, threadi);
        img_ids = cellstr(img_ids);
        save(fullfile(savedir, compose("Exp%03d_thread%d_all_codes.mat",Expi,threadi)),...
            "codes_all", "img_ids", "code_geni")
    catch exception
    % Handle the error
    disp([Expi, threadi])
    disp(['An error occurred: ', exception.message]);
    end
end
end