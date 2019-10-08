channel = 16;

%%
score_mat = zeros(11,11,5);
cnt_mat = zeros(11,11);
for i =  -5:5
    for j =  -5:5
        cur_fn = sprintf('norm_%d_PC2_%d_PC3_%d', norm, i*ang_step, j*ang_step);
        img_idx = find(contains(Trials.imageName, cur_fn));
        cnt_mat(i+6, j+6) = length(img_idx);
        psths = rasters(channel, :, img_idx);
        scores = squeeze(mean(psths(1, 50:150, :))- mean(psths(1, 1:50, :)) );
        score_mat(i+6, j+6, 1:length(img_idx)) = scores;
    end
end
figure(10)
imagesc(-90:18:90, -90:18:90, sum(score_mat,3)./cnt_mat)
ylabel("PC 2 degree")
xlabel("PC 3 degree")
title("Tuning map on PC2 3 subspace")
shading flat
axis image
colorbar
%%
colorseq = brewermap(11, 'BrBG');
figure(1);clf;hold on 
for i = 1 % -5:5
    for j = -5:5
        cur_fn = sprintf('norm_%d_PC2_%d_PC3_%d', norm, i*ang_step, j*ang_step);
        img_idx = find(contains(Trials.imageName, cur_fn));
        % cnt_mat(i+6, j+6) = length(img_idx);
        psths = rasters(channel, :, img_idx);
        plot(1:200, mean(squeeze(psths),2), 'Color', colorseq(j+6, :))
        pause
    end
end

% plot(1:200,squeeze(psths))
