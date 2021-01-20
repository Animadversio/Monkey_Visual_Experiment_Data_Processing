% Small utility code to visualize change from Pasupathy patch to faces
% D = load("E:\OneDrive - Washington University in St. Louis\ref_img_fit\Pasupathy\pasu_fit_code.mat");
D = load("N:\Users\Binxu\pasu_fit_code.mat");
pasu_code = D.pasu_code;
D = load("N:\Users\Binxu\genes4binxu.mat");
D = D.genes4binxu;
evol_code = D.all(D.gen==max(D.gen),:);

G = FC6Generator();
%%
i=50;j=1;
interp_codes = LERP(pasu_code(i,:), evol_code(j,:), linspace(0,1,11));
imgs = G.visualize(interp_codes);
for i = 1:size(imgs,4), imwrite(imgs(:,:,:,i), fullfile("N:\Users\Binxu",compose("Lin%02d.png", i))); end
savefast(fullfile("N:\Users\Binxu\interp_codes_lin.mat"), 'interp_codes')

interp_codes = SLERP(pasu_code(i,:), evol_code(j,:), linspace(0,1,11));
imgs = G.visualize(interp_codes);
for i = 1:size(imgs,4), imwrite(imgs(:,:,:,i), fullfile("N:\Users\Binxu",compose("Sph%02d.png", i))); end
savefast(fullfile("N:\Users\Binxu\interp_codes_sph.mat"), 'interp_codes')
