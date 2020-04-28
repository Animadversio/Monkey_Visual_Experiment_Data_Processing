pasu_path = "imgs/pasupathy-wg-f-4-ori/";
imlist = dir(fullfile(pasu_path,"*.jpg"));
imlist = struct2cell(imlist);
imgnms = cellstr(imlist(1,:));
%% Collect all Pasupathy images and their 
global samp_num
Bdr_i = 2;samp_num = 150;
Bd_col = {};
FFT_col = {};
for imgi = 1:length(imgnms)
    img = imread(fullfile(pasu_path, imgnms{imgi}));
    img_bw = img(:,:,1) < 185;
    B=bwboundaries(img_bw,8); % get the boundary
    Xcoor = B{2}(1:end, 2); % get the boundary pixel index
    Ycoor = B{2}(1:end, 1);
    Zcoor = Xcoor + j * Ycoor;
    N = length(Xcoor);
    i_qr = linspace(1,N,samp_num);
    z_samp = interp1(1:N, Zcoor, i_qr);
    Ucoef = fft(z_samp);
    Bd_col{imgi} = B{2};
    FFT_col{imgi} = Ucoef;
end
%% Visualize the Spectrum of smooth contour
figure(60);clf;loglog(1:74,abs(Ucoef(2:75)),'-');hold on;
loglog(1:74,abs(Ucoef(end:-1:end-73)));
xlabel("F component")
ylabel("Amplitude")
% figure;hold on,plot(log(1:74),log(abs(Ucoef(2:75))));plot(log(1:74),log(abs(Ucoef(end:-1:end-73))))
% [r,m,b] = regression(log(abs(Ucoef(2:75))),log(1:74));
%% Regress the amplitude of coefficient onto the 
lmdl = fitlm([log(1:74),log(1:74)],...
       log([abs(Ucoef(2:75)),abs(Ucoef(end:-1:end-73))]))
%%
scales = exp(7.8 + -1.6 * log(1:74));

scales_all = zeros(1,samp_num); % scale of each FFT component
scales_all(1) = exp(7.8 + 1.6);
scales_all(2:75) = exp(7.8 + -1.6 * log(1:74));
scales_all(end:-1:76) = exp(7.8 + -1.6 * log(1:75));

figure;hold on
plot(real(Ucoef(2:75))./scales)
plot(imag(Ucoef(2:75))./scales)
plot(real(Ucoef(end:-1:end-73))./scales)
plot(imag(Ucoef(end:-1:end-73))./scales)
%%
reallmdl = fitlm([log(1:74),log(1:74)],...
       log(abs([real(Ucoef(2:75)),real(Ucoef(end:-1:end-73))])))
imaglmdl = fitlm([log(1:74),log(1:74)],...
       log(abs([imag(Ucoef(2:75)),imag(Ucoef(end:-1:end-73))])))
%% Fit the curve with a few component
imgi = 181;
UAmp = abs(FFT_col{imgi});
UAng = angle(FFT_col{imgi});
N = samp_num;
comp_n = 12;
FFTidx = [1:comp_n + 1, N-comp_n+1:N];
%1/length(UAmp)
theta = [0/N:1/N:1];%[0:length(UAmp)-1]/length(UAmp);%[0:0.001:1]%
freq = [0:N-1]';
Xfit = sum(UAmp(FFTidx)' .* cos(2*pi*freq(FFTidx) * theta + UAng(FFTidx)'), 1)/N;
Yfit = sum(UAmp(FFTidx)' .* sin(2*pi*freq(FFTidx) * theta + UAng(FFTidx)'), 1)/N;

figure
plot(Bd_col{imgi}(:,2),Bd_col{imgi}(:,1));hold on 
plot(Xfit,Yfit,'red')
axis image equal
set(gca,"YDir","reverse")
%% Plot the contour
comp_n = 12;
FFTidx = [1:comp_n + 1, N-comp_n+1:N];
[Xfit,Yfit] = contourGen(0.5 * FFT_col{100}(FFTidx) + 0.5 * FFT_col{105}(FFTidx));
figure;hold on
%plot(Bd_col{imgi}(:,2),Bd_col{imgi}(:,1));hold on 
plot(Xfit,Yfit,'color',[1,0,1])
[Xfit,Yfit] = contourGen(FFT_col{100}(FFTidx));
plot(Xfit,Yfit,'red')
[Xfit,Yfit] = contourGen(FFT_col{105}(FFTidx));
plot(Xfit,Yfit,'blue')
axis image equal
set(gca,"YDir","reverse")
%% Fill the contour to get image
W=224; H=224;
figure(57);clf;hold on;set(57,"position",[150,300,W,H])
%plot(Bd_col{imgi}(:,2),Bd_col{imgi}(:,1));hold on 
plot(Xfit,Yfit,'color',[0,0,0])
axis image equal off
set(gca,"YDir","reverse",'position',[0 0 1 1],'units','normalized')
xlim([0,300]);ylim([0,300])
F=getframe(57);
img_F = frame2im(F);
%%
fill_img_F = imfill(img_F(:,:,1)<128,floor([mean(Xfit),mean(Yfit)]));
figure, imshow(fill_img_F)
%%
[filled_img] = fillShapeGen(FFT_col{105}(FFTidx));
figure, imshow(filled_img)
%% Try to evolve from it! 
net = alexnet;
addpath ..\..\CMAES_optimizer_matlab\
%%
batch_img = [];
for i = 1:40
[filled_img] = fillShapeGen(FFT_col{i}(FFTidx));
filled_rgb = repmat(filled_img,1,1,3);
batch_img = cat(4,batch_img,filled_rgb);
end
acts=activations(net, batch_img, 'fc6');
squeeze(mean(acts,[1]))
%%
my_layer = "conv5";iChan=10;xi=7;yi=7;
my_layer = "fc6";iChan=10;xi=1;yi=1;
n_gen=100;Visualize=true;scatclr="b";

comp_n = 10;
FFTidx = [2:comp_n + 1, N-comp_n+1:N];
init_gene = [real(FFT_col{i}(FFTidx))./scales_all(FFTidx), imag(FFT_col{i}(FFTidx))./scales_all(FFTidx)];
Optimizer = CMAES_simple(length(init_gene),init_gene,struct("init_sigma",0.25,"Aupdate_freq", 200));
h = figure(1);clf;
h.Position = [210         276        1201         645];
annotation(h,'textbox',...
    [0.485 0.467 0.154 0.50],'String',split(printOptionStr(options),','),...
    'FontSize',14,'FitBoxToText','on','EdgeColor','none')
genes = init_gene;
codes_all = [];
scores_all = [];
generations = [];
mean_activation = nan(1,n_gen) ;
for iGen = 1:n_gen
    % generate pictures
    coef_n = floor(size(genes,2)/2); % this is noe the same as comp_n!it's 2 * comp_n +1
    Ucoefs = (genes(:,1:coef_n) + 1j * genes(:,1+coef_n:end)).*scales_all(FFTidx);
    batch_img = [];
    for imgj = 1:size(Ucoefs,1)
        [filled_img] = fillShapeGen(Ucoefs(imgj,:));
        filled_rgb = repmat(filled_img,1,1,3);
        batch_img = cat(4,batch_img, filled_rgb);
    end
    % get activations
    acts = activations(net,255*batch_img,my_layer,'OutputAs','Channels');
    act_unit = squeeze( acts(xi,yi,iChan,:) ) ;
    %disp(act_unit')
    % Record info 
    scores_all = [scores_all; act_unit]; 
    codes_all = [codes_all; genes];
    generations = [generations; iGen * ones(length(act_unit), 1)];
    % pass that unit's activations into CMAES_simple
    % save the new codes as 'genes'
    [genes_new,tids] = Optimizer.doScoring(genes, act_unit, true, struct());
    if Visualize
    set(0,"CurrentFigure",h)
    % plot firing rate as it goes
    subplot(2,2,1)
    mean_activation(iGen) = mean(act_unit) ;
    scatter(iGen*ones(1,length(act_unit)),act_unit,16,...
        'MarkerFaceColor',scatclr,'MarkerEdgeColor',scatclr,...
        'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
    plot(iGen, mean(act_unit) ,'r.','markersize',20)
    xlim([0, n_gen])
    ylabel("scores")
    xlabel("generations")
    hold on
    subplot(2,2,3)
    code_norms = sqrt(sum(genes.^2, 2));
    scatter(iGen*ones(1,length(code_norms)),code_norms,16,...
        'MarkerFaceColor',scatclr,'MarkerEdgeColor',scatclr,...
        'MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.4)
    plot(iGen, mean(code_norms) ,'r.','markersize',20)
    xlim([0, n_gen])
    ylabel("code norm")
    xlabel("generations")
    if class(Optimizer) == "ZOHA_Sphere_lr"
    title(num2str(Optimizer.mulist(Optimizer.istep + 1)))
    end
    hold on
    [meanPic] = fillShapeGen(mean(Ucoefs,1));
    set(0,"CurrentFigure",h);subplot(2,2,2);cla
    imshow(meanPic);
    axis image off
    
    [mxscore, mxidx]= max(act_unit);
    [max_img] = fillShapeGen(Ucoefs(mxidx,:));
    set(0,"CurrentFigure",h);subplot(2,2,4);cla
    imshow(max_img);
    title(num2str(mxscore))
    axis image off
    drawnow
	end
% 	if SaveImg
%         image_name = sprintf('%s_%03d_%02d_%02d_%02d_%02d.jpg',my_layer,iChan,t_unit,nrows,ncols,iGen) ;
%         imwrite( meanPic ,  ...
%             fullfile(my_final_path, my_layer, sprintf('%02d',t_unit), image_name ) , 'jpg')
%     end
    genes = genes_new;
end
%%

%%
function [filled_img] = fillShapeGen(Ucoef)
[Xfit,Yfit] = contourGen(Ucoef);
W=227; H=227;
if ~ishandle(58), figure(58);end 
set(0,"CurrentFigure",58)
set(58,"position",[150,300,W,H],"Visible","off");clf;
plot(Xfit,Yfit,'color',[0,0,0])
axis image equal off
set(gca,"YDir","reverse",'position',[0 0 1 1],'units','normalized')
xlim([-150,150]);ylim([-150,150])
F=getframe(58);
img_F = frame2im(F);
seed = [clip(floor(mean(Xfit)./300.*227),1,W), ...
        clip(floor(mean(Yfit)./300.*227),1,H)];
filled_img = imfill(img_F(:,:,1)<128, seed);
end

function [Xfit,Yfit] = contourGen(Ucoef)
global samp_num
N = samp_num;
UAmp = abs(Ucoef);
UAng = angle(Ucoef);
comp_n = floor(length(Ucoef)/2);
% FFTidx = [1:comp_n+1, N-comp_n+1:N];
FFTidx = [2:comp_n+1, N-comp_n+1:N];
theta = [0/N:1/N:1];%[0:length(UAmp)-1]/length(UAmp);%[0:0.001:1]%
freq = [0:N-1]';
Xfit = sum(UAmp' .* cos(2*pi*freq(FFTidx) * theta + UAng'), 1)/N;
Yfit = sum(UAmp' .* sin(2*pi*freq(FFTidx) * theta + UAng'), 1)/N;
end