global samp_num
samp_num = 150;
%% real, imag paramterization
my_layer = "conv5";iChan=10;xi=7;yi=7;n_gen=100;
Visualize=true;scatclr="b";
comp_n = 10;
FFTidx = [2:comp_n + 1, N-comp_n+1:N];
init_gene = [real(FFT_col{183}(FFTidx))./scales_all(FFTidx), imag(FFT_col{183}(FFTidx))./scales_all(FFTidx)];
Optimizer = CMAES_simple(length(init_gene),init_gene,struct("init_sigma",0.25,"Aupdate_freq", 200));
h = figure(5);clf;
h.Position = [210         276        1201         645];
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
%% Complex number parametrization... not good
my_layer = "conv5";iChan=10;xi=7;yi=7;n_gen=100;
Visualize=true;scatclr="b";
comp_n = 10;
FFTidx = [2:comp_n + 1, N-comp_n+1:N];
init_gene = FFT_col{183}(FFTidx)./scales_all(FFTidx);
Optimizer = CMAES_simple(length(init_gene),init_gene,struct("init_sigma",0.25,"Aupdate_freq", 200));
h = figure(5);clf;
h.Position = [210         276        1201         645];
genes = init_gene;
codes_all = [];
scores_all = [];
generations = [];
mean_activation = nan(1,n_gen) ;
for iGen = 1:n_gen
    % generate pictures
    coef_n = floor(size(genes,2)); % this is noe the same as comp_n!it's 2 * comp_n +1
    Ucoefs = (genes(:,1:coef_n)).*scales_all(FFTidx);
    batch_img = [];
    for imgj = 1:size(Ucoefs,1)
        [filled_img] = fillShapeGen(Ucoefs(imgj,:));
        %filled_rgb = repmat(filled_img,1,1,1);
        batch_img = cat(4,batch_img, filled_img);
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
    code_norms = sqrt(sum(abs(genes).^2, 2));
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
%% Amp and Angle parametrization, Seems much better
% my_layer = "conv5";iChan=10;xi=7;yi=7;
my_layer = "fc6";iChan=10;xi=1;yi=1;
n_gen=100;Visualize=true;scatclr="b";
imgi = 58;
comp_n = 5;FFTidx = [2:comp_n + 1, N-comp_n+1:N];
init_gene = [abs(FFT_col{imgi}(FFTidx))./scales_all(FFTidx), angle(FFT_col{imgi}(FFTidx))];
Optimizer = CMAES_simple(length(init_gene),init_gene,struct("init_sigma",.5,"Aupdate_freq", 200));
h = figure(7);clf;
h.Position = [210         276        1201         645];
annotation(h,'textbox',...
    [0.485 0.467 0.154 0.50],'String',split(printOptionStr(Optimizer.opts),','),...
    'FontSize',14,'FitBoxToText','on','EdgeColor','none')
genes = init_gene;
codes_all = [];
scores_all = [];
generations = [];
mean_activation = nan(1,n_gen) ;
for iGen = 1:n_gen
    % generate pictures
    coef_n = floor(size(genes,2)/2); % this is noe the same as comp_n!it's 2 * comp_n +1
    Ucoefs = (genes(:,1:coef_n)).*scales_all(FFTidx) .* exp(1j .* genes(:,coef_n+1:end));
    batch_img = [];
    for imgj = 1:size(Ucoefs,1)
        [filled_img] = fillShapeGen(Ucoefs(imgj,:));
        batch_img = cat(4,batch_img, filled_img);
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
    title(compose("AlexNet-%s-Ch%d-%d,%d",my_layer,iChan,xi,yi))
    hold on
    subplot(2,2,3)
    code_norms = sqrt(sum(abs(genes).^2, 2));
    scatter(iGen*ones(1,length(code_norms)),code_norms,16,...
        'MarkerFaceColor',scatclr,'MarkerEdgeColor',scatclr,...
        'MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.4)
    plot(iGen, mean(code_norms) ,'r.','markersize',20)
    xlim([0, n_gen])
    ylabel("code norm")
    xlabel("generations")
    title(compose("Angle, amplitude space %d components",comp_n))
    hold on
    [meanPic] = fillShapeGen(mean(Ucoefs,1));
    set(0,"CurrentFigure",h);subplot(2,2,2);cla
    imshow(meanPic);
    title(num2str(mean(act_unit)))
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
[filled_img] = fillShapeGen(Ucoefs(2,:));
figure();imshow(filled_img);
%%
function [filled_rgb] = fillShapeGen(Ucoef)
[Xfit,Yfit] = contourGen(Ucoef);
% pgon = polyshape(Xfit,Yfit);
W=227; H=227; I = ones(H,W,3);
filled_img = poly2mask(W/2+Xfit/300*W,H/2+Yfit/300*H,W,H);
filled_rgb = single(repmat(filled_img,1,1,3));
% insertShape(I,'FilledPolygon')

% if ~ishandle(58), figure(58);set(58,"position",[150,300,W,H],"Visible","off");end 
% set(0,"CurrentFigure",58);clf;
% plot(pgon,'FaceColor','k','EdgeColor','none','FaceAlpha',1)%[0,0,0]
% set(gca,"YDir","reverse",'position',[0 0 1 1],'units','normalized')
% axis image equal off
% xlim([-150,150]);ylim([-150,150])
% F=getframe(58);
% filled_img = frame2im(F);

% img_F = frame2im(F);
% seed = [clip(floor(mean(Xfit)./300.*227),1,W), ...
%         clip(floor(mean(Yfit)./300.*227),1,H)];
% filled_img = imfill(img_F(:,:,1)<128, seed);
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