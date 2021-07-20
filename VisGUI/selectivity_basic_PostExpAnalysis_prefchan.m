function selectivity_basic_PostExpAnalysis_prefchan(Sel_bhv2, imgdir)
if nargin == 0
Sel_bhv2 = string(ls([datestr(now,'yymmdd'),'*.bhv2']))
return
end
if nargin < 2
    imgdir = pwd; 
end
if strcmp(imgdir,'')
    imgdir = pwd; 
    fprintf("Use image directory %s\n",imgdir)
end

% Sel_bhv2 = string(ls([datestr(now,'yymmdd'),'*.bhv2']))
% imgdir = "C:\Users\Ponce lab\Documents\ml2a-monk\generate_BigGAN\2021-04-16-12-05-37\CCFactor_vgg16\img";
%%
tic
[B,MLConfig,TrialRecord] = mlread(Sel_bhv2);
fprintf('loaded file, took %d s\n', toc)
% rand_seed = round( rand(1)*100 );
%%
pics_all_withPath = vertcat(TrialRecord.User.stimulusID{:}) ;
[~,pics_all] = cellfun(@fileparts,pics_all_withPath,'Uni',false) ; 
pics_unique = unique(pics_all) ;
pics_n = length(pics_unique) ; 

xy_all = TrialRecord.User.xy ; 
xy_unique = unique(xy_all,'rows') ; 
xy_n = size(xy_unique,1);
assert(xy_n==1)

spike_times_cell = TrialRecord.User.allEventInfo ; 

all_chans = cellfun(@(x) x(:,2) , spike_times_cell ,'Uni',false) ; 
all_chans = vertcat(all_chans{:}) ; 
chans_unique = unique(all_chans) ; 
chans_unique(chans_unique == 257) = [] ; 
chans_n = length(chans_unique) ; 
TrialRecord.User.chans_n = chans_n ; 
TrialRecord.User.chans_unique = chans_unique ; 

TrialRecord.User.bkgwindow = [0 40];
TrialRecord.User.evkwindow = [51 200];
CHANNUM = 64;MAXUNUM=4;
responseTensor = nan(CHANNUM,MAXUNUM+1) ; % basically channel by image id
resp_trial_col = {}; % cell(CHANNUM,MAXUNUM+1);
for iXY = 1:xy_n
    t_xy = xy_unique(iXY,:) ; 
    S_xy_i = xy_all(:,1) == t_xy(1) & xy_all(:,2) == t_xy(2) ;
    
    for iPic = 1:pics_n
        t_pic = pics_unique{iPic} ; 
        S_pic_i = strcmp(pics_all, t_pic) ;
        S_final = S_xy_i & S_pic_i ; 
        
        ts = spike_times_cell(S_final) ; 

        tStimStartCode = TrialRecord.User.stimulusStartCode(S_final) ;
        
        [spike_rates_mat,spike_rate_tsr_trials] = selectivity_basic_PostExperimentAnalysis_returnMat(ts ,tStimStartCode,TrialRecord) ; 
        responseTensor(:,:,iPic,iXY) = spike_rates_mat ; 
        resp_trial_col{iPic,iXY} = spike_rate_tsr_trials;
        responseTensor_sem(:,:,iPic,iXY) = std(spike_rate_tsr_trials,1,3,'omitnan') ./ sqrt(sum(~isnan(spike_rate_tsr_trials),3));
        
    end % of iPic
    
end % of iXY
%% Text reprsentation 
iCh = TrialRecord.User.prefChan;
iUnit = 1;
rspvec_pfch = squeeze(responseTensor(iCh,iUnit+1,:));
semvec_pfch = squeeze(responseTensor_sem(iCh,iUnit+1,:));
for iPic = 1:pics_n
    imgnm = pics_unique{iPic} ; 
    fprintf("%s %.2f (%.2f)\n", imgnm, rspvec_pfch(iPic), semvec_pfch(iPic));
end
%%
[imgfps, mapper] = map2fullpath(pics_unique, imgdir);
%%
parts = split(Sel_bhv2,".");
imgnms = string(ls(fullfile(imgdir,"*.png")));
matnm = fullfile(compose("%s_prefchan_rsp.mat",parts(1)));
save(matnm,'rspvec_pfch','semvec_pfch',...
    'responseTensor','responseTensor_sem','resp_trial_col','imgfps')
%% Interactive Analysis of selectivity experiment
h1 = figure;set(h1,'pos',[519        -950        1365         440])
% h2 = figure;
global ax2 ax1
% imgfps = strcat(pics_unique,'.png');
figure(h1)
dcm = datacursormode(h1);
ax1 = subplot(121);set(ax1,'pos',[0.05    0.10    0.43    0.85])
barH = shadedErrorBar([],rspvec_pfch,semvec_pfch);
xlabel("image id"); ylabel("Activation: rsp - bsl")
set(barH.mainLine,'ButtonDownFcn',@selectpoint,...(~,~)disp('line')
   'HitTest','on')
barH.patch.HitTest = 'off';
barH.mainLine.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('img',pics_unique);
barH.mainLine.DataTipTemplate.Interpreter = 'none';
ax2 = subplot(122);set(ax2,'pos',[0.54    0.075    0.43    0.85]);hold on
set(ax2, 'UserData', struct("imgfps",string(imgfps),"showed_ids",[], "img_col",[]))
dcm.UpdateFcn()
end
%%
% show_imgs(ax2, 2)
function [coordSelected, minIdx] = selectpoint(hObj, event)
global ax2
% disp(event)
% disp(hObj)
X = hObj.XData; 
Y = hObj.YData; 
coordinates = [reshape(X,[],1),reshape(Y,[],1)];
pt = event.IntersectionPoint; 
x = pt(:,1);y = pt(:,2);
dist = pdist2([x,y],coordinates);            %distance between your selection and all points
[~, minIdx] = min(dist);            % index of minimum distance to points
coordSelected = coordinates(minIdx,:);
disp(coordSelected)
dt = datatip(hObj,x,y,'DeleteFcn',@(dtob,~)show_imgs(ax2, int32(coordSelected(1))));
show_imgs(ax2, int32(coordSelected(1)))
end
function show_imgs(ax2, imgid)
disp(imgid)
for imgi = reshape(imgid,1,[])
    if any(ax2.UserData.showed_ids==imgi) 
        idx = find(ax2.UserData.showed_ids==imgi);
        ax2.UserData.showed_ids(idx) = [];
        ax2.UserData.img_col(idx) = [];
        continue
    else
        ax2.UserData.showed_ids(end+1) = imgi;
        img = imread(ax2.UserData.imgfps(imgi));
        ax2.UserData.img_col{end+1} = img;
    end
end
axes(ax2);cla(ax2)
montage(ax2.UserData.img_col)
% imshow(img)
% title(ax2.UserData.imgfps(imgi),'interp','none')
end

function [imgfps, mapper] = map2fullpath(picnm_arr, imgdir)
imgnm_wsfx = deblank(string([ls([imgdir+"\*.png"]);ls([imgdir+"\*.jpg"]);ls([imgdir+"\*.jpeg"])]));
[~,imgnms,sfxs] = fileparts(imgnm_wsfx);
mapper = containers.Map(imgnms,fullfile(imgdir, imgnm_wsfx));
imgfps = cellfun(@(I)mapper(I),picnm_arr,'uni',0);
end