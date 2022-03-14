%% Static visualization 
function [h1,ax1,ax2] = vis_selectivity_static(Snat,sorted,chan,unit,order, K, figh)
% thin wrapper around the interactive function to produce arguments for it.
% Parameters:
%   sorted: Boolean, sort the responses in ascending order (default false)
%   chan: Int, channel number to plot
%   unit: Int, unit number in that channel. 
%   order: "top" (default) or "bottom". 
%   K: Top K images to show in 2nd panel.. 
if nargin < 2, sorted = false; end
if nargin < 3, do_prefchan =true; 
else, do_prefchan=false; end 
if nargin < 4, unit = 1;end
if nargin < 5, order = "top";end
if nargin < 6, K = 16; end
if nargin < 7, figh = figure(); end
if do_prefchan
    unit_idx = Snat.units.pref_chan_id;
else
    unit_idx = find((Snat.units.spikeID == chan) & ...
                    (Snat.units.unit_num_arr == unit));
end
rspvec = Snat.resp.meanMat(unit_idx, :);
semvec = Snat.resp.semMat(unit_idx, :);
rspcol = cellfun(@(M) M(unit_idx,:)', Snat.resp.trial_col,'uni',0);
bsl = Snat.resp.bslmean(unit_idx);
imgnm_uniq = Snat.stim.imgname_uniq;
imgfps = Snat.stim.imgfps;
imgN = numel(rspvec);
[rsp_sort, sortidx] = sort(rspvec);
if sorted % sort the response from low to high to see the distribution
    rspvec = rspvec(sortidx);
    semvec = semvec(sortidx);
    rspcol = rspcol(sortidx);
    imgnm_uniq = imgnm_uniq(sortidx);
    imgfps = imgfps(sortidx);
end
idxcol = arrayfun(@(i) i*ones(numel(rspcol{i}),1),[1:numel(rspcol)]','uni',0);
idx_all_vec = cat(1,idxcol{:});
rsp_all_vec = cat(1,rspcol{:});
h1 = figure(figh);set(h1,'pos',[188     89    1365       655]);clf
ax1 = subplot(121);set(ax1,'pos',[0.05    0.10    0.43    0.85]);hold on
barH = shadedErrorBar([],rspvec,semvec);
scatter(idx_all_vec,rsp_all_vec,25,'bo');
line([0,numel(rspvec)],[bsl, bsl],'linestyle','-.','linewidth',1.5)
xlim([0,numel(rspvec)+1])
xlabel("image id",'fontsize',14); 
ylabel("Activation: rsp - bsl",'fontsize',14);

if sorted
    if strcmp(order,"top")
    idx2show = imgN:-1:imgN-K+1;
    elseif strcmp(order,"bottom")
    idx2show = 1:K;
    end
else
    if strcmp(order,"top")
    idx2show = sortidx(imgN:-1:imgN-K+1);
    elseif strcmp(order,"bottom")
    idx2show = sortidx(1:K);
    end
end
axes(ax1);hold on
SCT = scatter(idx2show,rspvec(idx2show),36,'ro');
SCT.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('img',imgnm_uniq(idx2show));
SCT.DataTipTemplate.Interpreter = 'none';

ax2 = subplot(122);set(ax2,'pos',[0.54    0.075    0.43    0.85]);hold on
% axes(ax2);cla(ax2)
montage(imgfps(idx2show))
% show_imgs(ax2, idx2show);
titlestr = compose("%s chan %s\nImage center [%.1f %.1f] size %.1f deg",Snat.meta.ephysFN,...
    Snat.units.unit_name_arr(unit_idx), Snat.stim.impos(1), Snat.stim.impos(2), Snat.stim.imsize_deg);
sgtitle(h1,titlestr)
end

% function show_imgs(ax2, imgid)
% % disp(imgid)
% for imgi = reshape(imgid,1,[])
%     if any(ax2.UserData.showed_ids==imgi) 
%         idx = find(ax2.UserData.showed_ids==imgi);
%         ax2.UserData.showed_ids(idx) = [];
%         ax2.UserData.img_col(idx) = [];
%         continue
%     else
%         ax2.UserData.showed_ids(end+1) = imgi;
%         img = imread(ax2.UserData.imgfps(imgi));
%         ax2.UserData.img_col{end+1} = img;
%     end
% end
% axes(ax2);cla(ax2)
% montage(ax2.UserData.img_col)
% % imshow(img)
% % title(ax2.UserData.imgfps(imgi),'interp','none')
% end