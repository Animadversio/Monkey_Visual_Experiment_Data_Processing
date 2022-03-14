function [h1,ax1,ax2]=vis_selectivity_interact_core(pics_unique,imgfps,rspvec,semvec,bsl,react,figh)
% Example: 
%    [h1,ax1,ax2]=vis_selectivity_interact_core(Snat.stim.imgname_uniq,Snat.stim.imgfps,...
%                     Snat.resp.meanMat_pref,Snat.resp.semMat_pref,...
%                     Snat.resp.bslmean(Snat.units.pref_chan_id));
% Parameter:
%    react: if false, the interactivity will be disabled. 
if nargin < 6, react=true; end
if nargin < 7, figh=figure(); end
h1 = figure(figh);set(h1,'pos',[188     89    1365       655]);clf
figure(h1)
dcm = datacursormode(h1);
ax2 = subplot(122);set(ax2,'pos',[0.54    0.075    0.43    0.85]);hold on
set(ax2, 'UserData', struct("rspvec",rspvec,"semvec",semvec,...
    "imgnm_uniq",string(pics_unique),"imgfps",string(imgfps),...
    "showed_ids",[], "img_col",[]))
ax1 = subplot(121);set(ax1,'pos',[0.05    0.10    0.43    0.85])
barH = shadedErrorBar([],rspvec,semvec);
line([0,numel(rspvec)],[bsl, bsl],'linestyle','-.','linewidth',1.5)
xlim([0,numel(rspvec)+1])
xlabel("image id"); ylabel("Activation: rsp - bsl")
if react
   set(barH.mainLine,'ButtonDownFcn',@(h,evt)selectpoint(h,evt,ax2),...
      'HitTest','on')
end
barH.patch.HitTest = 'off';
barH.mainLine.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('img',pics_unique);
barH.mainLine.DataTipTemplate.Interpreter = 'none';
dcm.UpdateFcn()
end