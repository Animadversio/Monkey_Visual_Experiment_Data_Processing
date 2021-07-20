function h = violinplot_cell(value_col, label_arr, varargin)
% Wrapper or another api for violinplot. 
% Support input cell array of vector values and array of labels. 
%
% Example: 
% cc_col = {cat(1,SUMUcorrvec_col{:,1}),...
%           cat(1,areacorrvec_col{:,1}),...
%           cat(1,SUMUcorrvec_col{:,2}),...
%           cat(1,areacorrvec_col{:,2}),...
%           cat(1,SUMUcorrvec_col{:,3}),...
%           cat(1,areacorrvec_col{:,3})};
% label_arr = ["SU-MU V1", "All V1",...
%              "SU-MU V4", "All V4",...
%              "SU-MU IT", "All IT"];
% violinplot_cell(cc_col, label_arr,'showData',true,'GroupOrder',cellstr(label_arr))
% 
Yvec = [];
labvec = [];
for i=1:numel(value_col)
    Yvec = cat(1,Yvec, reshape(value_col{i},[],1));
    labvec = cat(1, labvec, repmat(label_arr(i),numel(value_col{i}),1));
end
violinplot(Yvec,labvec,varargin{:})
end