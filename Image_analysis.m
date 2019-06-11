clear
%%
image_dir = 'D:\Monkey_Data\2019-06-Evolutions\beto-190610a\BACKUP_3';
code_dir = 'D:\Monkey_Data\2019-06-Evolutions\beto-190610a';
% imglist = ls([image_dir,'\block*.jpg']);
% codelist = ls([code_dir,'\block*.mat']);
imglist = dir(fullfile(image_dir, 'block*.jpg'));
img_fn = {imglist.name};
img_fullfn = fullfile(image_dir, {imglist.name});
codelist = dir(fullfile(code_dir, 'block*.mat'));
code_fn = {codelist.name};
code_fullfn = fullfile(code_dir, {codelist.name});
%%
score_arr = [];
id_arr = {};
code_arr = [];
for i =1:numel(code_fullfn)
   data = load(code_fullfn{i});
   code_arr = cat(1, code_arr, data.codes);
   id_arr = cat(2, id_arr, data.ids);
   if i~=1
       score_arr = cat(1, score_arr, data.scores);
   end
   if i == numel(code_fullfn)
       score_arr = cat(1, score_arr, nan([length(data.ids), 1]));
   end
end

%%
tic
code_dist = pdist(code_arr);
code_coord3 = mdscale(code_dist, 3);
toc
%%
figure()
scatter3(code_coord3(:,1),code_coord3(:,2),code_coord3(:,3),9,1:4030)
xlabel("MDS1")
xlabel("MDS2")
xlabel("MDS3")
