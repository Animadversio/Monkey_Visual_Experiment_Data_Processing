%% Visualize FC6 space.



G = FC6Generator();
%%
saveDir = "O:\TraineeDepartTalk";
vid = VideoWriter(fullfile(saveDir,'FC6Vis'),'MPEG-4');
vid.FrameRate = 12;%vid.CompressionRatio = 20;
open(vid)
targets = 4*randn(15,4096);
for i = 1:size(targets,1)-1
interp_vec = SLERP(targets(i,:),targets(i+1,:),linspace(0,1,32)) ;
imgs = G.visualize(interp_vec);
for j = 1:size(imgs,4)
   writeVideo(vid,imgs(:,:,:,j))
end
end
close(vid)
%%
which G.visualize_movie
