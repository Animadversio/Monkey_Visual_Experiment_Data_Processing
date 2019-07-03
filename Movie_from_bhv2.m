BHV_fn = 'D:\Alfa_freeviewing\190628_Alfa_freeViewing_basic.bhv2'; % 'C:\Users\ponce\OneDrive\Desktop\OneDrive_Binxu\OneDrive\190628_Beto_generate_parallel.bhv2';
output_dir = 'D:\Alfa_freeviewing\190628_Alfa_freeViewing_basic_webcam';
concat_fn = 'D:\Alfa_freeviewing\190628_Alfa_freeViewing_basic_webcam.avi';
mkdir(output_dir)
[data,MLConfig,TrialRecord] = mlread(BHV_fn);
%% Write different trials to different files 
for trial_i = 126:length(data)
    [video,ts]=mlreadwebcam(trial_i,1,BHV_fn);
    writer = VideoWriter(fullfile(output_dir, sprintf('Trial%04d.avi', trial_i)));
    if length(ts) > 2
        writer.FrameRate = 1000 * (length(ts)-1)/(ts(end)-ts(1)); %  estimated frame rate 
    else
        writer.FrameRate = 15; % if too few frames, use this average value. 
    end
    open(writer)
    writeVideo(writer,video);
    close(writer);
end
%% Write a concatenated full video 
writer = VideoWriter(concat_fn);
writer.FrameRate = 15;
open(writer);
for trial_i = 1:length(data)
    [video,ts]=mlreadwebcam(trial_i,1,BHV_fn);
    writeVideo(writer,video);
end
close(writer);