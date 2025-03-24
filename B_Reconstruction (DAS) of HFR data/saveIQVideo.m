% Script to create a video from IQ frames. 
% Can be used both in batch and separately
%% STEP 1: load data (copied from RunPIV_Hadi_Wietske_Jelle.m)
% we have binary files and an info (.mat) file that we need to load. 
% please adjust so your data is loaded correctly; acquisitions should be
% called "Frames"
if exist('callfrombatchscript', 'var') && callfrombatchscript
    Frames = IQData_matrix;
%     x= ReconVars.x;
%     z = ReconVars.z;
    if svd_receive
        title_ = strcat(model, " ", voltage, ", PC ", num2str(min_reg_val), " to ", num2str(max_reg_val), " (of ", num2str(size(Sigma,2)), ")");
    else
        title_ = strcat(model, " ", voltage);
    end
    sname_video = [sname_IQData, '_logcomp_video.avi'];
else
    startdir = fullfile('E:\Processed Jelle\');
    [fname,dirname] = uigetfile(startdir);
    filepath = fullfile (dirname, fname);
    binpath = fullfile(dirname, [fname(1:find(fname=='_',1,'last')-1), '.bin']);
    
    load(filepath); 
    fid = fopen(binpath, 'r');
    Frames = fread(fid,prod(SaveShape),['*' DataType]);
    Frames = reshape(Frames,SaveShape); 
    
    sname_video = [dirname, fname(1:end-9), 'logcomp_video_manual.avi'];
    nr_frames_video = 300;
    % Your Frames need to be arranged as x,z,angles,time
    Frames = reshape(Frames, size(Frames,1), size(Frames,2), 1, []);
    
    disp('Data Loaded');
    title_ = "HFR data";

    x=ReconInfo_DAS.x_start:ReconInfo_DAS.pas:ReconInfo_DAS.x_end;
    z=ReconInfo_DAS.z_start:ReconInfo_DAS.pas:ReconInfo_DAS.z_end;
    displayInfo.dynRange = 45; % Overwrite
end
 
%% Step 2: Create video of log compressed IQdata.

close all
maxvalue = max(max(abs(Frames(:))));
Frames_log_comp = 20*log10(abs(Frames(:,:,:))/maxvalue);

% The video will be saved in the same folder where the IQData was obtained.
v = VideoWriter(sname_video);
v.FrameRate = 10;

open(v);

f1=figure();ax=gca;
f1.Position = [100, 100, 800, 400];
% Only for the first frame

% Since imagesc doesn't work well with hold on, we have to set some
% properties
set(ax,'NextPlot','replacechildren');
axis equal
 ax.XLim = [x(1), x(end)].*1e3;
ax.YLim = [z(1), z(end)].*1e3;

caxis([-displayInfo.dynRange 0])
xlabel("x (mm)")
ylabel("z (mm)")

ax.YDir = 'reverse';
c = colorbar; 
colormap('gray')
ylabel(c,'Intensity [dB]');
set(gcf, 'Color', 'w')
set(gca, 'FontSize', 16)
if isnan(nr_frames_video)
    nr_frames_video = size(Frames_log_comp,3);
end
for frame_i=1:nr_frames_video

    imagesc(x.*1e3,z.*1e3, Frames_log_comp(:,:,frame_i)', [-displayInfo.dynRange 0])
    time_ms = P.PRP_HFR*frame_i*1e-3; % Frame time in milliseconds
    title(strcat(title_, ", Fr. ", num2str(frame_i), ", t= ", sprintf('%.1f',time_ms), " ms"),'Interpreter', 'none')
    drawnow
%     pause(0.1)
    videoframe = getframe(gcf);
    writeVideo(v,videoframe)
end

close(v)
