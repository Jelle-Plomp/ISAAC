clear all; clc
% Parameters to be set before running script
gitrepo_topdir = 'C:\Users\PlompJ\OneDrive - University of Twente\Matlab_git_br_pub\'; % replace
startdir = fullfile([gitrepo_topdir '\data_examples\IQData']); %Where to start data selection

do_image_enhancement=0; % Apply image enhancement?

% Apply Hilbert transform and/or log compression on the reconstructed
% frames?
use_hilbert = false;
use_log_comp = false;
% Load a mask (if false, user will be asked to draw it)
loadMask = true;

%% Change directory
workdir = [gitrepo_topdir 'C_Load_reconstructed_HFRData\'];
cd(workdir)
addpath([workdir 'ImageEnhancement_Cheung'])
addpath([workdir 'Masking'])

%% Selecting data
[fname_IQ,dirname_IQ] = uigetfile({'*_IQData_info.mat'},"Select IQ data info file", startdir );

%% Load IQ data
filepath = fullfile (dirname_IQ, fname_IQ);
binpath = fullfile(dirname_IQ, [fname_IQ(1:find(fname_IQ=='_',1,'last')-1), '.bin']);

load(filepath); 

fid = fopen(binpath, 'r');
Frames = fread(fid,prod(SaveShape),['*' DataType]);
Frames = reshape(Frames,SaveShape); 
fclose(fid);


% Rearrange as as x,z,angles,time
Frames = permute(Frames, [2 1 3]);
Frames = reshape(Frames, size(Frames,1), size(Frames,2), 1, []);


if use_hilbert
    Frames = abs(hilbert(Frames));
end
use_log_comp = false;
if use_log_comp
    Frames = 20.*log10(Frames/max(Frames,[],'all'));
end
disp('Data Loaded');
%% Load or draw a mask
maskdir = extractBefore(dirname_IQ, 'IQData\');
if loadMask
   [fname_mask,dirname_mask] = uigetfile(maskdir);
    load(fullfile(dirname_mask,fname_mask))
else
    % Draw a mask on reconstructed data (on which no clutter filter has been applied)
    mean_frame_Data = load([dirname_IQ, extractBefore(fname_IQ, 'MHz') 'MHz_averageframe_with_info.mat']);
    mean_frame = mean_frame_Data.mean_frame;
    mean_frame_log=20*log10(mean_frame'./max(mean_frame,[],"all"));
    % Draw mask
    [imMask, ROI, mPoly] = drawMask(mean_frame_log);
    % Save the mask
    maskpath = [maskdir, extractBefore(fname_IQ, '_IQData_info'), '_mask.mat'];
    save(maskpath,'imMask','ROI','mPoly');
end

%% Image enhancement (see separate script main_do_image_enhancement)

if do_image_enhancement
    frame_plot = 3;
    test_plot_regularisers = true;
    std_mm = 1; % Std with with the lumen is smoothed in the attenuation correction by cheung.
    main_do_image_enhancement
end
%% Some information from the data that could be relevant to perform PIV
dx = ReconInfo_DAS.pas.*1e3;  %[mm]
dz = ReconInfo_DAS.pas.*1e3;  %[mm]
dt = P.PRP_HFR*1e-6; %[s] pulse repetition time

% Centre frequency of transducer
f=P.TransFrequency_true

% Frequency used to calculate wavelengths in Verasonics variables and used
% to determine the sampling rate (this is 2x P.TransFrequency_true) 
f_fake = Trans.Frequency;

% Origin of coordinates
origin = [ReconInfo_DAS.z_start*1e3-0.5*dz, ReconInfo_DAS.x_start*1e3-0.5*dx]; % Position of top left corner in mm

% Speed of sound used for reconstruction   
c = ReconInfo_DAS.c;
