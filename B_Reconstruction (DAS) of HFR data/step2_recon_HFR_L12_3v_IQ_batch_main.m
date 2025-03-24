% Script to process a batch of HFR data and save the IQ data.
% Uses following files:
%   step2_recon_HFR_L12_3v_IQ_batch_reconpart.m
%       reconstruct_DAS_get_vars.m
%   saveIQData_DAS.m
%   saveIQVideo.m

% Inputs:
%   A directory where (multiple) HFR _info.mat and _RcvData.bin files are
%   stored.
%   A directory where the results can be stored (files in this directory
%   with the same name will be overwritten).

% Outputs:
%   _IQData.bin containing the data
%   _IQData.info containing some info about the file
%   _IQData.png/_IQData.fig containing a hilbert transformed, log
%   compressed image of the first frame in the IQ data.
%   _IQData_logcomp_vido.avi A video of the data

% 27-01-2023 Jelle Plomp
clear all; close all; clc
%% Change current directory to git repo top dir add paths
current_dir = pwd;
workfolder = 'C:\Users\PlompJ\OneDrive - University of Twente\Matlab_git_br_pub\';
cd(workfolder)

allpaths = genpath(workfolder);
allpaths = strsplit(allpaths, ';');

for i=1:length(allpaths)
    if contains(allpaths{i},'B_Reconstruction (DAS) of HFR data') || contains(allpaths{i},'delay-and-sum-reconstruction')
        addpath(allpaths{i}) 
    end
end
%% Data loading
step2_recon_HFR_L12_3v_IQ_batch_dataloading_A_general

%% Some options
% SVD filtering
svd_receive = input("Apply SVD on received data? (true/false)");
if svd_receive
    % SVD tests 9-2-2023
%     min_reg_val_array = [1:2:20];
%     max_reg_val_array = [60:20:200];
    
    min_reg_val_array = 20;
    max_reg_val_array = 1000; % If NaN, we define it after loading first dataset
%     Nfrs_recon = 50; % To save some time, only 50 frames are reconstructed. SVD is performed over all frames.

else % These are not actually applied but it allows running the for-loop below even without applying svd.
    min_reg_val_array = 0;
    max_reg_val_array = 0;
end

use_raddec = true; % Use decay (compensating for radial scattering 1/r)?
max_angle = 45; 
att = 0;
near_field_sat = 400;
f_c=0.5; %MHz, cut-off freq !! kan dit misschien hoger nu f=15MHz?

dynRange = 55; %will overwrite displayInfo.dynrange
overwrite_resolution = true; % If this is set to true, then the resolution saved in ReconInfo_DAS will be replaced
overwrite_zlimits = true; % limits of reconstruction are changed

recon_choices.demod_rf = true;
recon_choices.demod_rf_method = 'rf2iq'; % rf_hilb or rf2iq

overwrite_Nfrs_recon = true;
Nfrs_recon_set = NaN; % NaN=Reconstruct all frames in the file.
            
%% Perform DAS reconstruction and save "IQ" data
% Note that the IQ data from the DAS reconstruction is not complex data.
callfrombatchscript = true;
for min_reg_val = min_reg_val_array
    
    for max_reg_val = max_reg_val_array
        if svd_receive
           display(['Svd receive data with min_reg_val = ', ...
               num2str(min_reg_val), ' and max_reg_val = ', ...
               num2str(max_reg_val)])
        end
        for i_file=1:length(files_info)
            fname_info = files_info(i_file).name;
            % Check if IQ data already exists
            data_dir_results_IQ = [data_dir_results, '\IQData'];
            if ~isfolder(data_dir_results_IQ)
               mkdir(data_dir_results_IQ)
            end
            f = find(fname_info=='_');
            fname = fname_info(1:f(end)-1);
            if svd_receive
                sname_IQData_pre = [data_dir_results_IQ, '\', fname, '_RcvSVD_min_', num2str(min_reg_val), '_max_', num2str(max_reg_val)];
            else
                sname_IQData_pre = [data_dir_results_IQ, '\', fname];
            end
            if recon_choices.demod_rf
                sname_IQData_pre = [sname_IQData_pre, recon_choices.demod_rf_method];
            end
            sname_IQData = [sname_IQData_pre, '_IQData'];
            if isfile([sname_IQData '_info.mat'])
                warning([sname_IQData ' already exists, so skipped'])
                continue
            end
            display(['Reconstructing file ', num2str(i_file), ' of', num2str(length(files_info))])
            
            % Reconstruction of IQ data
            if overwrite_Nfrs_recon
                Nfrs_recon = Nfrs_recon_set;
            end
            step2_recon_HFR_L12_3v_IQ_batch_reconpart
            %%
            % Edit the dynamic rane
            if exist('dynRange', 'var')
                displayInfo.dynRange = dynRange;
            end
            
            % Save the IQ data
            saveIQData_DAS
            
            % Save a video of the data
            nr_frames_video = min(300,Nfrs_recon); % Don't save all 3000 frames to video. If NaN, then all frames will be saved to video.
            saveIQVideo
            
            % Clear variables 
            clearvars -except fname data_dir_results dirname Nfrs_recon...
                use_raddec max_angle att near_field_sat f_c files_info...
                callfrombatchscript dynRange overwrite_resolution svd_receive...
                min_reg_val max_reg_val min_reg_val_array max_reg_val_array...
                overwrite_zlimits recon_choices overwrite_Nfrs_recon...
                Nfrs_recon_set
        end
    end
end