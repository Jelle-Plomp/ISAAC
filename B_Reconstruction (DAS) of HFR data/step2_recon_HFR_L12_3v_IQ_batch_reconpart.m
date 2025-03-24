% Script for DAS reconstruction of data of the L12-3v probe

% Script by Jelle Plomp (11-1-2023), algorithm based on code by Guillaume
% Lajoinie and others (Charlotte Nawijn, Macy Vreman, ...).

% This particular version is adapted for batch processing.

% Description of some of the inputs:
%   acquisition_correct: if true, it assumes acquisition was done with the
%       step_1_MAIN_acquisition_HFR_L12_3v script, which results in a
%       ReconInfo_DAS struct. If false, this struct should be defined
%       manually under "For measurements not performed using the new script,
%       define ReconInfo_DAS"
%   Nfrs_recon: number of frames to reconstruct from the RF data.
%   recon_choices
%       .demod_rf: apply demodulation on rf data
%       .demod_rf_method: use 'hilbert' or 'rf2iq' (rf2iq worked best)


% The script assumes only one aperture and one angle is used.
% (The input superframes are re-ordered such that all frames are just in
% subsequent order).

% The reconstruction algorithm was edited such that everything that remains
% constant over all frames is only calculated once. The calculations are
% done in the function file reconstruct_DAS_get_vars.m.


%% Some parameters (which we may also want to change from the batch script)
if ~exist("Nfrs_recon") 
    Nfrs_recon = 20; % Number of frames to reconstruct.
end
if ~exist("acquisition_correct")
    acquisition_correct = true; % false if data to be reconstructed was not acquired using the step_1 script adapted to this reconstruction.
end
%% Loading Data
f = find(fname_info=='_');
fname = fname_info(1:f(end)-1); % last '_'
fdata = [fname,'_Rcvdata.bin'];
finfo = [fname,'_info.mat'];
if ~exist('load_infofile_externally', 'var') || ~load_infofile_externally
    load(fullfile(dirname, finfo),'Resource','Trans','TX','TW','P','TGC',...
        'TPC','ReconInfo_DAS','RcvDataType','RcvDataShape','Receive','displayInfo');
end
tic
disp('Loading RcvData');


fid = fopen(fullfile(dirname,fdata));
RcvData0 = fread(fid,prod(RcvDataShape),['*' RcvDataType]);
RcvData0 = reshape(RcvData0,RcvDataShape);
fclose(fid);
disp('RcvData Loaded');
toc

%% In case you want to change the speed of sound, change below.
Resource.Parameters.speedOfSound = Resource.Parameters.speedOfSound;

ReconInfo_DAS.c = Resource.Parameters.speedOfSound;

if exist('overwrite_resolution', 'var') && overwrite_resolution
    % Calculate the resolution to be half the wavelength (based on
    % Trans.frequency, which may be the 'fake' frequency)
    ReconInfo_DAS.pas = 0.5 * ReconInfo_DAS.c/(Trans.frequency*1e6);
    
end
if exist('overwrite_zlimits', 'var') && overwrite_zlimits
    % Calculate the resolution to be half the wavelength (based on
    % Trans.frequency, which may be the 'fake' frequency)
    ReconInfo_DAS.z_start = 24e-3;
    ReconInfo_DAS.z_end = 34e-3;
    
end
%% For measurements not performed using the new script, define ReconInfo_DAS
if acquisition_correct == false
    P.HFRBufFrames = P.numFrames;
    P.numAcqsSuperFrame = P.numAcqs;
    P.VSXEndReceiveSample = Resource.RcvBuffer(2).rowsPerFrame/P.numAcqsSuperFrame;

    ReconInfo_DAS.pas = 0.1e-3;
    ReconInfo_DAS.x_start = -15e-3;
    ReconInfo_DAS.x_end =  15e-3;
    ReconInfo_DAS.z_start = 0;
    ReconInfo_DAS.z_end = P.endDepth_mm*1e-3;
end
%% Reshape the RcvData
Npt = Receive.endSample-Receive.startSample+1;
Nfrstot = P.HFRBufFrames*P.numAcqsSuperFrame;

if isnan(Nfrs_recon) 
    Nfrs_recon = Nfrstot; % Reconstruct all frames.
end
if acquisition_correct == false
    % Re-order the data and get rid of empty elements.
    order = Trans.Connector(Receive.aperture:Receive.aperture+Resource.Parameters.numRcvChannels-1,1);
    RcvData = double(RcvData0(:,order,:)); % Assuming we have frames stacked in 3rd dim
else
    % Here we assume the data is already sorted in the HFR saving.
    RcvData = double(RcvData0); % Saving is done in int16 (original buffer data type), for processing, double is better
end

% Reshape (get rid of superframes, just stack each frame on top of each
% other).
RcvData_sorted = zeros(P.VSXEndReceiveSample, size(RcvData,2), Nfrstot);
for superframe_j = 1:P.HFRBufFrames
    for frame_i=1:P.numAcqsSuperFrame
        start_frame_i = (frame_i-1)*P.VSXEndReceiveSample+1;
        end_frame_i = start_frame_i+P.VSXEndReceiveSample-1;
        frame_num_new = (superframe_j-1)*P.numAcqsSuperFrame+frame_i;
        RcvData_sorted(:, :,frame_num_new) = RcvData(start_frame_i:end_frame_i, :,superframe_j);
    end
end
clearvars RcvData RcvData0
%% Get variables needed for reconstruction (constant for all frames)
ReconVars = reconstruct_DAS_get_vars_sparsematrix(Trans, ReconInfo_DAS, Receive, Resource, TW, size(RcvData_sorted(:, :,1)), max_angle, use_raddec, att,(recon_choices.demod_rf && strcmp(recon_choices.demod_rf_method, 'rf2iq')));
clear ProcTimes
% The variables are 'unpacked' from the struct here since that makes the 
% reconstruction faster.
att_vec = ReconVars.att_vec';
f_axis = ReconVars.f_axis;
lx = ReconVars.lx;
lz =  ReconVars.lz;
DAS_matrix_sparse = ReconVars.DAS_matrix_sparse;
clearvars ReconVars
% Used for plotting:
x=ReconInfo_DAS.x_start:ReconInfo_DAS.pas:ReconInfo_DAS.x_end;
z=ReconInfo_DAS.z_start:ReconInfo_DAS.pas:ReconInfo_DAS.z_end;

%% Save the average of frames 1 to 100 (without SVD)
% This mean image can for example be used to draw a mask for the PIV analysis.
frames_orig = zeros(lx,lz,min(100, Nfrstot));
for frame_i = 1:min(100, Nfrstot)% if total number of frames smaller than 100, use total
    RF_temp = preproc_RcvFrame(RcvData_sorted(:,:,frame_i), att_vec,near_field_sat,f_axis,f_c);
    frames_orig(:,:,frame_i)=do_reconstruction(RF_temp, DAS_matrix_sparse, recon_choices, Receive.decimSampleRate*1e6, TW.Parameters(1)*1e6, lx,lz);
end
frames_orig_mean = mean(frames_orig,3);
savename_frame = [extractBefore(sname_IQData, 'MHz') 'MHz_averageframe_with_info'];
save_mean_image(frames_orig_mean, x, z, displayInfo, finfo, files_info(i_file), savename_frame,P,ReconInfo_DAS,Receive,Trans,TX)
%% SVD filtering of Receive data (Added 8-2-2023)
SVD_ensemblesize = 3000; % If the data has a lot of frames, might want
% to do SVD in multiple ensembles since this is computationally very
% expensive. Was not used for the paper, always used 1 large ensemble.
if exist('svd_receive', 'var') && svd_receive
    % Reconstruct fifth frame without SVD
    RF_temp = preproc_RcvFrame(RcvData_sorted(:,:,5), att_vec,near_field_sat,f_axis,f_c);
    img_orig=do_reconstruction(RF_temp, DAS_matrix_sparse, recon_choices, Receive.decimSampleRate*1e6, TW.Parameters(1)*1e6, lx,lz);
    if size(RcvData_sorted,3)>SVD_ensemblesize && rem(size(RcvData_sorted,3),SVD_ensemblesize)~=0
        error("Number of frames should be a multiple of the SVD ensemble size")
    elseif size(RcvData_sorted,3)>SVD_ensemblesize && rem(size(RcvData_sorted,3),SVD_ensemblesize)==0
        % Do the svd in ensembles
        N_ensem = size(RcvData_sorted,3)/SVD_ensemblesize;
        warning([num2str(N_ensem) ' ensembles of (' num2str(SVD_ensemblesize) ' frames are used for SVD'])
    elseif  size(RcvData_sorted,3)==SVD_ensemblesize % Do SVD over all frames
        N_ensem = 1;
        warning(['1 ensemble of (' num2str(SVD_ensemblesize) ' frames is used for SVD'])
    elseif size(RcvData_sorted,3)<SVD_ensemblesize % Do SVD over all frames
        N_ensem = 1;
        SVD_ensemblesize = size(RcvData_sorted,3);
        warning(['Instead of original SVD ensemble size (' num2str(SVD_ensemblesize) '), the total number of frames (' num2str(size(RcvData_sorted,3)) ') is used as ensemble size'])
    end
    display("Converting to single")
    bf_temp_all = single(RcvData_sorted(:,:,:));
    display("Converting to single is done")
    
    for ensem_z = 1:N_ensem
        % SVD part
        bf_temp = bf_temp_all(:,:,(ensem_z-1)*SVD_ensemblesize+1:ensem_z*SVD_ensemblesize);
        bf_temp=squeeze(bf_temp(:,:,:));
        size_bf_temp = size(bf_temp);
        % Casorati matrix: size = nr. values in frame x nr. frames
        S = reshape(bf_temp,[size(bf_temp,1)*size(bf_temp,2) size(bf_temp,3)]);    
        S_mean = mean(S, 2);
        S = S - S_mean; % Centering (i.e. remove mean of pixel across frames)    
        % perform SVD
        tic_svd = tic;
        [U,Sigma,V] = svdecon(S);
        U_cell{1,ensem_z} = U;
        Sigma_original_cell{1,ensem_z} = Sigma;
        V_cell{1,ensem_z}=V;
    
    % %     [U,Sigma,V] = svd(S, "econ");
        
        toc_svd = toc(tic_svd);
        disp(['SVD ensemble (' num2str(ensem_z) '/' num2str(N_ensem) ') took' num2str(toc_svd) ' seconds.'])
    end 
    clearvars S
        
    if isnan(min_reg_val) | isnan(max_reg_val) % If the SVD thresholds have not been defined yet, give user prompt.
        start_with_prompt = 1;
        
    else
        start_with_prompt = 0;
   
    end
        % Plot component strength
        
        figure(1);clf(1);
        subplot(1,4,[1:2]); hold on
        for ensem_z = 1:N_ensem
            Sigma_temp = Sigma_original_cell{1,ensem_z};
            plot(diag(Sigma_temp*Sigma_temp')/(size(Sigma_temp,1)-1))
        end
        grid on; grid minor
        ylabel(strcat("Feature strength \sigma\sigma'/", num2str(size(Sigma_temp,1)-1)))
        xlabel("Feature nr.")
        title("Feature strength")
        subplot(1,4,3);
        hold on
        for ensem_z = 1:N_ensem
            Sigma_temp = Sigma_original_cell{1,ensem_z};
            plot(diag(Sigma_temp*Sigma_temp')/(size(Sigma_temp,1)-1))
        end
        grid on; grid minor
        ylabel(strcat("Feature strength \sigma\sigma'/", num2str(size(Sigma_temp,1)-1)))
        xlabel("Feature nr.")
        title("Feature strengts PC 1 to 100")
        xlim([0 100])
        subplot(1,4,4)
        hold on
        for ensem_z = 1:N_ensem
            Sigma_temp = Sigma_original_cell{1,ensem_z};
            plot(diag(Sigma_temp*Sigma_temp')/(size(Sigma_temp,1)-1))
        end
        grid on; grid minor
        ylabel(strcat("Feature strength \sigma\sigma'/", num2str(size(Sigma_temp,1)-1)))
        xlabel("Feature nr.")
        title("Feature strengths PC end")
        xlim([size(Sigma_temp,1)-200 size(Sigma_temp,1)])
        satisfied_with_SVD = 0;
        while satisfied_with_SVD ==0
            if start_with_prompt
                pause(2)
                prompt = {'min_reg_val (min.=1)','max_reg_val (max. = nr. frames)'};
                dlgtitle = 'Input SVD limits';
                fieldsize = [1 45; 1 45];
                definput = {'1','3000'};
                answer = inputdlg(prompt,dlgtitle,fieldsize,definput);
                min_reg_val = max([1 str2double(answer{1})]);
                max_reg_val = min([size(bf_temp,3) str2double(answer{2})]);
            end
            % =================== Test frame based on 1st ensemble ========
            display(['Svd receive data with min_reg_val = ', ...
               num2str(min_reg_val), ' and max_reg_val = ', ...
               num2str(max_reg_val)])

            % remove modes
            tic_1=tic;
            Sigma = Sigma_original_cell{1,1};
            remove_arr = [1:min_reg_val-1,max_reg_val+1:size(Sigma,1)];
            
            Sigma(remove_arr,remove_arr) = 0;
            
            % reconstruct original matrix
            s_svd = U_cell{1,1}*Sigma*V_cell{1,1}';
            s_svd_all(:,:,:) = reshape(s_svd,size_bf_temp); 
            toc_1 = toc(tic_1);
            disp(['Removal of components and reconstructing matrix for ensemble ' num2str(1) ' took ' num2str(toc_1) ' seconds.'])
            %{
            % Displaying svd results
            figure()
            imshow(s_svd_all(:,:,1),[])
            
            
            figure()
            plot(diag(Sigma_original*Sigma_original')/(size(Sigma_original,1)-1))
            ylim([0 100000])
            xlim([0 500])
            grid on; grid minor
            ylabel("Feature strength \sigma\sigma'/499")
            xlabel("Feature nr.")
            %}
            clearvars bf_temp s_svd




            % Reconstruct the fifth frame to check if you are satisfied with the SVD 
            RF_temp = preproc_RcvFrame(s_svd_all(:,:,5), att_vec,near_field_sat,f_axis,f_c);
            img_svd_test=do_reconstruction(RF_temp, DAS_matrix_sparse, recon_choices, Receive.decimSampleRate*1e6, TW.Parameters(1)*1e6, lx,lz);

            % Show the image and original image
            figure(2);clf(2);
            subplot(1,2,1)
            if ~recon_choices.demod_rf
                frame_orig = abs(hilbert(img_orig))';
                frame_svd = abs(hilbert(img_svd_test))' ;
            else
                frame_orig = abs(img_orig)';
                frame_svd = abs(img_svd_test)';
            end
    
            imagesc(x,z, 20*log10(frame_orig./max(frame_orig, [], 'all' ))); colormap('gray')
            axis equal;xlim([x(1) x(end)]);ylim([z(1) z(end)])
            title(" Original")
            subplot(1,2,2)
            
            imagesc(x,z,20*log10(frame_svd./max(frame_svd, [], 'all' ))); colormap('gray')
            axis equal; xlim([x(1) x(end)]);ylim([z(1) z(end)])
            title(['Image after SVD with min = ' num2str(min_reg_val) ' , max = ' num2str(max_reg_val) ])
            f2 = figure(2);f2.WindowState = 'maximized';
            if start_with_prompt
                % Ask user if now satisfied with result
                prompt = {'Satisfied with this SVD result? (0/1)'};
                dlgtitle = 'Check the results';
                fieldsize = [1 45];
                definput = {'0'};
                answer = inputdlg(prompt,dlgtitle,fieldsize,definput);
                satisfied_with_SVD = str2double(answer{1});
            else
                satisfied_with_SVD = 1;
            end
            if ~satisfied_with_SVD
                start_with_prompt =1;
            end
        end
        ensem_z=1;
        tic_assign= tic;
        RcvData_sorted(:,:,(ensem_z-1)*SVD_ensemblesize+1:ensem_z*SVD_ensemblesize)=s_svd_all;
        toc_assign = toc(tic_assign);
        disp(['Assigning reconstructed data to RcvData_sorted for ensemble ' num2str(ensem_z) ' took ' num2str(toc_assign) ' seconds.'])
        
        if N_ensem>1
            % ========SVD for all ensembles using one set of limits for SVD
            for ensem_z = 2:N_ensem
                tic_1=tic;
                Sigma = Sigma_original_cell{1,ensem_z};
                remove_arr = [1:min_reg_val-1,max_reg_val+1:size(Sigma,1)];
                
                Sigma(remove_arr,remove_arr) = 0;
                
                % reconstruct original matrix
                s_svd = U_cell{1,ensem_z}*Sigma*V_cell{1,ensem_z}';
                s_svd_all(:,:,:) = reshape(s_svd,size_bf_temp); 
                toc_1 = toc(tic_1);
                disp(['Removal of components and reconstructing matrix for ensemble ' num2str(ensem_z) ' took ' num2str(toc_1) ' seconds.'])
                tic_assign= tic;
                RcvData_sorted(:,:,(ensem_z-1)*SVD_ensemblesize+1:ensem_z*SVD_ensemblesize)=s_svd_all;
                toc_assign = toc(tic_assign);
                disp(['Assigning reconstructed data to RcvData_sorted for ensemble ' num2str(ensem_z) ' took ' num2str(toc_assign) ' seconds.'])
                
            end
        end
        clearvars s_svd_all;
        
end
clearvars S_mean U_cell V_cell Sigma_original 

%% Reconstruction for all frames


% Preallocation
images_IQ = cell(1,Nfrs_recon);

% Waitbar for parallel function
D = parallel.pool.DataQueue;
global h
h = waitbar(0, ['Reconstruction progress for file ', num2str(i_file), ' of ', num2str(length(files_info))]);
afterEach(D, @nUpdateWaitbar);
global p
p=1; % Progress counter
global Nfrs_recon

for i_frame = 1:Nfrs_recon

    % ====================================================================
    % Preprocess a single frame of the sorted receive data (see function at
    % end of this file)
    RF_temp = preproc_RcvFrame(RcvData_sorted(:,:,i_frame), att_vec,near_field_sat,f_axis,f_c);
    % =====================================================================
    % Reconstruct frame (see function at end of this file)
    images_IQ{i_frame} = do_reconstruction(RF_temp, DAS_matrix_sparse, recon_choices, Receive.decimSampleRate*1e6, TW.Parameters(1)*1e6, lx,lz);

    % Update progress bar
    send(D, i_frame);
end
toc
close(h); clearvars p


function nUpdateWaitbar(~)
    global p
    global Nfrs_recon
    global h
    waitbar(p/Nfrs_recon, h);
    p = p + 1;
end

function RF_temp = preproc_RcvFrame(RcvData_sorted_frame_i, att_vec,near_field_sat,f_axis,f_cutoff)
    % apply TGC
   RF_temp=double(RcvData_sorted_frame_i.*att_vec); 
    % cut off values of the RF-line that are higher/lower than the
   % near field saturation value specified
   RF_temp(RF_temp > near_field_sat) = near_field_sat;
   RF_temp(RF_temp <- near_field_sat) = -near_field_sat;
   % high-pass filter to get rid of near-field error due to
   % electronics
   RF_temp_FT = fft(RF_temp-mean(RF_temp)); %fourier transform of fluctuations&oscillations
   RF_temp_FT(floor(length(RF_temp_FT)/2):end, :) = 0; %half of the frequencies are set to zero
   RF_temp_FT(f_axis<=f_cutoff, :) = 0; %all values below the cut off frequency are set to zero
   RF_temp = 2.*real(ifft(RF_temp_FT)); %amplitude is doubled since we lost half of the frequencies by taking only the positive/real part
end

function image = do_reconstruction(RF_temp, DAS_matrix_sparse, recon_choices, fs, f_centre_demod, lx,lz)
    % fs (=sampling rate) and f_centre_demod should be in Hz
    if recon_choices.demod_rf && strcmp(recon_choices.demod_rf_method, 'rf_hilb')
        RF_temp = hilbert(RF_temp);
    elseif recon_choices.demod_rf && strcmp(recon_choices.demod_rf_method, 'rf2iq')
        RF_temp = rf2iq(RF_temp,fs, f_centre_demod);
    end
    image = reshape(DAS_matrix_sparse*reshape(RF_temp, [size(RF_temp,1)*size(RF_temp,2),1]), lx,lz);
    if recon_choices.demod_rf
        image = abs(image);
    end
end