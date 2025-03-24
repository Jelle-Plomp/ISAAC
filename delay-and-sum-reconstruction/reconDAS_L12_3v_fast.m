function ImageData = reconDAS_L12_3v_fast(RcvData)
    % This function is used to reconstruct images during live-view and 
    % during the premeasurement (iterative procecure).
    
    % Works with L12-3 (probably also other linear transducers), assuming only 1 aperture is used.
    
    % Assume RcvData has shape (bufferlength * n_elements)
    
    % Jelle Plomp. 2024.
    
    % Only Receive(end) is loaded, meaning that all the values used from
    % the Receive struct should be the same for all fields in the Receive
    % struct.
    

    tic_1=tic;
    % Parameters for reconstruction
    att = 0; % Attenuation coefficient
    use_raddec = true; % Use decay (compensating for radial scattering 1/r)
    max_angle = 20;
    near_field_sat = 400;
    f_c=0.5; %MHz, cut-off frequence


    ProcTimes = evalin('base', 'ProcTimes');
    frame_num = length(ProcTimes)+1;
    ReconInfo_DAS = evalin('base', 'ReconInfo_DAS');
    if ReconInfo_DAS.first_frame
        % Evaluate structs in base work space
        TW = evalin('base','TW');
        Receive = evalin('base', 'Receive(end)');
        Trans = evalin('base', 'Trans');
        Resource = evalin('base', 'Resource');
        nr_elements_receive = Resource.Parameters.numRcvChannels;

        % Order the element data
        Npt = Receive.endSample-Receive.startSample+1;
        order = Trans.Connector(Receive.aperture:Receive.aperture+nr_elements_receive-1,1); % Only use the elements we use in receive
        RcvData_sorted = double(RcvData(1:Npt,order,:));

        % Obtain variables which remain constant over all frames
        ReconVars = reconstruct_DAS_get_vars_sparsematrix(Trans, ReconInfo_DAS, Receive, Resource, TW, size(RcvData_sorted(:, :,1)), max_angle, use_raddec, att,(ReconInfo_DAS.demod_rf && strcmp(ReconInfo_DAS.demod_rf_method, 'rf2iq')));
        % Add order to ReconVars and asssign to workspace
        ReconVars.order = order;
        assignin('base', 'ReconVars', ReconVars)
        ReconInfo_DAS.first_frame = false;
        assignin('base', 'ReconInfo_DAS', ReconInfo_DAS)
    else
        ReconVars = evalin('base', 'ReconVars');
        % Order the element data
        RcvData_sorted = double(RcvData(1:ReconVars.Npt,ReconVars.order,:));

    end
    

    %% Data filtering and image reconstruction
    % Preprocess the frame of the sorted receive data (see function at
    % end of this file)
    tic_preproc = tic;
    RF_temp = preproc_RcvFrame(RcvData_sorted(:,:,1), ReconVars.att_vec',near_field_sat,ReconVars.f_axis,f_c);
    ProcTimes(frame_num).RFFilter = toc(tic_preproc);
    % =====================================================================
    % Reconstruct frame (see function at end of this file)
    tic_demod_and_DAS = tic;
    ImageData = do_reconstruction(RF_temp, ReconVars.DAS_matrix_sparse, ReconInfo_DAS, ReconVars.Fs, ReconVars.F_centre_true, ReconVars.lx,ReconVars.lz);
    ProcTimes(frame_num).time_demod_and_DAS = toc(tic_demod_and_DAS);

 
    %% Processing time 
    
    ProcTimes(frame_num).time_in_reconDAS = toc(tic_1);
    % Total time between frames (including acquisition)
    if frame_num>2
        ProcTimes(frame_num).framerate = toc(ProcTimes(frame_num-1).endtime);
    end
    ProcTimes(frame_num).endtime = tic;
    
    assignin('base', 'ProcTimes', ProcTimes);
    
end

function RF_temp = preproc_RcvFrame(RcvData_sorted_frame, att_vec,near_field_sat,f_axis,f_cutoff)
    % apply TGC
   RF_temp=double(RcvData_sorted_frame.*att_vec); 
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
        image = abs(image)';
    else
        image = abs(hilbert(image'));
    end
end