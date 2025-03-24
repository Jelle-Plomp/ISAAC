function PreM_main(ImData)
    % PREM_MAIN describes the steps of the iterative procedure for active
    % attenuation correction (ISAAC). 

    %{

        Jelle Plomp
        Documentation updated 01-02-2023.

        Input ImData will be an image (before log 
        compressing it), the input is defined in the
        Setup via Process().Parameters{'srcbuffer', } where srcbuffer can
        be 'image' or 'imageP' (both have the same size)
            (These refer to ImageBuffer and ImagePBuffer)

        Outputs:
            History: will contain for each image in the image buffer (of
            the premeasurement) some information, such as the apodisation,
            voltage, error, etc. 

            TPC and TX will be updated and both assigned to the matlab
            workspace as well as using set&run/update&run to also apply the
            changes in the acquisition hardware.

            ProcTimes struct is updated with the time spent in this
            function.
        
    %} 
    
    tic_PreM=tic;
    
    % Evaluate variables from the PreM struct 
    PreM = evalin('base', 'PreM');
   
    %% Obtain the ROI (draw or get it from workspace or file)
    doesExist = isfield(PreM,'ROI');
    
    if ~doesExist % Case: ROI does not exist in PreM struct (i.e. in the first iteration)
     
        PreM.displayInfo = evalin('base','displayInfo');
        % Ask if ROI should be loaded from file
        load_from_file = inputdlg("Do you want to load an existing ROI from file? (true/false)");
        if strcmp(load_from_file{1}, 'true')
            P = evalin('base','P');
            [fname,dirname] = uigetfile({'*_Mask.mat'}, 'Select a mask file', P.savefiledir);
            disp('Loading ROI from file');
            load(fullfile(dirname, fname),'ROI', 'poly');
            PreM.ROI = ROI;
            PreM.poly = poly;
            P.filename_ROImask = fname; 
            assignin('base', 'P', P)
            
        else 
            
            % Ask user to draw ROI
           [PreM.ROI, PreM.poly] = PreM_drawROI(PreM.x,PreM.z,20*log10(ImData/max(max(ImData))), PreM.displayInfo.dynRange, 1);

          
        end
        if ~isfield(PreM, 'xel')
            PreM.xel = evalin('base', 'ReconVars.xel');
            PreM.zel = evalin('base', 'ReconVars.zel');
        end
        if strcmp(PreM.method, 'rayleigh')
            % Initialise rayleigh integral
            Trans = evalin('base','Trans');
            TW = evalin('base','TW');
            PreM.rayleigh.waveform_250MHz = TW(1).Wvfm1Wy;
            line_in_el_x_range = PreM.ROI.x>=PreM.xel(1) & PreM.ROI.x<=PreM.xel(end);
        
            PreM.rayleigh.x_bounded = PreM.ROI.x(line_in_el_x_range==1);
            [PreM.rayleigh.F, PreM.rayleigh.P0, PreM.rayleigh.P1, PreM.rayleigh.P0_min_range, PreM.rayleigh.P1_max_range] = PreM_v_rayleigh_initialisation(PreM.rayleigh.x_bounded, mean(PreM.ROI.z).*ones(1,length(PreM.rayleigh.x_bounded)), PreM.xel, PreM.zel, PreM.ReconInfo_DAS.c, Trans, PreM.rayleigh.waveform_250MHz);
        else
            error('No other attenuation correction models were included in this code version.')
        end
        assignin('base',"PreM", PreM) % so that the ROI and poly are added in base as well.
    end
   
       
    %% Obtain the relevant data (average lineprofile) from the ROI
    ROI_Data = ImData(PreM.ROI.Position_pix(2):PreM.ROI.Position_pix(2)+PreM.ROI.Position_pix(4), PreM.ROI.Position_pix(1):PreM.ROI.Position_pix(1)+PreM.ROI.Position_pix(3));
    ROI_averageline = mean(ROI_Data,1);
    ROI_max = max(ROI_averageline);
    % Some parameters 
    % PreM.PI_contrl.std_mm  indicates the standard deviation of the filter used in mm 0.5;
    std_datapoints = (PreM.std_mm*1e-3)/PreM.ReconInfo_DAS.pas;
    
    % Get a smooth profile from the ROI average line
    ROI_smooth= imgaussfilt(ROI_averageline, std_datapoints);
    
    if strcmp(PreM.method, 'rayleigh')
        line_in_el_x_range = PreM.ROI.x>=PreM.xel(1) & PreM.ROI.x<=PreM.xel(end);
        
        I_L_used = ROI_smooth(line_in_el_x_range==1);

    end
    %% History
    % Get the current apodization
    TX = evalin('base','TX');
    apod_current = TX(1).Apod;
    % Get current Voltage from TPC.hv
    TPC = evalin('base', 'TPC');
    Voltage_current = TPC(1).hv;
    % Save the current apodization and Voltage (as a reference when
    % studying the results of the measurements)
    History = evalin('base','History'); % This struct contains all the values such that we can look back at it later
    History(end) = struct('Voltage', Voltage_current, 'Apod', apod_current,...
        'Imax_ROI', ROI_max, 'I_ROI_line', ROI_averageline, 'I_L_used', I_L_used,...
        'P', NaN, 'I', NaN, 'I_clip',NaN,'error_total', NaN, 'el_clip', NaN, 'Nit_stable', NaN);
    
    
    %% Set the SetPoint for the PI controller
    if length(History)==1 % In the first iteration
        if strcmp(PreM.PI_contrl.sp_version, 'median_high10perc')
            % Obtain highest 10% of values in the I_L_used
            I_L_used_sorted = sort(I_L_used, 'descend');
            ROI_highest_10perc = I_L_used_sorted(1:ceil(length(I_L_used)*0.1));
            PreM.PI_contrl.sp = ones(1,length(I_L_used)).*median(ROI_highest_10perc); % Set the set point as the median of the 10% highest values.
        else
            error("No valid method assigned in PreM.PI_contrl.sp_version")
        end
        assignin('base',"PreM", PreM)
    end
    %% Check if current image is satisfactory
    if PreM.usestopping
        error('A stopping criterium has not yet been defined but can be implemented in PreM_checkCondition. Set PreM.usestopping to false to continue without it.')
        if length(History)==1
%             I_L_used_prev = NaN;
            apod_prev = NaN;
            Nit_stable_in = 0;
            
        else
%             I_L_used_prev = History(end-1).I_L_used;
            apod_prev = History(end-1).Apod;
            Nit_stable_in=History(end-1).Nit_stable;
        end

        [continue_iters,History(end).Nit_stable] = PreM_checkCondition(apod_current, apod_prev,  Nit_stable_in, PreM);

    elseif PreM.ES.Nit_max == length(History)
        continue_iters = false;
    else
        continue_iters = true;
    end
    %% Calculate new apodization 
    if continue_iters
        if strcmp(PreM.method, 'rayleigh')
            if strcmp(PreM.PI_contrl.errormethod,'clipintegralerror')
                [apod_new, History, pressure_req] = PreM_v_rayleigh_calcApod_clipintegralerror(PreM, History);
            else
                error("With PreM.method 'rayleigh', PreM.PI_contrl.errormethod can only be set to 'clipintegralerror'")
            end

        end
        
        % Assigning outputs to base workspace for next acquisition
        TX(1).Apod = apod_new;
        assignin('base', 'TX', TX)

    end
    
    %% Plotting
    if ~continue_iters
        apod_new_plot = nan;
        pressure_req = nan;
    else
        apod_new_plot = apod_new;
    end
    if doesExist
        ROI_averageline_it1 = History(1).I_ROI_line;
    else
        ROI_averageline_it1 = NaN;
    end
    % Calculate history statistics if needed for plotting
    if ~isfield(PreM, 'GUI') || strcmp(PreM.GUI, 'default') 
        History_stats = NaN;
    elseif isfield(PreM, 'GUI') && strcmp(PreM.GUI, 'extended')
        it = length(History);
        if it==1
            % Initialise the fields and assign to base. Since in the first
            % iteration, there is nothing to compare to, no calculations
            % are performed yet.
          
            History_stats.I_change_rel = nan(1, PreM.ES.Nit_max);
            History_stats.I_change_rel(1) = 0;
            History_stats.Apod_change = nan(1, PreM.ES.Nit_max);
            History_stats.Apod_change(1) = 0;
            History_stats.E_cur_mean = nan(1, PreM.ES.Nit_max);
            if strcmp(PreM.method, 'rayleigh')
                History_stats.I_clip_max = nan(1, PreM.ES.Nit_max);
    %             History_stats.E_tot_mean_abs = zeros(1, PreM.ES.Nit_max);
                History_stats.I_clip_mean = nan(1, PreM.ES.Nit_max);
            end
        else
            History_stats = evalin('base', 'History_stats');
            History_stats.I_change_rel(it) = mean(abs((History(end).I_L_used-History(end-1).I_L_used)./PreM.PI_contrl.sp));
            History_stats.Apod_change(it) = mean(abs((History(end).Apod-History(end-1).Apod)));
        end
        if strcmp(PreM.method, 'rayleigh')
            History_stats.I_clip_max(it) = max(abs(History(end).I_clip));
            History_stats.I_clip_mean(it) = mean(History(end).I_clip);
        end
            error_cur = (PreM.PI_contrl.sp - History(end).I_L_used)./PreM.PI_contrl.sp;
        History_stats.E_cur_mean(it) = mean(error_cur);
        assignin('base', "History_stats", History_stats)
    end
    
    if ~isfield(PreM, 'GUI') || strcmp(PreM.GUI, 'default') 
        PreM_GUI_default(PreM, ImData, ROI_averageline,  ROI_averageline_it1,  apod_new_plot, continue_iters, doesExist)

    elseif isfield(PreM, 'GUI') && strcmp(PreM.GUI, 'extended') 
         PreM_GUI_extended(PreM, ImData, ROI_averageline,  ROI_averageline_it1,  apod_new_plot, pressure_req, continue_iters, doesExist, History_stats)

    end
    
    %% Initialise next row in History struct, and assing to base
    if continue_iters
        History(end+1) = struct('Voltage', NaN, 'Apod', NaN, 'Imax_ROI', NaN, ...
            'I_ROI_line', NaN, 'I_L_used', NaN,...
            'P', NaN, 'I', NaN,'I_clip', NaN, 'error_total', NaN, 'el_clip', NaN, 'Nit_stable', NaN); 
    end
    assignin('base', 'History', History)
    %% update&Run
    % To actually use the new TX, we need to update it using control
    % This control_i prevents overwriting any present control arguments (if
    % this is not done, then some of the GUI buttons may fail)
    Control = evalin('base', 'Control');
    if isempty(Control.Command)
        control_i = 1;
    else
        control_i =2;
    end

    if continue_iters
                  
        Control(control_i).Command = 'update&Run';
        Control(control_i).Parameters = {'TX'};
    else  % Stop iterations
        
        Resource = evalin('base', 'Resource');
        event_nr_live = evalin('base', 'event_nr_live');
        current_startevent = Resource.Parameters.startEvent;
        if current_startevent == 1
            Resource.Parameters.startEvent = event_nr_live;
        else
            error("Something must have gone wrong with assigning startEvent" )
        end
        assignin('base', 'Resource', Resource);
        
        Control(control_i).Command = 'update&Run';
        Control(control_i).Parameters = {'Parameters'};
        control_i = control_i+1;
        if PreM.enablePartialCompensation
           % Calculate TX arrays for Partial compensation
           % Note: in case of manual stop, this is called in the callback
           % function (see step1_MAIN_acquisition_HFR_L12_3v.m)
           if strcmp(PreM.method, 'rayleigh')
                TX = calcTX_PartialCompensation_v_rayleigh(TX);
           end
           assignin('base', 'TX', TX)
           Control(control_i).Command = 'update&Run';
           Control(control_i).Parameters = {'TX'};
        end
        
    end
    assignin('base','Control', Control);
    
    ProcTimes = evalin('base', 'ProcTimes');
    ProcTimes(end).PreM_main = toc(tic_PreM);
    assignin('base', 'ProcTimes', ProcTimes);
    
end
