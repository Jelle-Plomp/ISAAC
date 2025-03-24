function [Apod, History_updated, pressure_req] = PreM_v_rayleigh_calcApod_clipintegralerror(PreM, History)
    % PreM_v_rayleigh_calcApod_clipintegralerror Version 13-05-2024: A PI
    % controller that clips the integral error.

    % Apod is the resulting apodization (value per Transducer element)
    
    % PreM.xel = element x positions (m)
    % PreM.PI_contrl = Proportional integral controller parameters
    % PreM.std_mm standard deviation used to smooth apodisation

    plottrue = false;
  
    I_L_used = History(end).I_L_used;

    % Getting total error of previous iteration
    if length(History)==1 % In the first iteration
        I_prev = zeros(1,length(PreM.rayleigh.x_bounded));
        I_noclip_prev = zeros(1,length(PreM.rayleigh.x_bounded));
    else
        I_prev = History(end-1).I_clip;
        I_noclip_prev = History(end-1).I;
    end
    % Error in current iteration
    error_i = (PreM.PI_contrl.sp - I_L_used)./PreM.PI_contrl.sp;
    History(end).P = PreM.PI_contrl.Kp.*error_i;
    
    % Just for reference, we also see what the integral error would have
    % been without ever clipping it
    History(end).I = I_noclip_prev + PreM.PI_contrl.Ki.*error_i;
    % Error accumulation
    
    
    I_cur = I_prev + PreM.PI_contrl.Ki.*error_i;
    lim_low = zeros(1,length(PreM.rayleigh.P0));
    lim_up = PreM.rayleigh.P1_max_range-PreM.rayleigh.P0;
    History(end).I_clip = min(max(I_cur,lim_low),lim_up);
    History(end).el_clip = (I_cur~=History(end).I_clip); % Store for which elements the signal is clipped (just for reference, not used in any calculations)
    
    % Requested pressure
    pressure_req = History(end).P + History(end).I_clip + PreM.rayleigh.P0;

    % Calculate apodisation that would be needed to obtain the requested 
    % pressure. 

    % The value F*An - pressure_req is minimised, as well as weight.*(An-imgaussfilt(An,std_el))
   
    apod_lowerlimit = 0.2;
    n_el = length(PreM.xel);
    el_pitch = (PreM.xel(2)-PreM.xel(1))*1000;
    std_el = PreM.std_mm/el_pitch; % Standard deviation in nr. elements.
    myfun_rmse =  @(An, F) horzcat(An*F - pressure_req, PreM.rayleigh.weight.*(An-imgaussfilt(An,std_el)));
    optim_targetvals_rmse = horzcat(zeros(1,length(PreM.rayleigh.x_bounded)), zeros(1,n_el));
    % Curve fitting. A boundary of 0.2 to 1 is opposed on An_out 
    % (apod can only be between 0.2 and 1).
    An_out = lsqcurvefit(myfun_rmse, History(end).Apod, PreM.rayleigh.F, ...
        optim_targetvals_rmse,apod_lowerlimit*ones(1,n_el), ones(1,n_el));
    
    % Just to ensure apod is actually clipped between 0.2 and 1, we clip
    % once more. But should not change anything.
    Apod = min(max(An_out, apod_lowerlimit), 1); %
    History_updated = History;

    %% Plotting
    if plottrue
        plotprofiles(1).profile = control_signal_clipped; plotprofiles(1).name = 'control signal S_c (after clip)';
        plotprofiles(2).profile = Apod_raw; plotprofiles(2).name = 'Apod_{raw} = Apod_init + S_c';
        plotprofiles(3).profile = Apod_smooth;  plotprofiles(3).name = 'Apod smooth';
    
        colour = 'b';
        figure(102)
        for i=1:3
            subplot(3,1,i)
            plot(PreM.xel.*1e3,plotprofiles(i).profile, 'Color', colour)
            xlim([-15 15])
            title(plotprofiles(i).name)
            if i~=1
                ylim([0 1])
            end
        end
        figure(101)
    end
             
    
end
