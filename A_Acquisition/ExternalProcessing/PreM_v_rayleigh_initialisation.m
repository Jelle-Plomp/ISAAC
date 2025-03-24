function [F, P0, P1, P0_min_range, P1_max_range] = PreM_v_rayleigh_initialisation(xi, zi, xel, zel, recon_c,Trans, waveform250MHz)
        % PreM_v_rayleigh_initialisation calculates the matrix F from the
        % article (Imaging Behind the Plaque: Improved Blood Flow
        % Quantification Using an Iterative Scheme for Active Attenuation
        % Correction - UMB 2025), as well as the pressure profiles obtained
        % with all apodization values at 0.2 and all at 1.0.

        % Jelle Plomp. 2024.
        
        % Get time to waveform peak (peak positive pressure)
        t_wave0 = [0:length(waveform250MHz)-1]./(250*1e6);
        [~,t_idx_max] = max(waveform250MHz);
        t_peak = t_wave0(t_idx_max);
        % Define the system that describes our problem
        T = zeros(1,length(xi)); % time that wave arrives at each queried pixel
        D = zeros(length(xi),length(xel)); % for each pixel, the time to each element
        F = zeros(length(xi),length(xel));
        for i=1:length(xi)
            [T(i), D(i,:)] = calc_ti_d_iel(xi(i), zi(i), xel, zel, recon_c, t_peak, Trans, 'linear');
            F(i,:) = get_wave_amp(T(i),D(i,:),recon_c, waveform250MHz);
        end
        % Transpose
        F = F';
        % Calculate the pressure obtained with all apodisations at 0.2
        P0 = (0.2*ones(1,length(xel)))*F;
        % Calculate pressure obtained with all apodisations at 1
        P1 = ones(1,length(xel))*F;
        % Calculate minimum of P0 within x range -10 to 10mm
        P0_min_range= min(P0(xi>-0.01&xi<0.01));
        % Calculate maximum of P1 within x range -10 to 10 mm
        P1_max_range= max(P1(xi>-0.01&xi<0.01));
        
end
            
function [ti,d_iel, dt_peak]=calc_ti_d_iel(xi,zi, xel,zel,c, dt_peak, Trans, linear_or_curved)
    
    % Calculate travel times from all elements to position (xi, zi)
    % Output will be an array of the same size as the number of elements
    plot_on = false;
    d_iel = sqrt((xi-xel).^2+(zi-zel).^2); % distance from point to each element
    
    if strcmp(linear_or_curved, 'linear')
        if length(unique(zel))>1
            error("For a linear array, all element z coordinates should be equal")
        end
        ti = zi(1)./c+  dt_peak ; 
        % A lens correction is ommitted here. If it were to be included, 
        % it would also need to be included in the travel time d_in/c in 
        % t_query (function get_wave_amp), so then they would cancel out.
    elseif strcmp(linear_or_curved, 'curved')
        error('No implementation for curved arrays yet')
    end
    if plot_on
        figure(101);clf(101); plot(d_iel)
    end
end

function amp_out=get_wave_amp(ti, din_array, c, waveform250MHz)
    plot_on = false;
    % Get wave amplitude contribution from an array of elements that fires
    % waveform250MHz, assuming spherical emission from the elements.
    % ti: query time, i.e. the time at which it is assumed that the
    % wavefront arrives at this location (follows from DAS, is a single time). 
    % The ti should already include a time correction factor, just like
    % with DAS.
    
    % waveform250MHz sampled at 250 MHz

    % din_array : distance from each element to point of interest
    % c: reconstruction speed of sound

    t_wave0 = [0:length(waveform250MHz)-1]./(250*1e6);
    t_wave0_spline = [t_wave0(1):1/(501*1e6):t_wave0(end)];
    waveform_spline = spline(t_wave0, waveform250MHz, t_wave0_spline);

    t_query = -din_array/c + ti; 

    % if ti = din_array/c, t_query will be zero, so we get the first value
    % in the waveform. If ti = dt_peak + dt_lenscorr, we should get the maximum amplitude from the closest element. 
    [~,t_idx] = min(abs(t_wave0_spline-t_query'),[],2); % Get minimum along second dimension 
 
    % if t_idx==1 || t_idx==length(t_wave0), then the output will also be 0
    % since the first and last values of waveform250MHz should be zero.
    amp_out = waveform_spline(t_idx)./din_array; % compensate for radial decay
    
    if plot_on
        figure(7)
        subplot(2,2,1);plot(t_query); title("t_query=ti-din_array/c");xlabel("Element n"); ylabel("Distance (m)")
        hold on; plot([1 128],[t_wave0_spline(1) t_wave0_spline(1)],'r'); plot([1 128],[t_wave0_spline(end) t_wave0_spline(end)],'r');legend("t_query", "t_wave0_spline limits")
        subplot(2,2,2);hold on; plot(t_wave0_spline(t_idx)); title("Time in waveform to be sampled");xlabel("Element n"); ylabel("Time (\mus)")
    end 

end