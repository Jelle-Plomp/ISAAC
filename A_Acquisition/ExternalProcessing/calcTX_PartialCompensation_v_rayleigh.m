function TX = calcTX_PartialCompensation_v_rayleigh(TX)
    % CALCTX_PARTIALCOMPENSATION_v_rayleigh uses the input profile TX1.Apod
    % It assumes a linear relation between apodisation and pressure,
    % meaning we can just half the (added) apodisation (=apodisaiton-0.2)
    % to assume that we get half the extra pressure.
    
    
    % This method assumes that the value of 0.2 is the baseline value.
    % TX2 is assumed to be the profile wiht 0.2.
    if sum(TX(2).Apod==0.2) ~= 128
        error("calcTX_PartialCompensation requires the Apodisation profile to be initialised with value of 0.2" )
    end
    
    %
    ratio = [0.2 0.4 0.5 0.6 0.8];
    if length(TX) ~= 7
        error("Current method assumes length of TX to be 7")
    end
    
     for j = 1:5  
        Apod_minus_baseline = TX(1).Apod-0.2;

        Apod = ratio(j).*Apod_minus_baseline + 0.2;

        TX(j+2).Apod = min(max(Apod, 0.2), 1);
    end

end