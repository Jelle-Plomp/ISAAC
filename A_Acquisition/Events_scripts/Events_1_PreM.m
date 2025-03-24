% Event structure for the premeasurement.
%{ 
    
%}
for i = 1:PreM.ES.Nit_max
    if i==1 % Dummy acquisition
        % A single dummy acquisition is included to prevent artefacts (that
        % would otherwise be present in the first frame).
        Event(n).info = 'PreMeasurement: acquire dummy aperture';
        Event(n).tx = 1;
        Event(n).rcv = i+rcvNum_prem_start; 
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = [2,nsc];
            SeqControl(nsc).command = 'transferToHost'; % nsc count used because there must be a unique SeqControl structure for each 'transferToHost' command in the sequence
            nsc = nsc + 1;                               
        n = n+1;
    else
           
        
        Event(n).info = 'Sync'; 
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0; 
        Event(n).seqControl = 7; 
        n=n+1; 

        Event(n).info = 'PreMeasurement: acquire (first) aperture';
        Event(n).tx = 1;
        Event(n).rcv = i+rcvNum_prem_start-1; % Always put in buffer 2. -1 such that the "dummy" frame (for i==1) is overwritten 
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2; % Time to next acquisition
        n = n+1;
      
        % Transfer to host after acquiring the frame
        Event(n-1).seqControl = [2, nsc]; % Time to next acq, transfertohost
        SeqControl(nsc).command = 'transferToHost'; % nsc count used because there must be a unique SeqControl structure for each 'transferToHost' command in the sequence
           nsc = nsc + 1;                              
        
        Event(n).info = 'Sync'; % To ensure we actually reconstruct the last acquired frame
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0; 
        Event(n).seqControl = 7; 
        n=n+1; 
        
        Event(n).info = 'Reconstruct';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0; 
        Event(n).process = P_nr_DASpreM; % our own DAS reconstruction
        Event(n).seqControl = 0;
        n = n+1;

        Event(n).info = 'Call PreMeasurement';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = P_nr_PreM;  % Call Premeasurement main function
        Event(n).seqControl = [3,5]; % return to matlab, and (just to be sure the change applies) reset the TPC to TPC(1)
        % Note: this setting the TPC profile was a remainder of a version
        % of the algorithm where the voltage was also adapted. For changing
        % the apodisation, the TPC probably does not need to be set.
        n=n+1; 

        Event(n).info = 'Wait';
        Event(n).tx = 0;
        Event(n).rcv = 0;
        Event(n).recon = 0;
        Event(n).process = 0; 
        Event(n).seqControl = 10; % Wait 10 milliseconds
        n=n+1; 
    end
end