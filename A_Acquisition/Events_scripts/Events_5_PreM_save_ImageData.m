% Events to save the imagedata and other history of the premeasurement
% (PreM), which is the iterative procedure.

% Written by Jelle Plomp.

%We don't need a 'waitForTransferComplete' command before saving, since in the 
% premeasurement upon each frame the data is being transferred anyway.


Event(n).info = 'Save Data to File'; 
Event(n).tx = 0;         % no transmit
Event(n).rcv = 0;        % no rcv
Event(n).recon = 0;      % reconstruction
Event(n).process = P_nr_PreMsave;    % external processing - SaveImageData
Event(n).seqControl = 8; % wait 
n = n+1;

Event(n).info = 'Sync'; 
Event(n).tx = 0;         % no transmit
Event(n).rcv = 0;        % no rcv
Event(n).recon = 0;      % reconstruction
Event(n).process = 0;    % external processing - SaveImageData
Event(n).seqControl = 7; % sync
n = n+1;

Event(n).info = 'Reset StartEvent'; 
Event(n).tx = 0;         
Event(n).rcv = 0;        
Event(n).recon = 0;      
Event(n).process = [P_nr_resetStartEvent, 3] ;   
Event(n).seqControl = 0;
n = n+1;

Event(n).info = 'Jump back to live view'; % This event is probably not reached since the previous event resets the start event.
Event(n).tx = 0;         % no transmit
Event(n).rcv = 0;        % no rcv
Event(n).recon = 0;      % reconstruction
Event(n).process = 0;    
Event(n).seqControl = 4; % return to beginning of live view
n = n+1;