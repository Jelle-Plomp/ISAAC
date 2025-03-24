%% Events: Saving HFR events
% Does HFR acquisition using the apodisation stored in TX(1)
% Written by Jelle Plomp.

Event(n).info = 'start timer';
Event(n).tx = 0;
Event(n).rcv =0;
Event(n).recon = 0;
Event(n).process = P_nr_timer_end; 
Event(n).seqControl = 0;
n = n+1;

% Acquire all frames defined in RcvBuffer
for i = 1:P.HFRBufFrames
    for j = 1:P.numAcqsSuperFrame
        Event(n).info = 'HFR - acquisition TX1 V_c.';
        Event(n).tx = 1;
        Event(n).rcv = rcvNum_HFR_start+(i-1)*P.numAcqsSuperFrame+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 6; %Time to next acqusition 
        n = n+1;
    end
    % Set last acquisitions SeqControl for transferToHost.
    Event(n-1).seqControl = [6,nsc];
        SeqControl(nsc).command = 'transferToHost'; % transfer all acqs in one super frame
        nsc = nsc + 1;
end

Event(n).info = 'wait for final transfer';
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 0;   
Event(n).seqControl = nsc; % wait command
SeqControl(nsc).command = 'waitForTransferComplete';  %
SeqControl(nsc).argument = nsc-1;
nsc=nsc+1;
n = n+1;

Event(n).info = 'Save Data to File'; 
Event(n).tx = 0;         % no transmit
Event(n).rcv = 0;        % no rcv
Event(n).recon = 0;      % no reconstruction
Event(n).process = P_nr_saveHFR;    % external processing - SaveRcvData
Event(n).seqControl = 8; % wait 
n = n+1;

Event(n).info = 'Sync'; 
Event(n).tx = 0;         % no transmit
Event(n).rcv = 0;        % no rcv
Event(n).recon = 0;      % reconstruction
Event(n).process = 0;    % no external processing
Event(n).seqControl = 7; % sync
n = n+1;

Event(n).info = 'end timer';
Event(n).tx = 0;
Event(n).rcv =0;
Event(n).recon = 0;
Event(n).process = P_nr_timer_end;
Event(n).seqControl = 0;
n = n+1;

Event(n).info = 'Reset StartEvent'; 
Event(n).tx = 0;         
Event(n).rcv = 0;        
Event(n).recon = 0;      
Event(n).process = [P_nr_resetStartEvent, 3] ;   % Set&Run and then return to matlab. I think after going back into runAcq, the set&run will be applied.
Event(n).seqControl = 0;
n = n+1;

Event(n).info = 'Jump back to first event'; % Not sure if this will be called after the Reset StartEvent (but then we already jump back so that would not be a problem)
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 0; 
Event(n).seqControl = 4; % jump back to the begining of live view
n=n+1;