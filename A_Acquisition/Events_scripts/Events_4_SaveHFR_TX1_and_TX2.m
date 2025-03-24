% Events: Acquisition and Saving HFR (first using TX1, then using TX 2)
% The two acquisition sequences are first completed, after which the data
% is saved. 

% Written by Jelle Plomp.
%% Part 1: acquisition TX1 
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
        Event(n).info = 'HFR - acquisition TX1 V_cur.';
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


%% Part 2: acquisition TX2
% Acquire all frames defined in RcvBuffer
for i = 1:P.HFRBufFrames
    for j = 1:P.numAcqsSuperFrame
        Event(n).info = 'HFR - acquisition TX2 V_cur.';
        Event(n).tx = 2;
        Event(n).rcv = rcvNum_HFR_D2_TX2_Vend+(i-1)*P.numAcqsSuperFrame+j;
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

%% Saving frames from both buffers

Event(n).info = 'Save buffer 3 data to File'; 
Event(n).tx = 0;         % no transmit
Event(n).rcv = 0;        % no rcv
Event(n).recon = 0;      % no reconstruction
Event(n).process = P_nr_saveHFR_D1_TX1_Vend;    % external processing - SaveRcvData
Event(n).seqControl = 8; % wait 
n = n+1;


Event(n).info = 'Save buffer 4 data to File'; 
Event(n).tx = 0;         % no transmit
Event(n).rcv = 0;        % no rcv
Event(n).recon = 0;      % no reconstruction
Event(n).process = P_nr_saveHFR_D2_TX2_Vend;    % external processing - SaveRcvData
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
Event(n).process = [P_nr_resetStartEvent, 3] ;   % Set&Run and then return to matlab. After going back into runAcq, the set&run will be applied.
Event(n).seqControl = 0;
n = n+1;

Event(n).info = 'Jump back to first event'; % Not sure if this will be called after the Reset StartEvent (but then we already jump back so that would not be a problem)
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 0; 
Event(n).seqControl = 4; % jump back to the begining of live view
n=n+1;