% Events: Saving HFR for all different TX.
% When the partial compensation is on, this means there will be 7 HFR
% recordings:
% TX1: A_ISAAC
% TX2: A_0
% TX3 to TX7: A_pc (20,40,50,60,80%)

% Written by Jelle Plomp.
%% Part 1: acquisition TX1 and current voltage
Event(n).info = 'start timer';
Event(n).tx = 0;
Event(n).rcv =0;
Event(n).recon = 0;
Event(n).process = P_nr_timer_end; 
Event(n).seqControl = 0;
n = n+1;

for TX_i = 1:length(TX)
    % Acquire all frames defined in RcvBuffer
    for i = 1:P.HFRBufFrames
        for j = 1:P.numAcqsSuperFrame
            Event(n).info = ['HFR - acquisition TX ', num2str(TX_i)];
            Event(n).tx = TX_i;
            Event(n).rcv = PreM.PartialCompensation.RcvStartNums(TX_i)+(i-1)*P.numAcqsSuperFrame+j; 
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
end


%% Saving frames from all buffers

for TX_i = 1:length(TX)
    Process_num_saveTX_i = PreM.PartialCompensation.ProcessNumbers(TX_i);
    Bufnum_TX_i = PreM.PartialCompensation.BufferNumbers(TX_i);
    Event(n).info = ['Save buffer ', num2str(Bufnum_TX_i), ' data to File (TX', num2str(TX_i), ')']; 
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % no reconstruction
    Event(n).process = Process_num_saveTX_i;    % external processing - SaveRcvData
    Event(n).seqControl = 8; % wait 
    n = n+1;

end
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