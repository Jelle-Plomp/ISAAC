%% Events:Live view
% Events for a live view that uses a reconstruction defined in an external
% process.
for i = 1:Resource.RcvBuffer(1).numFrames
    Event(n).info = 'Live: acquire aperture';
    Event(n).tx = 1;
    Event(n).rcv = P.na*(i-1) + 1; % Corresponds to indices set in Receive
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = [2,nsc];
       SeqControl(nsc).command = 'transferToHost'; % nsc count used because there must be a unique SeqControl structure for each 'transferToHost' command in the sequence
       nsc = nsc + 1;
    n = n+1;

    Event(n).info = 'Reconstruct';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = P_nr_DASlive; % Our own DAS reconstruction
    Event(n).seqControl = 0; 

    n = n+1;
    Event(n).info = 'Display';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = P_nr_display; % Display the reconstructed image
    Event(n).seqControl = 0;  
    
    if floor(i/3) == i/3     % Return to Matlab every 3rd frame
        %(only after a return to matlab, the program will respond to the change in GUI button)
        Event(n).seqControl = 3;
    end

    n = n+1;
end

% end
Event(n).info = 'Jump back to first event';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 4;
n=n+1;