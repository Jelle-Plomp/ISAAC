function  SaveRcvDataHFR(RcvData)
    % Function to save the received data.
    % Jelle Plomp 2024.
    tic
    
    % Define parameters to save
    Resource = evalin('base','Resource');
    P = evalin('base','P'); 
    TPC = evalin('base','TPC');
    Trans = evalin('base','Trans'); 
    PData = evalin('base','PData');
    Receive = evalin('base','Receive(end)'); % only saving last Receive (instead of thousands) makes it much faster
    TX = evalin('base','TX');
    TW = evalin('base','TW');
    ReconInfo_DAS = evalin('base','ReconInfo_DAS');
    displayInfo = evalin('base', 'displayInfo');
    TGC = evalin('base','TGC');
    P.ProbeVoltage = TPC(1).hv; % TPC 2 = HFR saving (should only be the case when multiple TPC profiles are created, which we do not do)
    
    pause(0.2) % Pause to wait for RCVdata to transfer to our function.
    %% Check that the re-ordered RcvData does not contain zero frames
    % There were some issues in previous version with saving the data. Data 
    % saving would already be started when the RcvData had not yet been fully 
    % transferred to this function, resulting in loss of data. This
    % was resolved by including the pause statement above. The checks below
    % are a remainder from the debugging process.
    num_empty_superframes = 0;
    num_partly_empty_superframes = 0;
    for i_super = 1:P.HFRBufFrames
        maxima_Rcv_per_channel = max(abs(RcvData(:,:,i_super)));
        idx_zero_per_channel = maxima_Rcv_per_channel == 0;
        num_zeros_in_superframe = sum(idx_zero_per_channel);
        if num_zeros_in_superframe==128
            num_empty_superframes = num_empty_superframes+1;
        elseif num_zeros_in_superframe>0
            num_partly_empty_superframes = num_partly_empty_superframes +1;
        end
    end
    if num_empty_superframes>0 || num_partly_empty_superframes>0
        warning(strcat("There are ", num2str(num_empty_superframes), ...
            " completely empty superframes and ", num2str(num_partly_empty_superframes),...
            " partly empty superframes in the RcvData for this acquisition"))
    end
    %% select only usefull part of RcvData to save and sort it
    RcvData = RcvData(1:Receive.endSample,:,:);
    order = Trans.Connector(Receive.aperture:Receive.aperture+Resource.Parameters.numRcvChannels-1,1); %[1:size(RcvData,2)]'; % Trans.Connector;
    RcvData = RcvData(:,order,:); % Assuming we have frames stacked in 3rd dim
    
    %% Check that the re-ordered RcvData does not contain zero frames
    num_empty_superframes = 0;
    num_partly_empty_superframes = 0;
    for i_super = 1:P.HFRBufFrames
        maxima_Rcv_per_channel = max(abs(RcvData(:,:,i_super)));
        idx_zero_per_channel = maxima_Rcv_per_channel == 0;
        num_zeros_in_superframe = sum(idx_zero_per_channel);
        if num_zeros_in_superframe==128
            num_empty_superframes = num_empty_superframes+1;
        elseif num_zeros_in_superframe>0
            num_partly_empty_superframes = num_partly_empty_superframes +1;
        end
    end
    
    
    if num_empty_superframes>0 || num_partly_empty_superframes>0
        warning(strcat("There are ", num2str(num_empty_superframes), ...
        " completely empty superframes and ", num2str(num_partly_empty_superframes),...
        " partly empty superframes in the re-ordered RcvData for this acquisition"))
    
    
        warning("Saving of acquired HFR data cancelled since there were (partly) empty frames")
    else
        %% Saving
        Resource.DisplayWindow=[];
        RcvDataShape = size(RcvData);  
        RcvDataType = class(RcvData);
    
        sfNow = datestr(now,30);
        sfnow_date = sfNow(1:8);
        sfnow_time = sfNow(10:end);
        Trans.abb_HM = Trans.name(1:3);
        tube_model = ['T' num2str(P.Tube_nr) '_PM' P.Model_nr];
    
    
        savefilename = [P.savefiledir tube_model '_' sfnow_date '_' sfnow_time '_' Trans.name '_HFR_' num2str(P.ProbeVoltage) 'V_fc_' num2str(P.TransFrequency_true) 'MHz' ];
    
        % save parameters
        save([savefilename,'_info.mat'],'Resource','Trans','TX','TW','P','TGC',...
            'TPC','ReconInfo_DAS','RcvDataType','RcvDataShape','Receive',...
            'displayInfo','-v6'); %,'ImgDataType','ImgDataShape');
    
        fileID = fopen([savefilename,'_RcvData.bin'],'w');
        fwrite(fileID,RcvData,RcvDataType);
        fclose(fileID); 
    
        timer1 = toc;
        fprintf('Done! RcvData saving time was  %.5g \n',timer1)
    end

end
