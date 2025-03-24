function SavePreM_L12_3v(ImageData)
%SAVEPREM_L12_3V Saves the ImageData and the history struct. By inputting the same data
% here as is used in the premeasurement (PreM), it is ensured that after the
% measurement one can look back at the data based on which the apodization,
% voltage and local gain were changed.
tic


% Parameters to be saved with the data
Resource = evalin('base','Resource');
P = evalin('base','P'); 
ProcTimes = evalin('base','ProcTimes'); 
TPC = evalin('base','TPC');
Trans = evalin('base','Trans'); 
PData = evalin('base','PData');
Receive = evalin('base','Receive(end)'); % only saving last Receive (instead of thousands) makes it much faster
TX = evalin('base','TX');
TW = evalin('base','TW');
ReconInfo_DAS = evalin('base','ReconInfo_DAS');
TGC = evalin('base','TGC');
ReconVars = evalin('base','ReconVars');
incl_Hist_stats = evalin('base', "exist('History_stats', 'var')");
if incl_Hist_stats
    History_stats = evalin('base', 'History_stats');
end

x = ReconVars.x;
z = ReconVars.z;
P.ProbeVoltage = TPC(1).hv; % TPC 2 = HFR saving (should only be the case when multiple TPC profiles are created, which we do not do)
% History struct
History = evalin('base','History');
PreM = evalin('base', 'PreM');
% Replace NaN values in History(end) by the final voltage and apodisation.
if isnan(History(end).Voltage)
    History(end).Voltage = P.ProbeVoltage;
    History(end).Apod = TX(1).Apod;
end
% Cut off all the empty frames
ImageData = ImageData(:,:,:,1:length(History)-1); % Last History element is of a not-yet acquired frame, hence the -1 here.
ImageDataShape = size(ImageData);  
ImageDataType = class(ImageData);

% Savefilename
sfNow = datestr(now,30);
sfnow_date = sfNow(1:8);
sfnow_time = sfNow(10:end);
Trans.abb_HM = Trans.name(1:3);
tube_model = ['T' num2str(P.Tube_nr) '_PM' P.Model_nr];


savefilename = [P.savefiledir tube_model '_' sfnow_date '_' sfnow_time '_' Trans.name '_' num2str(P.ProbeVoltage) 'V_fc_' num2str(P.TransFrequency_true) 'MHz' ];
% save parameters
fileID = fopen([savefilename,'_PreM_ImageData.bin'],'w');
fwrite(fileID,ImageData,ImageDataType);
fclose(fileID); 
% Save the name of the file, such that this info can be used for the HFR
% saving as well
if incl_Hist_stats
    save([savefilename,'_infoPreM.mat'],'History','PreM','History_stats','ProcTimes', 'x', 'z','Resource','Trans','TX','TW','P','TGC',...
    'TPC','ReconInfo_DAS','ImageDataType','ImageDataShape','Receive', '-v6'); %,'ImgDataType','ImgDataShape');

else
    save([savefilename,'_infoPreM.mat'],'History','PreM','ProcTimes', 'x', 'z','Resource','Trans','TX','TW','P','TGC',...
    'TPC','ReconInfo_DAS','ImageDataType','ImageDataShape','Receive', '-v6'); %,'ImgDataType','ImgDataShape');
end
P.filename_ImageDataPreM = savefilename;
assignin('base','P',P)
timer1 = toc;
fprintf('Done! Premeasurement saving time was  %.5g \n',timer1)

end

