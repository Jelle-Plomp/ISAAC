function save_mean_image(mean_frame, x, z, displayInfo, finfo, files_info_i_struct, savename_frame,P,ReconInfo_DAS,Receive,Trans,TX)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%

%
% Save parameters
Data = 'IQData_DAS';
save([savename_frame,'.mat'],'mean_frame','Data','P', 'ReconInfo_DAS','displayInfo','Receive','Trans','TX','savename_frame','-v6');

% extract some info
voltage = extractBetween(finfo, "HFR_", "_info"); voltage = voltage{1};
if contains(voltage, '_fc_')
    voltage = extractBefore(voltage, '_fc_');
end
yearnum = year(datetime(files_info_i_struct.datenum, 'ConvertFrom','datenum'));
model = extractBefore(finfo, ['_', num2str(yearnum)]);

%% Save an image of a frame (or mean frame)
% We save the image for easy insight in what the data represents
% x= ReconVars.x;
% z = ReconVars.z;

example_img_log = 20.*log10(mean_frame'./max(max(mean_frame)));

h14 = figure(14); clf(14); 
colormap("gray")

imagesc(x.*1e3,z.*1e3,example_img_log); hold on  % Initiate empty frame
ylabel('lateral distance (mm)')
xlabel('axial distance (mm)')

title_ = strcat(model, " ", voltage);

title(title_, 'Interpreter', 'none')
colorbar
clim([-displayInfo.dynRange 0])
axis equal
axis([x(1) x(end) z(1) z(end)].*1e3)


saveas(h14,[savename_frame, '.png'])

end