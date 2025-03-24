% Script to save IQ data created with step2_recon_HFR... in the format
% required for echoPIV processing.


% Where to store the IQ data
% data_dir_results_IQ = [data_dir_results, '\IQData'];

% if exist('svd_receive', 'var') && svd_receive
%     sname_IQData = [data_dir_results_IQ, '\', fname, '_RcvSVD_min_', num2str(min_reg_val), '_max_', num2str(max_reg_val) '_IQData'];
% else
%     sname_IQData = [data_dir_results_IQ, '\', fname, '_IQData'];
% end
%% Concatenate the frames
IQData_matrix = zeros([size(images_IQ{1}) size(images_IQ,2)]);
for i=1:Nfrs_recon
    IQData_matrix(:,:,i) = images_IQ{i};
end

disp(['File name to be saved is: ' sname_IQData]);
%% Saving
% Save IQData matrix
fileID = fopen([sname_IQData,'.bin'],'w');
DataType = class(IQData_matrix); SaveShape = size(IQData_matrix);
fwrite(fileID,IQData_matrix,DataType);
fclose(fileID);

% Save parameters
Data = 'IQData_DAS';
save([sname_IQData,'_info.mat'],'Data','P', 'ReconInfo_DAS','displayInfo','DataType','SaveShape','Receive','Trans','TX','sname_IQData','-v6');

close all

% Extract some info from the filename
voltage = extractBetween(finfo, "HFR_", "_info"); voltage = voltage{1};
if contains(voltage, '_fc_')
    voltage = extractBefore(voltage, '_fc_');
end
yearnum = year(datetime(files_info(i_file).datenum, 'ConvertFrom','datenum'));
model = extractBefore(fname, ['_', num2str(yearnum)]);

%% Save 1 log compressed image along with the data
% We save the image for easy insight in what the data represents
% x= ReconVars.x;
% z = ReconVars.z;
example_img_hilb = abs(hilbert(IQData_matrix(:,:,1)'));
example_img_log = 20.*log10(example_img_hilb./max(max(example_img_hilb)));

h13 = figure(13); clf(13); 
% set(h13,'position',[50 0 400 800])
colormap("gray")

imagesc(x.*1e3,z.*1e3,example_img_log,[-displayInfo.dynRange 0]); hold on  % Initiate empty frame
ylabel('lateral distance (mm)')
xlabel('axial distance (mm)')
if exist('svd_receive', 'var') && svd_receive
    title_ = strcat(model, " ", voltage, ", PC ", num2str(min_reg_val), " to ", num2str(max_reg_val), " (of ", num2str(size(Sigma,2)), ")");
else
    title_ = strcat(model, " ", voltage);
end
title(title_, 'Interpreter', 'none')
colorbar
caxis([-45 0])
axis equal
axis([x(1) x(end) z(1) z(end)].*1e3)

savefig(h13,[sname_IQData, '.fig'])
saveas(h13,[sname_IQData, '.png'])

disp([fname ' is reconstructed and saved']);

