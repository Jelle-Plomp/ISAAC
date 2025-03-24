%% A general approach to data loading, using uigetdir
% Other option (file
% step2_recon_HFR_L123v_IQ_batch_dataloading_B_finalresults) is
% specifically designed to load files based on information from an excel
% sheet.

% Select folder from which to load the HFR data
dirname = uigetdir('title','Select folder in which HFR data is stored');
files = dir(dirname);
files_name = {files.name};
files_isHFRinfofile = ~cellfun('isempty', regexp(files_name, '_info.mat'));
files_info=files(files_isHFRinfofile);

% Only reconstruct 1 file (so we skip the "batch" processing of the entire
% folder)
select_nr_files = 1; % If 0, then all files are processed
if select_nr_files > 0
     fnames_sel = uigetfile({'*_info.mat'}, 'Select info file of HFR data', dirname, 'MultiSelect','on');
     if isa(fnames_sel, 'char') % in case a single file is selected
         f_selection = find(~cellfun('isempty', regexp({files_info.name}, fnames_sel)));
     else
        for i=1:length(fnames_sel)
            f_selection(i) = find(~cellfun('isempty', regexp({files_info.name}, fnames_sel{i})));
        end
     end
     % Remove all the other structs from the infofiles
     files_info = files_info(f_selection);
end
%%
% Get the name of the experiments folder (such that we can use it to create
% a folder with the same name in the results)
[~,exp_folder_name] = fileparts(extractBefore(dirname,'\Data'));
% Select folder in which to store results
data_dir_results = uigetdir('title','Select folder to store results');
if contains(data_dir_results, exp_folder_name)
    if isequal(dirname, data_dir_results)
        error('Select another folder to store results')
    end
else
    if ~isfolder([data_dir_results, '\', exp_folder_name])
        mkdir([data_dir_results, '\', exp_folder_name])
        data_dir_results = [data_dir_results, '\', exp_folder_name];
    end
end
display(['Results will be stored in', data_dir_results])