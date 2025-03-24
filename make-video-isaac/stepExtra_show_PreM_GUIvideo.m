%% Loading Data, creating videos (in batch)
clear all; close all;
% addpath('A_Acquisition\ExternalProcessing\')
% addpath([extractBefore(pwd, 'L12-3 attenuation correction') 'Visualisation_General\Colormaps'])
git_repo_topdir = 'C:\Users\PlompJ\OneDrive - University of Twente\Matlab_git_br_pub'; % Replace this

dirname = [git_repo_topdir '\Data_examples\iterative_procedure'];
files=dir(dirname); 
files_name = {files.name};
files_isPreMinfofile = ~cellfun('isempty', regexp(files_name, '_infoPreM.mat'));
files_PreMinfo=files(files_isPreMinfofile);
%%
selection = 1; % for just 1 file, just input 1 number.

for file_z = selection
    close all;
    fname = files_PreMinfo(file_z).name;
    f = find(fname=='_');
    fname = fname(1:f(end)-1); % last '_'
    
    fdata = [fname,'_PreM_ImageData.bin'];
    finfo = [fname,'_infoPreM.mat'];
    
  
    load(fullfile(dirname, finfo),'PreM','History','ProcTimes','History_stats','x','z', 'Resource','Trans','TX','TW','P','TGC',...
        'TPC','ReconInfo_DAS','ImageDataType','ImageDataShape','Receive');
    fid = fopen(fullfile(dirname,fdata));
    ImageData = fread(fid,prod(ImageDataShape),['*' ImageDataType]);
    ImageData = reshape(ImageData,ImageDataShape);
    fclose(fid);
    %% Get the tube thickness
    if str2num(P.Tube_nr) == 11
        d_ridge_array = [0.25, 0.50, 1.5];
    end
    d_ridge = d_ridge_array(str2num(P.Model_nr));
    %% Loading the mask
    load(fullfile(dirname, P.filename_ROImask))
    
    PreM.GUI = 'extended';
    PreM.displayInfo.dynRange = 55;
    calcHistory_stats = ~exist('History_stats', 'var');

    % =====================================================================
    % Pre-calculations just to calculate maximum and minimum errors and
    % apod changes
    % ================================
    if calcHistory_stats
        for iter_i = 1:length(History)
        
        
            if iter_i>1
                History_stats.Apod_change(iter_i) = mean(abs((History(iter_i).Apod-History(iter_i-1).Apod)));
            else
               History_stats.Apod_change = nan(1, length(History));
                History_stats.Apod_change(1) = 0;
                History_stats.E_tot_max = nan(1, length(History));
%               History_stats.E_tot_mean_abs = zeros(1, length(History));
                History_stats.E_tot_mean = nan(1, length(History));
                History_stats.E_cur_mean = nan(1, length(History));
                
            end
            History_stats.E_tot_max(iter_i) = max(abs(History(iter_i).error_total));
            History_stats.E_tot_mean(iter_i) = mean(History(iter_i).error_total);
            error_cur = (PreM.PI_contrl.sp - History(iter_i).ROI_128)./PreM.PI_contrl.sp;
            History_stats.E_cur_mean(iter_i) = mean(error_cur);
        end
    end
    History_stats.maxiter = length(History);
    if strcmp(PreM.method,'rayleigh')
        for j=1:length(History)-1
            History_stats.I_clip_max(j) = max(History(j).I_clip);
        end
        History_stats.I_clip_max_min = min(History_stats.I_clip_max);
        History_stats.I_clip_max_max = max(History_stats.I_clip_max);
    end
    History_stats.Apod_change_max = max(History_stats.Apod_change);
    History_stats.Apod_change_min = min(History_stats.Apod_change);
    History_stats.E_cur_mean_max= max(History_stats.E_cur_mean );
    History_stats.E_cur_mean_min= min(History_stats.E_cur_mean );
    History_stats.I_ROI_line_min = min([History(:).I_ROI_line]);
     History_stats.I_ROI_line_max = max([History(:).I_ROI_line]);
    
    % ====================================================================
    % Video
    % ====================================================================
%         History_stats.E_tot_mean_abs(it) = mean(abs(History(end).error_total));
    sname_video = [fullfile(dirname, [fname, '_PreMvideo.avi'])];
    v = VideoWriter(sname_video);
    v.FrameRate = 1; % just 1 frame per second
    
    open(v);
    History_stats_fields = fieldnames(History_stats);
     for iter_i = 1:length(History)
            % Use values from oriignal history_stats, up until current
            % iteration
            for field_idx = 1:length(History_stats_fields)
                if length(History_stats.(char(History_stats_fields{field_idx})))>1
                    field_preall = nan(1,length(History));
                    field_history_statsorig = History_stats.(char(History_stats_fields{field_idx}));
                    field_preall(1:iter_i) = field_history_statsorig(1:iter_i);
                    History_stats_iter_i.(char(History_stats_fields{field_idx})) = field_preall;
                else % for single values, value is the same for each iter
                    History_stats_iter_i.(char(History_stats_fields{field_idx})) = History_stats.(char(History_stats_fields{field_idx}));
                end
            end
        if iter_i==length(History)
            continue_iters = false;
            imagedata = nan;
        else
            continue_iters = true;
            imagedata = ImageData(:,:,:,iter_i);
            % subplot(2,3,2);cla
        end
        
        if strcmp(PreM.method,'rayleigh')
            if iter_i==1
                pressure_req = PreM.rayleigh.P0;
            else
                pressure_req = History(iter_i-1).P + History(iter_i-1).I_clip + PreM.rayleigh.P0;
            end
            PreM_GUI_presentation_v_rayleigh(PreM, imagedata, ...
            History(iter_i).I_ROI_line,History(1).I_ROI_line, ...
            History(iter_i).Apod, continue_iters, (iter_i~=1), History_stats_iter_i,pressure_req);
            

        end
        
        if exist('h1', 'var')
            delete(h1)
        end
%         h1=annotation('textbox',[0.4, 0.26, 0.2, 0.05], ...
%             'String',['Iteration  ' num2str(iter_i)],'EdgeColor','none', 'FontSize', 14);
        if iter_i == length(History)
            sgtitle('Final apodisation $A_{ISAAC}$','Interpreter','latex')
            annotation('textbox',[0.55, 0.75, 0.2, 0.05], ...
            'String','$A_{ISAAC}$','EdgeColor','none', 'FontSize', 14, 'Color', [252, 3, 252]./255,'Interpreter','latex');
        else
            sgtitle(['Iteration  $i$ = ' num2str(iter_i-1)],'Interpreter','latex')
        end
            
        if iter_i==1
            % Text description
            
            y_pos = 0.5;
            x_pos = 0.02;
            if str2num(P.Tube_nr) == 5
                annotation('textbox',[0.4, y_pos, 0.2, 0.05], ...
                'String',['t_{ridge} = ' num2str(t_ridge(str2num(P.Model_nr))) ' mm'] ,'EdgeColor','none', 'FontSize', 14);
            elseif str2num(P.Tube_nr) == 6
                annotation('textbox',[0.4, y_pos, 0.2, 0.05], ...
                'String',['T6, turbulent flow'] ,'EdgeColor','none', 'FontSize', 14);
            end
            % y_pos = y_pos-0.04;
            annotation('textbox',[x_pos, y_pos-0.03, 0.2, 0.05], ...
            'String',['$f_c$ = ' num2str(P.TransFrequency_true) ' MHz'] ,'EdgeColor','none', 'FontSize', 14,'Interpreter','latex');
            % y_pos = y_pos-0.04;
            annotation('textbox',[x_pos, y_pos-0.07, 0.2, 0.05], ...
            'String',[ num2str(History(1).Voltage) ' V' ] ,'EdgeColor','none', 'FontSize', 14,'Interpreter','latex');
            % y_pos = y_pos-0.04;
            % annotation('textbox',[0.75, y_pos, 0.3, 0.05], ...
            % 'String',['Smoothing with \sigma = ' num2str(PreM.std_mm) ' mm'],'EdgeColor','none', 'FontSize', 14);
            % y_pos = y_pos-0.04;
            annotation('textbox',[x_pos, y_pos+0.08, 0.2, 0.05], ...
            'String',['$K_P$ = ' num2str(PreM.PI_contrl.Kp)],'EdgeColor','none', 'FontSize', 14, 'Color', 'k','Interpreter','latex');
            annotation('textbox',[x_pos, y_pos+0.04, 0.2, 0.05], ...
            'String',['$K_I$ = ' num2str(PreM.PI_contrl.Ki)],'EdgeColor','none', 'FontSize', 14, 'Color', 'k','Interpreter','latex');
             
%             annotation('textbox',[0.4, y_pos, 0.3, 0.05], ...
%             'String',['SP version:  ' PreM.PI_contrl.sp_version],'EdgeColor','none', 'FontSize', 14, 'Color', 'k', 'Interpreter','none');
            annotation('textbox',[x_pos, 0.92, 0.2, 0.05], ...
            'String',['$d_{ridge}$ = ' sprintf('%.2f',d_ridge) ' mm'],'EdgeColor','none', 'FontSize', 14, 'Color', 'k','Interpreter','latex');

            % P0 and A0
            text_p0=annotation('textbox',[0.23, y_pos+0.05, 0.2, 0.05], ...
            'String','$P_0$','EdgeColor','none', 'FontSize', 14, 'Color', 'k','Interpreter','latex');
            text_a0=annotation('textbox',[0.5, y_pos+0.07, 0.2, 0.05], ...
            'String','$A_0$','EdgeColor','none', 'FontSize', 14, 'Color', 'k','Interpreter','latex');
        else
            delete(text_p0);delete(text_a0)
        end
    
        drawnow
        videoframe = getframe(gcf);
        writeVideo(v,videoframe)
    end
    
        
    close(v)
   clearvars History_stats_iter_i
end

function History_stats_empty=preallocate_history_stats(length_history)
      History_stats_empty.I_change_rel = nan(1, length_history);
    History_stats_empty.I_change_rel(1) = 0;
    History_stats_empty.Apod_change = nan(1, length_history);
    History_stats_empty.Apod_change(1) = 0;
    History_stats_empty.E_tot_max = nan(1, length_history);
    History_stats_empty.E_tot_mean = nan(1, length_history);
    History_stats_empty.E_cur_mean = nan(1, length_history);
end