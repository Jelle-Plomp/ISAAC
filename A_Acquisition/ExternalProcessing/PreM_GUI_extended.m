function PreM_GUI_extended(PreM, ImData, ROI_averageline, ROI_averageline_it1, apod_new, pressure_req,continue_iters, doesExist, History_stats)
%PREM_GUI plots the data during the premeasurement and gives the user the
%information needed to assess whether the active compensation is
%successfull.

% For debugging, run line below:
% PreM_GUI_extended(PreM, rand(100,100), ones(1,length(PreM.ROI.x)),ones(1,length(PreM.ROI.x)) ones(1,128), 1, 0)
% or PreM_GUI_extended(PreM, rand(100,100), ones(1,length(PreM.ROI.x)),ones(1,length(PreM.ROI.x)), ones(1,128),100.*ones(1,length(PreM.rayleigh.x_bounded)), 1, 0)

% Jelle Plomp. 2024.
    figure(101)
    if ~doesExist
        clf(101)
        set(gcf, 'Position', [100 100 1200 700])
        subplot(2,3,1)
       
        % Since imagesc doesn't work well with hold on, we have to set some
        % properties
        title('Latest image in buffer + ROI');
        ax = gca; 
        caxis([-PreM.displayInfo.dynRange 0]); c = colorbar; colormap('gray')
        c.Location = 'south';
        c.AxisLocation = 'out';
        ylabel(c,'Intensity [dB]'); xlabel("Width (mm)"); ylabel("Depth (mm)")
        axis equal
        ax.XLim = [PreM.x(1), PreM.x(end)].*1e3;
        ax.YLim = [PreM.z(1), PreM.z(end)].*1e3;
    
        c.Position = c.Position + [0 -0.2 0 0];

        ax.YDir = 'reverse';
        set(ax,'NextPlot','replacechildren') ; % This command keeps existing axes properties while replacing data.
        
        % Other subplots
        subplot(2,3,2)
        title('Average profile in ROI'); xlabel('x (mm)');ylabel('I(a.u.)');  
        xlim([PreM.x(1) PreM.x(end)].*1e3); grid on; grid minor
        subplot(2,3,5)
        title('Apodization per element'); xlabel('x (mm)'); ylabel('Apodisation'); xlim([PreM.x(1) PreM.x(end)].*1e3); ylim([0 1])
         grid on; grid minor
        % Extended GUI subplots
 
        subplot(2,3,3)
        
        yyaxis left; ylabel("mean(|I_{L,used,i} - I_{L,used,i-1}|/SP)") 
        yyaxis right; ylabel("mean(|Apod_{i} - Apod_{i-1}|)")
        xlabel('Iteration i')
        grid on
        grid minor
        title("Changes in improfile and apod")
        %set(gca, 'FontSize', 14)
        if strcmp(PreM.method,'rayleigh') 
            subplot(2,3,4)
            ylabel('Pi/\gamma')
            xlabel('x')
            title("Requested pressure")
            grid on
            grid minor
        end
        subplot(2,3,6)
        if strcmp(PreM.method,'rayleigh') 
            yyaxis left;  ylabel("I (max and mean)")
            yyaxis right; ylabel("mean(P_{i})")
        end

        xlabel('Iteration i')
        title("Error development")
        grid on
        grid minor
        %set(gca, 'FontSize', 14)
        
    end
    subplot(2,3,1)
    imagesc(PreM.x*1e3, PreM.z*1e3, 20*log10(ImData/max(max(ImData))), [-PreM.displayInfo.dynRange 0]); 
    rectangle('Position', PreM.ROI.Position, 'EdgeColor', [0.47,0.67,0.19], 'LineWidth', 1.5, 'Curvature', [0.1, 0.1]);
    
    % Plot the image profile
    subplot(2,3,2)
    
    hold on; cla;
    plot(PreM.ROI.x.*1e3, ROI_averageline, 'LineWidth', 1.5)
    if strcmp(PreM.method,'rayleigh') 
        plot(PreM.rayleigh.x_bounded.*1e3, PreM.PI_contrl.sp, '--k', 'LineWidth', 1.5)
    end
    if doesExist
        plot(PreM.ROI.x.*1e3, ROI_averageline_it1, '-k', 'LineWidth', 1.5)
    end
    hold off
    
    if continue_iters
        if strcmp(PreM.method,'rayleigh') 
            % plot requested pressure
            subplot(2,3,4)
            hold on
            plot(PreM.rayleigh.x_bounded.*10^3,pressure_req, '-', 'LineWidth', 1.5);
            hold off
        end
        %% Plot apodization
        subplot(2,3,5)
        hold on
        plot(PreM.xel.*10^3,apod_new, '-', 'LineWidth', 1.5);
        hold off
    end
    % Changes in apod, image profile and errors

    it = find(isnan(History_stats.I_change_rel),1,'first')-1;
    if it>1
        subplot(2,3,3); 
        xlim([0 it])
        yyaxis left; cla; hold on;
        plot([2:it],History_stats.I_change_rel(2:it), 'Marker', '.', 'MarkerSize',16, 'LineWidth', 2)
        yyaxis right; cla; hold on;
        plot([2:it],History_stats.Apod_change(2:it), 'Marker', '.', 'MarkerSize',16, 'LineWidth', 2)
    end
    subplot(2,3,6); 
    yyaxis left; cla; hold on
    xlim([0 it])
    if strcmp(PreM.method,'rayleigh') 
       plot(History_stats.I_clip_max, 'Marker', '.', 'MarkerSize',16, 'LineWidth', 2); hold on
       plot(History_stats.I_clip_mean, 'Marker', '.', 'MarkerSize',16, 'LineWidth', 2)

    end

    yyaxis right; cla; hold on
    plot(History_stats.E_cur_mean, 'Marker', '.', 'MarkerSize',16, 'LineWidth', 2)

end

