function PreM_GUI_presentation_v_rayleigh(PreM, ImData, ROI_averageline, ROI_averageline_it1, apod_new, continue_iters, doesExist, History_stats, pressure_req)
%PREM_GUI_PRESENTATION_V_RAYLEIGH plots the saved data from during ISAAC.

%The progress of various variables during the iterative procedure is 
% plotted. 

% Jelle Plomp.


    if ~strcmp(PreM.method,'rayleigh')
        error("only use this function when PreM.method == 'rayleigh'")
    end
    iter_i = find(isnan(History_stats.I_change_rel),1,'first')-1;
    if isempty(iter_i)
        iter_i=length(History_stats.I_change_rel);
    end
    colours_iter = jet(History_stats.maxiter);
    colours_iter(1,:) = [0 0 0];

    % Define which plot goes where in the GUI
    plot_idx_improfile = 5;
    plot_idx_image = 4;
    plot_idx_apod = 2;
    plot_idx_apodchanges = 3;
    plot_idx_error_development = 6;
    plot_idx_pressure = 1;
    figure(101); 
    if ~doesExist
        clf(101)

        set(gcf,'color','w');
        set(gcf, 'Position', [100 100 1200 700])
        % ====================== Image subplot ===========================
        subplot(2,3,plot_idx_image)
       
        % Since imagesc doesn't work well with hold on, we have to set some
        % properties
        title('Acquired image + ROI','Interpreter','latex');
        ax = gca; 
        clim([-PreM.displayInfo.dynRange 0]); colormap('gray');
        c = colorbar(ax,"west"); 
        % c.Location = 'west';
        c.AxisLocation = 'out';
        ylabel(c,'Intensity [dB]','Interpreter','latex'); xlabel("Width (mm)",'Interpreter','latex'); ylabel("Depth (mm)",'Interpreter','latex')
        axis equal
        ax.XLim = [PreM.x(1), PreM.x(end)].*1e3;
        ax.YLim = [PreM.z(1), PreM.z(end)].*1e3;
    
        c.Position = c.Position + [-0.08 -0.01 0 0];
        c.Label.Position(1) = 1;
        c.Label.Color = [1 1 1];
        ax.YDir = 'reverse';
        set(ax,'NextPlot','replacechildren') ; % This command keeps existing axes properties while replacing data.

        % ====================== Improfile subplot ========================
              subplot(2,3,plot_idx_improfile)
      
        title('Average profile in ROI','Interpreter','latex'); xlabel('$x_q$ (mm)','Interpreter','latex');ylabel('$I_L$ (a.u.)','Interpreter','latex');  
        xlim([PreM.x(1) PreM.x(end)].*1e3); grid on; grid minor
  
        % ====================== Apodisation subplot ========================
        subplot(2,3,plot_idx_apod)
    
        title('Apodisation $A_i$','Interpreter','latex'); xlabel('$x_e$ (mm)','Interpreter','latex'); ylabel('Apodisation $a_e$','Interpreter','latex'); xlim([PreM.x(1) PreM.x(end)].*1e3); ylim([0 1])
        grid on; grid minor

        % ====================== Apodisation changes subplot ========================
      
        subplot(2,3,plot_idx_apodchanges)
        ylabel("mean($|A_{i} - A_{i-1}|$)",'Interpreter','latex')

        xlabel('Iteration $i$','Interpreter','latex')
        grid on
        grid minor
        title("Changes in apodisation",'Interpreter','latex')
        ylim([History_stats.Apod_change_min History_stats.Apod_change_max])
        
        % ====================== Error development subplot ========================
        subplot(2,3,plot_idx_error_development)
        yyaxis right; ylabel("mean($E_{P}$)",'Interpreter','latex')
        yyaxis left;ylabel("max($E_{I,clip}$)",'Interpreter','latex')
   
        xlabel('Iteration $i$','Interpreter','latex')
        title("Error development",'Interpreter','latex')
        grid on
        grid minor
        %set(gca, 'FontSize', 14)
        % ====================== Requested pressure subplot ========================
        subplot(2,3,plot_idx_pressure)
        ylabel("$P_{req}$ (a.u.)", 'Interpreter','latex')
        
        xlabel('$x_q$ (mm)','Interpreter','latex')
        title("Required pressure $P_{req}$",'Interpreter','latex')
        grid on
        grid minor
        xlim([PreM.x(1) PreM.x(end)].*1e3); 
        ylim([0 max(PreM.rayleigh.P1)])
    end

    %% Plot image
    subplot(2,3,plot_idx_image)
   
    if continue_iters
        imagesc(PreM.x*1e3, PreM.z*1e3, 20*log10(ImData/max(max(ImData))), [-PreM.displayInfo.dynRange 0]); 
        rectangle('Position', PreM.ROI.Position, 'EdgeColor', [0.47,0.67,0.19], 'LineWidth', 1.5, 'Curvature', [0.1, 0.1]); 
    end
    %% Plot the image profile
    subplot(2,3,plot_idx_improfile)
    if continue_iters
        cla; hold on;
        plot(PreM.ROI.x.*1e3, ROI_averageline, 'LineWidth', 1.5, 'Color', colours_iter(iter_i,:))

        plot(PreM.ROI.x.*1e3, PreM.PI_contrl.sp, '--k', 'LineWidth', 1.5, 'Color', colours_iter(1,:))
        
        plot(PreM.ROI.x.*1e3, ROI_averageline_it1, '-k', 'LineWidth', 1.5, 'Color', colours_iter(1,:))
            
        hold off
        ylim([0 History_stats.I_ROI_line_max])
    end
    
    %% Plot apodization
    subplot(2,3,plot_idx_apod)
    
    hold on
    if continue_iters
        plot(PreM.xel.*10^3,apod_new, '-', 'LineWidth', 1.5,  'Color', colours_iter(iter_i,:));
    else
        plot(PreM.xel.*10^3,apod_new, '-', 'LineWidth', 1.5,  'Color',[252, 3, 252]./255);
    end 
    hold off
    
    %% Changes in apod
    if iter_i>1
        subplot(2,3,plot_idx_apodchanges); 
        xlim([0 History_stats.maxiter])
        cla; hold on;
%         plot([2:iter_i],History_stats.Apod_change(2:iter_i), 'Marker', '.', 'MarkerSize',16, 'LineWidth', 2)
        plot([1:iter_i-1],History_stats.Apod_change(2:iter_i),'LineWidth', 1)
        scatter([1:iter_i-1],History_stats.Apod_change(2:iter_i), 18, colours_iter(2:iter_i,:), 'filled')
        ylim([History_stats.Apod_change_min History_stats.Apod_change_max])
    end
    %% Error development
    subplot(2,3,plot_idx_error_development); 
    cla; hold on
    xlim([0 History_stats.maxiter])
    yyaxis left; cla
%     plot(History_stats.E_tot_max, 'Marker', '.', 'MarkerSize',16, 'LineWidth', 2); hold on
 
    plot([0:iter_i-1],History_stats.I_clip_max(1:iter_i), 'LineWidth',1)
    scatter([0:iter_i-1],History_stats.I_clip_max(1:iter_i), 18, colours_iter(1:iter_i,:), 'filled')% 'Marker', '.', 'MarkerSize',16, 'LineWidth', 2); hold on
    if PreM.PI_contrl.Ki>0
    ylim([History_stats.I_clip_max_min History_stats.I_clip_max_max])
    end
       % plot(History_stats.I_clip_mean, 'Marker', '.', 'MarkerSize',16, 'LineWidth', 2)
    
    yyaxis right; cla; hold on
%     plot(History_stats.E_cur_mean, 'Marker', '.', 'MarkerSize',16, 'LineWidth', 2)
    plot([0:iter_i-1],History_stats.E_cur_mean(1:iter_i), 'LineWidth', 1)
    scatter([0:iter_i-1],History_stats.E_cur_mean(1:iter_i), 18, colours_iter(1:iter_i,:), 'filled')

    xlim([0 History_stats.maxiter])
    
    ylim([History_stats.E_cur_mean_min History_stats.E_cur_mean_max])
    for subplot_i=[1,2,3,4,5,6]
        subplot(2,3,subplot_i); set(gca, 'FontSize', 14)
    end
    %% Pressure
    subplot(2,3,plot_idx_pressure); hold on
    if continue_iters
        plot(PreM.ROI.x.*1000, pressure_req, '-', 'LineWidth', 1.5,  'Color', colours_iter(iter_i,:))
    else
         plot(PreM.ROI.x.*1000, pressure_req, '-', 'LineWidth', 1.5,  'Color', [252, 3, 252]./255)
    end
    
end

