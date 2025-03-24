function PreM_GUI_default(PreM, ImData, ROI_averageline, ROI_averageline_it1, apod_new, continue_iters, doesExist)
    %PREM_GUI plots the data during the premeasurement and gives the user the
    %information needed to assess whether the active compensation is
    %successfull.
    
    % For debugging, run line below:
    % PreM_GUI_default(PreM, rand(100,100), ones(1,length(PreM.ROI.x)),ones(1,length(PreM.ROI.x)), ones(1,128), 1, 0)
    
    % Jelle Plomp. 2024. 
    figure(101)
    if ~doesExist
        clf(101)
        set(gcf, 'Position', [100 100 800 560])
        subplot(2,2,1); 
        
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
        subplot(2,2,2); 
      
        title('Average profile in ROI'); xlabel('x (mm)');ylabel('I(a.u.)');  
        xlim([PreM.x(1) PreM.x(end)].*1e3); grid on; grid minor
        subplot(2,2,4); 
       
        title('Apodization per element'); xlabel('x (mm)'); ylabel('Apodisation'); xlim([PreM.x(1) PreM.x(end)].*1e3); ylim([0 1])
         grid on; grid minor
    
    end
    subplot(2,2,1);
  
    imagesc(PreM.x*1e3, PreM.z*1e3, 20*log10(ImData/max(max(ImData))), [-PreM.displayInfo.dynRange 0]); 
    rectangle('Position', PreM.ROI.Position, 'EdgeColor', [0.47,0.67,0.19], 'LineWidth', 1.5, 'Curvature', [0.1, 0.1]);
    
    % Plot the image profile
    subplot(2,2,2); 
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
        %% Plot apodization
        subplot(2,2,4); 
        hold on
        plot(PreM.xel.*10^3,apod_new, '-', 'LineWidth', 1.5);
        hold off
    end
end

