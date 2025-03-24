function displayLogComp(imageData)
    %DISPLAYLOGCOMP Displays the input imageData after log-compressing it.
    %The input data is assumed to be hilbert transformed data of the DAS 
    % reconstruction.
    
    % Jelle Plomp. 2024.

    tic_disp = tic;
    x = evalin('base', 'ReconVars.x');
    z = evalin('base', 'ReconVars.z');

    hilb = imageData;
    
    dynRange = evalin('base', 'displayInfo.dynRange');
    persistent frame
    if isempty(frame)
        frame = 0;
    else 
        frame = frame+1;
    end

    f11=figure(11); 
    screensize= get(0,'screensize');
    f11.Position=[0.1*screensize(3) 0.5*screensize(4) 0.5*screensize(3) 0.4*screensize(4)];
    subplot(1,2,1);ax=gca;
    isEmptyAxes = isempty(get(ax, 'Children'));
    if isEmptyAxes % Only for the first frame
        title("Live View")
        % Since imagesc doesn't work well with hold on, we have to set some
        % properties
        set(ax,'NextPlot','replacechildren') ;
        ax.XLim = [x(1), x(end)].*1e3;
        ax.YLim = [z(1), z(end)].*1e3;
        axis equal
        ax.YDir = 'reverse';
        caxis([-dynRange 0])
        c = colorbar; 
        colormap('gray')
        ylabel(c,'Intensity [dB]');
        xlabel("x (mm)")
        ylabel("z (mm)")
    end

    cla
    logcomp = 20*log10(hilb/max(hilb(:)));
    imagesc(x.*1e3, z.*1e3, logcomp, [-dynRange 0]);
    hold on
    plot([x(1) x(end)].*1e3, [26 26], 'LineWidth', 0.5)
    plot([x(1) x(end)].*1e3, [32.8 32.8], 'LineWidth', 0.5)
    
    % In a second axis, the mean image intensity in a predefined region of
    % interest (see line 67) is plotted, along with a reference line at
    % 0.3.
    subplot(1,2,2); ax2=gca;
    isEmptyAxes = isempty(get(ax2, 'Children'));
    if isEmptyAxes % Only for the first frame
        title("Image intensity")
        % Since imagesc doesn't work well with hold on, we have to set some
        % properties
        set(ax2,'NextPlot','replacechildren') ;
        ax2.XLim = [x(1), x(end)].*1e3;
        ax2.YLim = [0 1];

        xlabel("x (mm)")
    end
    
    cla
    [~,idx_z29mm] = min(abs(z-0.029)); hold on
    plot(x.*1000,mean(hilb(idx_z29mm-10:idx_z29mm+10,:),1))
    plot(x.*1000,0.3.*ones(size(x)),'-r')
%     title(num2str(frame))
%     title(num2str(Resource.ImageBuffer(1).lastFrame))
    ProcTimes = evalin('base', 'ProcTimes');
    ProcTimes(end).displayLogComp = toc(tic_disp);
    assignin('base', 'ProcTimes', ProcTimes);
end
