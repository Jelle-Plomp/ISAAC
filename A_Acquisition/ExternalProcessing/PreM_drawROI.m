function [ROI, poly] = PreM_drawROI(x,z,image, dynRange, showfig)
% DRAWROI can be used to draw an ROI over an image and convert it to a
% rectangle (for use in the attenuation correction)
% Jelle Plomp, 13-1-2023
% Code to draw and save an ROI

% Starting point for developing this function: "select_ROIs_calcum_phantoms.m"
% by Majorie van Helvert.
%{
For the purpose of the attenuation correction, we need: 
    A polygon ROI drawn over a part of the tube that includes the shadow in
    the middle, and the two non-shadowed regions on the left and right.
        - This is the output poly (with a mask and coordinates)
        - The poly is not used for the compensation strategy. In a future
        version, the compensation could be done over the polygon region
        instead of the rectangle, but the algorithm would require some
        adjustments (currently the mean image intensity is obtained by
        averaging the pixel rows in the rectangle)
    
    A rectangle ROI that is fully within the vessel.
 
        - Based on this rectangle, the apodization can be adjusted. We can 
          average the intensity profile at different pixel rows (different 
          depths z) to get a smoother intensity profile based on which the
          apodization can then be adjusted. This assumes the shadow egdes 
          are vertical. 
        - This is the output ROI, with ROI.Position the rectangle Position
            in image coordinates, and ROI.Position_pix the same in pixel
            coordinates (this is convenient for data selection).
        
     
%}
    set(0, 'DefaultAxesFontSize', 12, 'defaultAxesFontname', 'Calibri Light')



    scnsize = get(0,'ScreenSize');
    pos     = [10 10 scnsize(3)/2 scnsize(4)/2]; %window position on screen
    hfig1   = figure('position', pos); 

    disp(['Please draw 4 points, starting top left and work counterclockwise']);

    h=imagesc(x*1e3, z*1e3, image, [-dynRange 0]); 
    xlabel('Width [mm]'); ylabel('Depth [mm]'); 
    title(['Please, draw ROI' ]);

    colormap gray
    c = colorbar; 
    ylabel(c,'Intensity [dB]'); 
    
    % Get the polygon mask and coordinates
    maskBorder      = impoly;
    poly.Mask          = maskBorder.createMask; % think imMask is enough
    poly.Coor           = maskBorder.getPosition;
    % Calculate "inner bounding box" (mask + coordinates)
    x_left = max([poly.Coor(1,1) poly.Coor(2,1)]); % min x of left 2 points
    y_top = max([poly.Coor(1,2) poly.Coor(4,2)]); % max y of top 2 points
    x_right = min([poly.Coor(3,1) poly.Coor(4,1)]); % min of right 2 points
    y_bottom = min([poly.Coor(2,2) poly.Coor(3,2)]); % min of bottom 2 points
    
    ROI.Position = [x_left, y_top, x_right-x_left, y_bottom-y_top];% Xmin, Ymin, Xwidth, Ywidth
    
    % Non-rounded for ROImask
    ROIBox_object= images.roi.Rectangle(gca,'Position', ROI.Position);
    ROI.Mask = createMask(ROIBox_object);
    % Convert axis coordinates to pixel coordinates
    x_left_pix = round(axes2pix(size(image,2),get(h,'XData'),x_left));
    x_right_pix = round(axes2pix(size(image,2),get(h,'XData'),x_right));
    y_top_pix = round(axes2pix(size(image,1),get(h,'YData'),y_top));
    y_bottom_pix = round(axes2pix(size(image,1),get(h,'YData'),y_bottom));
    
    ROI.Position_pix = [x_left_pix, y_top_pix, x_right_pix-x_left_pix, ...
                        y_bottom_pix-y_top_pix];
    close(hfig1)

    ROI.x = x(ROI.Position_pix(1):ROI.Position_pix(1)+ROI.Position_pix(3));
    ROI.z =  z(ROI.Position_pix(2):ROI.Position_pix(2)+ROI.Position_pix(4)); % Here we use z as it represents the Transducer coordinate frame

    % Save ROI and poly objects
    P = evalin('base', 'P');

    sfNow = datestr(now,30);
    sfnow_date = sfNow(1:8);
    sfnow_time = sfNow(10:end);
    tube_model = ['T' num2str(P.Tube_nr) '_PM' P.Model_nr];

    P.filename_ROImask = [tube_model '_' sfnow_date '_' sfnow_time '_Mask' ];
    
    save([P.savefiledir P.filename_ROImask], 'ROI', 'poly')
    assignin('base','P', P) % We store the savefilename of the mask.

    % show ROI
    if showfig
        figure()
        imagesc(x*1e3, z*1e3, image, [-dynRange 0]); 
        colormap gray
        c = colorbar; 
        ylabel(c,'Intensity [dB]'); 
        hold on 
        rectangle('Position', ROI.Position, 'EdgeColor', [0.47,0.67,0.19], 'LineWidth', 1.5, 'Curvature', [0.1, 0.1]);
        ps = polyshape(poly.Coor);
%         plot(ps, 'FaceAlpha', 0.5, 'FaceColor', [0.47,0.67,0.19], 'EdgeAlpha', 0);
        xlabel('Width [mm]'); ylabel('Depth [mm]'); 
        title('ROI box for attenuation correction');
        
        % Save figure of the ROI
        TPC = evalin('base', 'TPC');
        voltage = num2str(TPC(1).hv);
        if ~isfile([P.savefiledir P.filename_ROImask '_' voltage 'V_fig.png'])
            saveas(gca, [P.savefiledir P.filename_ROImask '_' voltage 'V_fig.png']);
            saveas(gca, [P.savefiledir P.filename_ROImask '_' voltage 'V_fig.fig']);
        end 
        
    end 
end