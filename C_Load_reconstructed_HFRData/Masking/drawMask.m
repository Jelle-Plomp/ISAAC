function [imMask, ROI, poly] = drawMask(im,drange)
    % Function was adapted from script by Jason Voorneveld.

    if (nargin==1)
        drange=[min(im(:)),max(im(:))];
    end
    scnsize = get(0,'ScreenSize');pos=[10 10 scnsize(3)/2 scnsize(4)/2]; %window position on screen
    disp('Draw the borders around where you want to measure. Make sure to close the polygon');
    hfig1=figure('position', pos); hAx = axes; imagesc(im,drange); axis equal
    colormap gray
    %pick rectangle region in source frame for tracking:
    maskBorder = impoly;
    imMask=maskBorder.createMask;
    poly = maskBorder.getPosition;
    temp = boundingBox(poly);
    ROI = [temp(1),temp(3),temp(2)-temp(1),temp(4)-temp(3)]; %Xmin, Ymin, Xwidth, Ywidth
    imrect(hAx,ROI);
    hold on
    rectangle('Position', ROI, 'EdgeColor', 'r'); %show ROI
    hold off
    close(hfig1); 

end