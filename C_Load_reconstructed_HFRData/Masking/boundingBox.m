function boundbox = boundingBox(coordinates)
%BOUNDINGBOX outputs the minimum bounding box of the points specified in
%coordinates
    xmin = min(coordinates(:,1));
    xmax = max(coordinates(:,1));
    ymin = min(coordinates(:,2));
    ymax = max(coordinates(:,2));
    boundbox = [xmin xmax ymin ymax];
end

