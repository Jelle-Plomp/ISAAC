function [image_resized, pixel_shortage_top, pixel_shortage_bottom] = resize_image(image, boundary_class,std_datapoints, show_output)
%RESIZE_IMAGE Summary of this function goes here
%   Works for single image as well as stack of images (stack dims: row,
%   col,[],frames)
    
% Written by Jelle Plomp, February 2024
    % The indices of the top and bottom edge pixels
    idx_top=find(boundary_class==1);
    idx_bottom=find(boundary_class==2);
    
    % Calculate how much the image needs to be enlarged to be able to apply
    % the mirroring
    % First, find highest and lowest point in mask
    [row_bottom, ~] = ind2sub(size(image), idx_bottom);
    row_max = max(row_bottom);
    [row_top, ~] = ind2sub(size(image), idx_top);
    row_min = min(row_top);
    %% Extend the image bounds
    filter_half_size = ceil(2*std_datapoints);
    pixel_shortage_top = max([0 filter_half_size-row_min]);
    pixel_shortage_bottom = max([0 filter_half_size-(size(image,1)-row_max)]);
    
    if pixel_shortage_bottom==0 && pixel_shortage_top ==0
        image_resized = image;
    else
        image_resized = zeros(size(image,1)+pixel_shortage_bottom+pixel_shortage_top, size(image,2), size(image,3), size(image,4));
        image_resized(pixel_shortage_top+1:end-pixel_shortage_bottom,:,1,:)=image;
    end
    if show_output
        figure(102);clf(102)
        if size(image,4)>1
            imagesc(image_resized(:,:,:,1))
            title("resized image (1st image in stack)")
        else
            imagesc(image_resized)
        end
        title("resized image")
    end

    
end

