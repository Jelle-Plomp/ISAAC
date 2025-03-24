function masked_image = vertical_mirrorreflectpadding(masked_image, boundary_class, filter_half_size, show_output)
%VERTICAL_MIRRORREFLECT reflects the masked image on the lower and upper
%edge of the mask.
    
% Written by Jelle Plomp, February 2024
    if show_output
        figure(103);clf(103)
        if size(masked_image,4)>1
            subplot(2,1,1); imagesc(masked_image(:,:,:,1)); title("Input masked image (1st frame of stack)")
        else
            subplot(2,1,1); imagesc(masked_image); title("Input masked image")
        end
    end
    % The indices of the top and bottom edge pixels
    idx_top=find(boundary_class==1);
    idx_bottom=find(boundary_class==2);
    
    for i=1:length(idx_top)
        [row_top, col_top] = ind2sub(size(masked_image), idx_top(i));
        for j=1:filter_half_size-1
            masked_image(row_top-j,col_top,:,:) = masked_image(row_top+j, col_top,:,:);
        end
    end
    
    for i=1:length(idx_bottom)
        [row_bottom, col_bottom] = ind2sub(size(masked_image), idx_bottom(i));
        for j=1:filter_half_size-1
            masked_image(row_bottom+j, col_bottom,:,:) = masked_image(row_bottom-j, col_bottom,:,:);
        end
    end
    
    if show_output
        if size(masked_image,4)>1
             subplot(2,1,2); imagesc(masked_image(:,:,:,1)); title("Output masked image after mirroring")
        
        else
            subplot(2,1,2); imagesc(masked_image); title("Output masked image after mirroring")
        end
    end
end

