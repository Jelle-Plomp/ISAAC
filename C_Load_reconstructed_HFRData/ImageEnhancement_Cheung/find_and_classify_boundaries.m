function [boundaries, boundary_class] = find_and_classify_boundaries(mask_binary, show_output)
    %FIND_AND_CLASSIFY_BOUNDARIES Finds boundaries using bwperim, and then
    %classifies them as top (1), bottom (2), left (3), right (4) or corner (5) segment.
    
    % Written by Jelle Plomp, February 2024
    boundaries=bwperim(mask_binary);

    if show_output
        figure(101);clf(101)
        subplot(1,2,1)
        imagesc(boundaries)
        axis equal
        title("Boundary pixels")
    end
    % indices of boundary elements
    boundary_idx = find(boundaries);
    boundary_class = zeros(size(boundaries));
    % Classification of the pixels in "boundaries"
    for i=1:length(boundary_idx)
    
        [row,col] = ind2sub(size(boundaries),boundary_idx(i));
        if row>1 && col>1 && row<size(boundaries,1) && col<size(boundaries,2)
            if boundaries(row-1,col) == 1 || boundaries(row+1,col) == 1 
                if boundaries(row,col-1) == 1 || boundaries(row,col+1) == 1 % Corner segment
                    % Don't mirror? Or save indices until the end (after other
                    % parts have been mirrored)
                    boundary_class(row,col)=5;
                else % Vertical segment
                    if mask_binary(row,col+1) %Left edge
                       boundary_class(row,col)=3; 
                    elseif mask_binary(row,col-1) % Right edge
                        boundary_class(row,col)=4; 
                    end
                end
            else % Horizontal edge segment 
                if mask_binary(row+1,col) %Top edge
                       boundary_class(row,col)=1; 
                elseif mask_binary(row-1,col) % Bottom edge
                        boundary_class(row,col)=2; 
                end
            end
        end
        % Corners of the image, can only be corners of the perimeter
        if (row==1 || row==size(boundaries,1)) && ((col==1 || col==size(boundaries,2)))
            boundary_class(row,col)=5;
        elseif row==1
            boundary_class(row,col)=1;
        elseif row==size(boundaries,1)
            boundary_class(row,col)=2;
        elseif col==1
            boundary_class(row,col)=3; 
        elseif col==size(boundaries,2)
             boundary_class(row,col)=4;
        end

    end

    if show_output
        subplot(1,2,2)
        imagesc(boundary_class)
        title("Boundary classification")
        axis equal
    end
end

