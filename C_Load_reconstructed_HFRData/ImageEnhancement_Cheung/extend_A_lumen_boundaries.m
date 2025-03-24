function [A_full, A_mask]= extend_A_lumen_boundaries(A_lumen,imMask,boundary_class)
%EXTEND_A_LUMEN_BOUNDARIES Extends the lower and upper boundaries of the
% smoothed intensity profile A_lumen.

% Jelle Plomp, September 2024.

% This calculates the values A_ub(x,y) and A_lb(x,y) as used in equation 4
% in [1], except here I chose to first construct the full map A_full
% (including A_lumen, A_ub and A_lb depending on the x,y coordinates),
% after which CI can be obtained at once (outside of this function).


% [1]   Attenuation Correction and Normalisation for Quantification of
%       Contrast enhancement in Ultrasound Images of Carotid Arteries, 
%       Cheung et al. 2025, Journal of US in Med & Biol.
%       https://doi.org/10.1016/j.ultrasmedbio.2015.02.010 
    
    A_full = A_lumen;
    A_mask = double(imMask);
    % Extending the upper boundary (for each frame in A_full)
    idx_top = find(boundary_class==1);
    for i = 1:length(idx_top)
        [row, col] = ind2sub(size(boundary_class), idx_top(i));
        for frame_i=1:size(A_full,4)
            A_full(1:row-1,col,:,frame_i) = A_full(row,col,:,frame_i);
            A_mask(1:row-1,col) = 1;
        end
    end
    % Extending the lower boundary (for each frame in A_full)
    idx_bottom = find(boundary_class==2);
    for i = 1:length(idx_bottom)
        [row, col] = ind2sub(size(boundary_class), idx_bottom(i));
        for frame_i=1:size(A_lumen,4)
            A_full(row+1:end,col,:,frame_i) = A_full(row,col,:,frame_i);
            A_mask(row+1:end,col) = 1;
        end
    end
    % Replace zero values in A_full with nan;
    A_mask(A_mask==0)=nan;
    A_full = A_full.*A_mask;
end

