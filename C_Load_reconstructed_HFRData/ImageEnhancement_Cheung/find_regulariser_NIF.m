function [reg_min,NIF, NIF_mean] = find_regulariser_NIF(I_lumen, A_lumen)
%FIND_REGULARISER_NIF finds the regulariser value which minimises the
%normalised intensity fluctuation (NIF) function (equation3 in [1]).

% The NIF is calculated for each frame in I_lumen and all the NIF curves
% are then averaged. The regulariser value at which this mean curve is
% minimum is given as the output.

% Input: 
%   I_lumen: lumen intensity (can be multiple frames)

% Output:
%   reg_min: regulariser value that miniminises the mean NIF
%   NIF: all NIF curves (for all frames)
%   NIF_mean: the mean of the NIF curves.

% [1]   Attenuation Correction and Normalisation for Quantification of
%       Contrast enhancement in Ultrasound Images of Carotid Arteries, 
%       Cheung et al. 2025, Journal of US in Med & Biol.
%       https://doi.org/10.1016/j.ultrasmedbio.2015.02.010 

    
        regularisers = [0.0001:0.0002:0.002, 0.01:0.01:0.2];
        NIF = zeros(size(I_lumen,4), length(regularisers));
        for reg_j = 1:length(regularisers)
            regulariser=regularisers(reg_j);
            CI_lumen = I_lumen./(A_lumen+regulariser);
            
            % Normalised intensity fluctuation per frame
            for frame_i=1:size(CI_lumen,4)
                NIF(frame_i,reg_j) = sqrt(sum((CI_lumen(:,:,:,frame_i)-...
                    mean(CI_lumen(:,:,:,frame_i),"all","omitnan")).^2,"all","omitnan"))...
                    /mean(CI_lumen(:,:,:,frame_i),"all","omitnan");
            end
    
        end
        NIF_mean = mean(NIF,1,'omitnan');
        % make a fit
        regularisers_interp = min(regularisers):0.001:max(regularisers);
        NIF_mean_interp = interp1(regularisers, NIF_mean, regularisers_interp);
        % Find regulariser for which NIF is minimal
        [~,reg_min_idx] = min(NIF_mean_interp);
        reg_min = regularisers_interp(reg_min_idx);
end

