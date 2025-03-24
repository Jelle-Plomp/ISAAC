% This script performs an image-based attenuation correction and normalisation as described
% by Cheung et al [1]. Some results are presented in the supplementary
% material of our paper Improved blood flow imaging using an iterative
% scheme for active attenuation correction.

% Jelle Plomp. 2024.

% Run this script from step3_load_data_for_analysis such that the following
% inputs are known:

% Inputs:
% Frames - the data, ordered as [x,z,-,time]
% imMask
% std_mm
% frame_plot
% test_plot_regularisers

% Output:
%   Frames_out
%   ... other


% [1]   Attenuation Correction and Normalisation for Quantification of
%       Contrast enhancement in Ultrasound Images of Carotid Arteries, 
%       Cheung et al. 2025, Journal of US in Med & Biol.
%       https://doi.org/10.1016/j.ultrasmedbio.2015.02.010 

% if frames contain negative values, means that hilbert has not been
% applied yet
if min(Frames,[],'all')<0
    Frames = abs(hilbert(Frames));
end
% Apply the mask
Frames_masked = Frames.*imMask;

std_datapoints = (std_mm*1e-3)/ReconInfo_DAS.pas;
filter_size=2*ceil(2*std_datapoints)+1; % See imgaussfilt.m          
filter_half_size = ceil(2*std_datapoints); 

% Find and classify the boundaries of the mask
[~, boundary_class] = find_and_classify_boundaries(imMask, 1);
% Resize the frame stack and mask in case it is too small to apply reflection padding
[Frames_resized, ~,~]= resize_image(Frames_masked, boundary_class,std_datapoints, 1);
[mask_resized, pix_added_top, pix_added_bottom]= resize_image(imMask, boundary_class,std_datapoints, 1);

% Reclassify boundaries (after resizing)
[~, boundary_class] = find_and_classify_boundaries(mask_resized, 1);

% Apply mirroring on the masked image.
I_lumen = vertical_mirrorreflectpadding(Frames_resized, boundary_class, filter_half_size, 1); 
I_lumen_samp = I_lumen(:,:,:,frame_plot);
% Obtain smoothed version of I_lumen
A_lumen = imgaussfilt(I_lumen, std_datapoints, "Padding", "symmetric","FilterDomain","spatial");
%{
% Verify that imgaussfilt works on the frames separately:               
A_lumen_1by1 = zeros(size(I_lumen));
for frame_i=1:size(I_lumen,4)
    A_lumen_1by1(:,:,:,frame_i)=imgaussfilt(I_lumen(:,:,:,frame_i), std_datapoints, "Padding", "symmetric","FilterDomain","spatial");
end
% Check differences
difference = sum(abs(A_lumen-A_lumen_1by1),"all")
%}

% Mask the result and re-mask the lumen
A_lumen = A_lumen.*mask_resized;
mask_nan = double(mask_resized); mask_nan(~mask_resized)=nan;
I_lumen = I_lumen.*mask_nan;
size(I_lumen)

% Go back to original frame size without extra padding
A_lumen = A_lumen(pix_added_top+1:end-pix_added_bottom,:,:,:);
I_lumen = I_lumen(pix_added_top+1:end-pix_added_bottom,:,:,:);
size(I_lumen)

% =============
% Attenuation correction using A_lumen
% ============

if test_plot_regularisers
    [reg_min,NIF, NIF_mean] = find_regulariser_NIF(I_lumen, A_lumen);
    % Plotting
    figure(5);clf(5); hold on
    for frame_i=1:2:min(size(NIF,1),40)
    plot(regularisers, NIF(frame_i,:),'-k', 'HandleVisibility','off')
    
    end
    plot(regularisers, NIF_mean, '-r','LineWidth',2, 'DisplayName', 'Mean NIF over all frames')
    plot(regularisers_interp, NIF_mean_interp, '-g','LineWidth',2, 'DisplayName', 'Interpolated mean NIF')
    ylabel("NIF")
    xlabel("Regulariser"); grid on; 
    set(gca,'FontSize',14)
else
    reg_min = 0.001;
end
% Use the best regulariser for the final calculation
CI_lumen = I_lumen./(A_lumen+reg_min);



% Boundary extension: make sure A_lumen can also be applied outside of
% the lumen
[~, boundary_class] = find_and_classify_boundaries(imMask, 0);

[A_full, A_mask]= extend_A_lumen_boundaries(A_lumen,imMask,boundary_class);
% Calculate mode (most occuring value in the lumen)
% normalisation by the mode (most occuring value) of all frames (was
% not specified by cheung if they do this per frame or over all frames)
bin_edges = 0:0.05:max(CI_lumen(~isinf(CI_lumen)),[],'all')+0.05;
[histcounts_N,~] = histcounts(CI_lumen, bin_edges);

[~, i_max] = max(histcounts_N);
I_mode = mean([bin_edges(i_max) bin_edges(i_max+1)]);


% Corrected intensity and normalised corrected intensity
CI = Frames./(A_full+reg_min);

CI_norm = CI./I_mode;

%% Plotting of al steps
figure(6);clf(6);
subplot(6,2,1)
imagesc(Frames(:,:,:,frame_plot)); colormap(gca,"gray"); title("Frame")
subplot(6,2,2)
imagesc(boundary_class); title("Classify boundaries"); 
subplot(6,2,3)
imagesc(Frames_resized(:,:,:,frame_plot));colormap(gca,"gray"); title("Frame resized")
subplot(6,2,5)
imagesc(I_lumen_samp); colormap(gca,"gray"); title("I_{lumen} (incl. vertical reflect)")
subplot(6,2,7)
imagesc(A_lumen(:,:,:,frame_plot)); colormap(gca,"gray"); title("A_{lumen} = Gaussian filtered I_{lumen}")
subplot(6,2,8)
imagesc(A_full(:,:,:,frame_plot)); colormap(gca,"gray"); title("A_{full} (boundary values extended)")
subplot(6,2,9)
imagesc(CI(:,:,:,frame_plot)); colormap(gca,"gray"); title(strcat("CI  = Frame./(A_{full} + ", num2str(reg_min), ")"))
subplot(6,2,11)
imagesc(CI_norm(:,:,:,frame_plot)); colormap(gca,"gray"); title("Normalised CI = CI/I_{mode}")
subplot(6,2,10)
histogram(CI_lumen, bin_edges); title("Histogram CI_{lumen} (all frames)")
subplot(6,2,12)
CI_lumen_norm = CI_lumen./I_mode;
histogram(CI_lumen_norm, bin_edges); title("Histogram normalised CI_{lumen} (all frames)")
sgtitle(strcat("Example results for frame ", num2str(frame_plot)))
for i=[1,3,5,7,8,9,11]
    subplot(6,2,i);
    if i==9
        clim([min(CI_lumen,[],'all') max(CI_lumen,[],'all')])
    elseif i==11
        clim([min(CI_lumen,[],'all') max(CI_lumen,[],'all')]./I_mode)
    else
        clim([min(I_lumen,[],'all') max(I_lumen,[],'all')])
    end
end

Frames_out = CI_norm;
