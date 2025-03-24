function ReconVars = reconstruct_DAS_get_vars(Trans, ReconInfo_DAS, Receive, Resource, TW, size_RFdata, max_angle, use_raddec, att)
%RECONSTRUCT_DAS_GET_VARS calculates all delays, attenuation compensation,
% multiplication factors and selection indices that remain constant over 
% all frames.

% Written by Jelle Plomp (10-01-2023)
% Algorithm based on previous code by Guillaume Lajoinie and others.

% The code was written for a linear transducer and using only 1 angle and
% 1 aperture.
%   Inputs:
%   Verasonics structs: Trans, Receive, Resource, TW
%   ReconInfo_DAS: struct containing relevant info for the DAS
%   reconstruction
%   size_RFdata = Npt x N elements
%   max_angle = maximum angle (degrees) bewteen pixel and element that for which element data is used to reconstruct a pixel.
%   use_raddec % Use decay (compensating for radial scattering 1/r)? (true/false)
%   att = attenuation coefficient.

tic
%% Get parameters related to transducer frequency and sampling frequency
Npt = Receive.endSample-Receive.startSample+1;
f = (Trans.frequency*1e6); %frequency in Hz 
% Note that we need Trans.frequency here, regardless of whether it
% matches the actual frequency. This is because we use this f to
% calculate the position of the elements, which Verasonics calculates
% using Trans.frequency.
lambda = ReconInfo_DAS.c/f; %wavelength
Fs = Receive.decimSampleRate; % decimSampleRate is set by software.

tline = 0:1/Fs:(Npt-1)/Fs; % [Jelle] Timestamp of the samples
tline = tline.*1e-6; % Microseconds to seconds (Since Fs is in MHz)

% Create x and z coordinate arrays for reconstruction
x = ReconInfo_DAS.x_start:ReconInfo_DAS.pas:ReconInfo_DAS.x_end;
z = ReconInfo_DAS.z_start:ReconInfo_DAS.pas:ReconInfo_DAS.z_end;

x = x';
z = z';

lx = length(x);
lz = length(z);

xel = Trans.ElementPos(Receive.aperture:Receive.aperture+Resource.Parameters.numRcvChannels-1,1).*lambda;
zel = Trans.ElementPos(Receive.aperture:Receive.aperture+Resource.Parameters.numRcvChannels-1,3).*lambda;
%% dt correction [Jelle]: correction for travelling through lens, and something with waveform position

dt_lens_corr = 2*Trans.lensCorrection/(f); 
% Trans.lensCorrection is in terms of wavelengths (so assumes
% Trans.frequency to be the frequency)
dt_peak = TW(1).peak/f; % Distance to peak of waveform
dt_corr = dt_lens_corr + dt_peak;
tline = tline + dt_corr;
    
    
%% Attenuation correction and a frequency axis for Fourier Filtering
ReconVars.att_vec = 10.^(att.*tline.*ReconInfo_DAS.c.*f.*1e-6.*1e2./20);
ReconVars.f_axis = linspace(0,Fs,Npt); %dit is de frequentie "as" bij een FT

%% Calculating variables related to selection of datapoints from RF data
d_origin_to_z = z; % 1xlz, for each z the distance to the origin (=z)

%Preallocation
[ReconVars.mult_factor, ReconVars.RF_loc_ind] = deal(NaN(length(xel),lx, lz));  %  Nelements x len(x) x len(z) 
ReconVars.count = zeros(lx, lz);

% Calculations
for j =  1:lz  
    for i =  1:lx  
        % Probably still some speed to gain by using matrcies rather
        % than for-loop.

        % Distance from this pixel to each of the elements:
        d_pixel_to_el = sqrt( (x(i)-xel).^2 + (z(j)-zel).^2  );
        %time of signal from sending to receiving:
        deltt=(d_origin_to_z(j)+d_pixel_to_el)./ReconInfo_DAS.c;
        % Angle between z-axis and the line connecting pixel and
        % element:
        angle_to_pixel = atand(abs(xel-x(i))/abs(z(j)));
        % Distance in number of points of returning signal per element:
        NN= round((deltt + dt_corr ) *Fs*1e6 );
        % Indicate whether the value at this NN index should be used (per
        % element)
        % The inner "and" checks if NN index is higher than 1 and shorter than the signal length
        % The outer "and" determines if the angle between pixel and element
        % is within the specified limit max_angle.
        isvalid= and(and(NN <= Npt, NN > 1), angle_to_pixel<=max_angle);
        % If radial decay is used, the isvalid (0/1) values are multiplied
        % by the distance from a pixel to each element. If no decay is
        % used, the isvalid values are used directly.
        
        % isvalid ensures that any non-valid value is multiplied by 0 so doesn't contribute
        % to the signal.
        if use_raddec
            % d_pixel_to_el is used to account for radial scattering, which
            % decays with 1/r (adopted from macy's code)
            ReconVars.mult_factor(:,i,j) = isvalid.*d_pixel_to_el;
        else
            ReconVars.mult_factor(:,i,j) =  isvalid;
        end
        ReconVars.count(i,j) = sum(isvalid); % The number of datapoints used to recosntruct this pixel
        % Indicate which datapoints from RF_loc to use
        % Here size_RFdata should have the same size as RF_loc in the 
        % reconstruction loop
        ReconVars.RF_loc_ind(:,i,j) = sub2ind(size_RFdata, NN, [1:size_RFdata(2)]');
    end
end
ReconVars.lx = lx;
ReconVars.lz = lz;
ReconVars.x = x;
ReconVars.z = z;
ReconVars.Npt = Npt;
ReconVars.xel = xel;
ReconVars.zel = zel;
toc
end

