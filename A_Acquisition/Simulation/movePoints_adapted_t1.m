function movePoints_adapted_t1
%
% Copyright 2001-2017 Verasonics, Inc.  All world-wide rights and remedies under all intellectual property laws and industrial property laws are reserved.  Verasonics Registered U.S. Patent and Trademark Office.
% Note: in 2023, Verasonics announced Permissive License, see LICENSE file
% in the repository root directory.

% Adapted by Jelle Plomp, Nov. 2022.

% movePoints moves the x location of points in the Media.Pts array in a
% sinusoidal pattern, and changes the reflectivity based on their position.
% This is meant to model an acoustic shadow (as an alternative for
% modelling actual attenuation)

persistent tval
delta = 2*pi/30; % change in tval between successive calls
amp = 10;        % amplitude of displacement in wavelengths

if isempty(tval)
    tval = 0;
end

if evalin('base','exist(''Media'',''var'')')
    Media = evalin('base','Media');
else
    disp('Media object not found in workplace.');
    return
end

if evalin('base','exist(''Media_att_region'',''var'')')
    Media_att_region = evalin('base','Media_att_region');
    x_left = Media_att_region(1);
    x_right = Media_att_region(2);
else
    disp('Attenuation region object not found in workplace.');
    return
end

% Compute amount of change in x position from previous.
old = amp*sin(tval);
tval = tval + delta;
if tval >= 2*pi
    tval = tval - 2*pi;
end
new = amp*sin(tval);
dif = new - old;

% Modify x position of all media points
Media.MP(:,1) = Media.MP(:,1) + dif;

% Modify reflectivity of all media points
Media.MP(:,4) = 0.04*rand(1,size(Media.MP,1)) + 0.08; % Reset with random values

in_region = Media.MP(:,1) > x_left & Media.MP(:,1) <x_right;
Media.MP(in_region,4) = Media.MP(in_region,4) *0.1;
    
assignin('base', 'Media', Media);

return