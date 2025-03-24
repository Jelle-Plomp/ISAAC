% This script is used to set a media object.
% The object and function movePoints_adapted.m were made to resembe moving
% bubbles in a tube, with a small part of the bubbles being in the acoustic
% shadow of the plaque. The acoustic shadow is just made by setting the
% reflectivity of points in that area lower, so it is not a physical model.

% Media points set by [x,y,z,reflectivity]
% Units = nr of wavelengths from the centre of the transducer

Media_att_region = [-20, 20];
% Now we define a very wide bar, to represent our microbubbles in a tube.
MP_vessel = rand(4000,4); % For all points, we have random speckle, but for a bar of size 2*halfwidth, we 
MP_vessel(:,1) = 2*80*(MP_vessel(:,1)-0.5); % Shift mean x position to zero
MP_vessel(:,2) = 0; % y=0 so everything in-plane
MP_vessel(:,3) = ((MP_vessel(:,3)*6)+ 23)*P.endDepth/40; % scale values in range [0,1] to new range (8mm height, at position 26 mm)

MP_vessel(:,4) = 0.04*MP_vessel(:,4) + 0.08;  % Random amplitude

% Media.MP = vertcat(MP_plaque, MP_vessel);
Media.MP = MP_vessel;
% Now lets add some scatterers below the block of scatterers
MP_i = size(Media.MP, 1);
Media.MP(MP_i+1,:) = [0,0, 0.9*P.endDepth, 1];
Media.MP(MP_i+2,:) = [-20,0, 0.9*P.endDepth, 1];
Media.MP(MP_i+3,:) = [-15,0, 0.95*P.endDepth, 1];
Media.MP(MP_i+4,:) = [-60,0, 0.9*P.endDepth, 1];
Media.MP(MP_i+5,:) = [60,0, 0.9*P.endDepth, 1];

% Set attenuation and function for moving the point.
Media.attenuation = -0.5;
Media.function = 'movePoints_adapted_t1';