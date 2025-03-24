%{
Main script for acquisition using the iterative scheme for active
attenuation correction.

If you use this software, please consider citing the paper
"Imaging Behind the Plaque: Improved Blood flow Quantification Using an
Iterative Scheme for Active Attenuation Correction" - Ultrasound in
Medicine and Biology 2025, Jelle Plomp, Ashkan Ghanbarzadeh-Dagheyan,
Michel Versluis, Guillaume Lajoinie and Erik Groot Jebbink.

Jelle Plomp 2024.
%}

clear all; close all; clc

%% Parameters
simulateMode = 1; % Simulation mode?
% File names and saving
filename = 'L12-3v_attcorr.mat'; % Should be called filename, this is recognized by VSX.
workfolder = 'C:\Users\local.la\Documents\Vantage-4.2.0-2001220500\Jelle\Git_publication'; % Directory where the repo is cloned/downloaded

% =========================================================================
% General
% =========================================================================
% Some inputs (used in the filenames when saving data)
P.Tube_nr = input('Tube number: ', 's');
P.Model_nr = input('Model number: ','s');
% Imaging parameters
P.startDepth = 0;   % Acquisition depth in wavelengths
P.endDepth_mm = 40; % The end depth to be defined in [mm]
P.TransFrequency_true = 7.3; % true transmit frequency in MHz, trans.frequency is 2x higher to increase sample rate
P.Trans_maxhighvoltage = 40; % max allowed voltage for live view and premeasurement

% =========================================================================
% HFR
% =========================================================================
P.PRP_HFR = 600; %  time (in microseconds) between frames in microseconds 111--> 9000 HZ, 167  -->11 6000Hz
if simulateMode 
    P.numAcqsSuperFrame = 10;
    P.HFRBufFrames = 2;
else
    P.numAcqsSuperFrame = 100; % no. of Acquisitions in a Receive frame (this is a "superframe")
    P.HFRBufFrames = 30;
end


% =========================================================================
% Live view 
% =========================================================================
P.PRP_Live = 10000; % 10 ms 
P.LiveRcvBufFrames = 200; 

% Reconstruction parameters (used in reconDAS (external processing
% function))
P.ReconResolution = 0.1*1e-3;
P.Recon_x_start = -15e-3;
P.Recon_x_end = 15e-3;
P.Recon_z_start = 10*1e-3;
P.Recon_z_end =40*1e-3;

% =========================================================================
% Premeasurement (=Iterative Procedure) (uses same framerate as liveview)
% =========================================================================
PreM.method = 'rayleigh'; % Only rayleigh is allowed. This refers to the 
% fact that the published version of ISAAC uses a linear acoustic model
% based on a rayleigh(-sommerfeld) model.
PreM.GUI = 'extended'; % default OR extended
% PreM definition (continues later in script after defining ReconInfo_DAS)

PreM.apod_init = 0.2; % Current version of ISAAC assumes this is the case.

PreM.std_mm = 1; % standard deviation in mm of the gaussian filter used to smooth results.

% PreM parameters related to the stopping criterium
PreM.usestopping = false; % If false, parameters below are not used and user has to stop premeasurement manually.
if PreM.usestopping
    error('Stability criterium was not implemented. Options below are examples be used for defining the PreM_checkCondition function.') 
    % Stability criterium (early stopping (ES)):
    PreM.ES.method = 'apod'; % 'apod' or 'roi' (which profile is used to check stability)
    PreM.ES.Nit = 50; % Number of stable iterations after which algorithm stops even if condition is not met.
    PreM.ES.maxchange = 0.04; % 0.1 = 10% --> less than 10% change counts as stable
end
PreM.ES.Nit_max = 50; % Maximum number of iterations after which premeasurements always stops
P.PreMBufFrames = PreM.ES.Nit_max; 

% Stopping criterium (not implemented or used):
PreM.condition = nan;
PreM.condition_arg = nan;

% PI controller
PreM.PI_contrl.Kp = 20;
PreM.PI_contrl.Ki = 100;
PreM.PI_contrl.sp = NaN; % Set-point value (i.e. target image value). Will be determined based on first measurement.
PreM.PI_contrl.sp_version = 'median_high10perc'; % 'median_high10perc'
PreM.PI_contrl.errormethod = 'clipintegralerror'; % Has to be 'clipintegralerror' since PreM.method=='rayleigh' only takes this.
% Enable saving with partial compensation?
PreM.enablePartialCompensation = false;
if PreM.enablePartialCompensation
    N_TX = 7; % In this case more TX profiles are set such that in each TX a different TX.Apod profile can be stored.
    PreM.PartialCompensation.BufferNumbers = [3,4,6:N_TX+3];
else
    N_TX = 2;
end

if strcmp(PreM.method,'rayleigh')
    PreM.rayleigh.weight = 1000; % The regularisation weight used in the solver.
end

% =========================================================================
% Other (no need to adapt - except if you want to change the angle)
% =========================================================================

% fps calculation assuming 1 angle and 1 aperture
P.fps = 10^6/(P.PRP_HFR);
acqtime = P.HFRBufFrames*P.numAcqsSuperFrame*1/P.fps;
fprintf('Aquiring HFR Frames at %.2g fps for %.2g seconds\n',P.fps, acqtime);

% Script is only ready for use with 1 angle, but we can still change the
% angle here:
dtheta=0; % Not in use
startAngle=0;

% This we do not really use since I assumed 1 angle and 1 aperture. 
np = 1;                 % I have assumed for the rest of the script that this will remain 1.
na = 1;                 % number of angles: 1, 3, 5 or 7
P.PRP_AA = 500; % Time between acquisition of angles and apertures (also in live view)
P.na=na;
P.np=np;

if P.na ~= 1 || P.np ~= 1
    error("This script is not ready (or meant to be used) for multiple angles/apertures")
end
% Empty struct to keep track of processing times
ProcTimes=struct();
%% Folders for data saving, and adding the correct paths
cd(workfolder)
if ~strcmp(workfolder(end),'\')
    workfolder= [workfolder '\'];
end
% Create Data and MatFiles folders 1 level above the working folder
P.savefiledir = [workfolder 'Data\'];
matfilesavedir = [workfolder 'MatFiles\'];
% Add the paths in A_acquisition and delay-and-sum-reconstruction
allpaths = genpath(workfolder);
allpaths = strsplit(allpaths, ';');

for i=1:length(allpaths)
    if contains(allpaths{i},'A_Acquisition') || contains(allpaths{i},'delay-and-sum-reconstruction')
        addpath(allpaths{i}) 
    end
end

% Create the savefiledir for Data and MatFiles 
if ~isfolder(P.savefiledir)
    mkdir(P.savefiledir) 
end
addpath(P.savefiledir)
if ~isfolder(matfilesavedir)
    mkdir(matfilesavedir)
end
addpath(matfilesavedir)
clearvars allpaths
%% Define system parameters.
Resource.Parameters.numTransmit = 128;  % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;  % number of receive channels.
Resource.Parameters.speedOfSound = 1480;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.connector = 1;
Resource.Parameters.simulateMode = simulateMode;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

%% Specify Trans structure array.
Trans.name = 'L12-3v';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
Trans.frequency = 2*P.TransFrequency_true;
Trans = computeTrans(Trans);  % L12-3v transducer is 'known' transducer so we can use computeTrans.
% Trans.lensCorrection = 6; % In wavelengths. In computeTrans, the value of LensCorrection is 1.183 mm
Trans.maxHighVoltage = P.Trans_maxhighvoltage;      % set maximum high voltage limit for pulser supply.

WL = Trans.spacingMm/Trans.spacing; % Wavelength in [mm]
P.endDepth = ceil((P.endDepth_mm/WL)/64)*64; % endDepth in number of wavelength. This should preferrably be a multiple of 128 samples.

P.SampleRate = 250/round(250/(Trans.frequency*4));  % Real sample rate to be use to fill the buffer depth
P.SamplePerWave = P.SampleRate/Trans.frequency;  %Number of samples per wave
%P.AnglesSensitivity = (90/(pi/2))*linspace(-pi/2,pi/2,length(Trans.ElementSens));

%% Media (from Simulation folder)
if simulateMode ==1
    media_object_Jelle_t1
end
%% TPC profiles
TPC(1).name = 'TPC general';
TPC(1).maxHighVoltage = P.Trans_maxhighvoltage; 

%% Specify PData structure array.
% A Dummy PData struct. Since we do not use VSX display, this is not in
% use. Without it however, Matlab crashes.
PData.PDelta = [Trans.spacing, 0, 0.5];
PData.Size(1) = 1; 
PData.Size(2) = 1;
PData.Size(3) = 1;      
PData.Origin = [0,0,0]; 


%% VSX calculations
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
P.LengthDiff = maxAcqLength - P.startDepth;
P.OriginalEndReceiveSample =  2*P.LengthDiff*P.SamplePerWave; % It is multiplied by 2 to capture the go & return length
P.division = P.OriginalEndReceiveSample/128; % check to see whether the required number of samples are in order of 128 samples block.
P.VSXEndReceiveSample = ceil(P.division)*128; % Required number of samples (by Verasonics) per acquasition to capture the new maximum depth --> number of samples in buffer for each acquasition will be automatically replaced by a value higher than the original sample number, which is dividable by 128
P.VSXmaxAcqLength = (P.VSXEndReceiveSample/(2*P.SamplePerWave) ) + P.startDepth; % Computing the maximum depth that will be captured by Verasonics
P.VSXendDepth = sqrt( P.VSXmaxAcqLength^2 - ((Trans.numelements-1)*Trans.spacing)^2);

%% Specify Resources.

% Buffer 1 - Live view 
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = P.VSXEndReceiveSample*np*na;   % this size allows for maximum range
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = P.LiveRcvBufFrames;       % number of 'super frames'

% Image buffer for live view reconstructions
Resource.ImageBuffer(1).numFrames = 10;
Resource.ImageBuffer(1).rowsPerFrame = length(P.Recon_z_start:P.ReconResolution:P.Recon_z_end);
Resource.ImageBuffer(1).colsPerFrame = length(P.Recon_x_start:P.ReconResolution:P.Recon_x_end);

% Rcvbuffer and imagebuffer 2 -> for frames of premeasurement
Resource.RcvBuffer(2).datatype = 'int16';
Resource.RcvBuffer(2).rowsPerFrame =  P.VSXEndReceiveSample*np*na;
Resource.RcvBuffer(2).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(2).numFrames = P.PreMBufFrames;     

Resource.ImageBuffer(2).numFrames = P.PreMBufFrames;
Resource.ImageBuffer(2).rowsPerFrame = length(P.Recon_z_start:P.ReconResolution:P.Recon_z_end);
Resource.ImageBuffer(2).colsPerFrame = length(P.Recon_x_start:P.ReconResolution:P.Recon_x_end);


% Rcvbuffer 3 -> HFR saving 
% Rcvbuffer 4 -> HFR saving of current voltage with TX(2) (equal
% apodisation)
% Rcvbuffer 5 -> HFR saving of initial voltage with TX(2) (equal
% apodisation)
for i=[3,4,5]
    Resource.RcvBuffer(i).datatype = 'int16';
    Resource.RcvBuffer(i).rowsPerFrame = P.VSXEndReceiveSample*P.numAcqsSuperFrame*np*na;
    Resource.RcvBuffer(i).colsPerFrame = Resource.Parameters.numRcvChannels;
    Resource.RcvBuffer(i).numFrames = P.HFRBufFrames;     %  frames for rcvDataLoop buffer.
end

if PreM.enablePartialCompensation
    for i = PreM.PartialCompensation.BufferNumbers
        Resource.RcvBuffer(i).datatype = 'int16';
        Resource.RcvBuffer(i).rowsPerFrame = P.VSXEndReceiveSample*P.numAcqsSuperFrame*np*na;
        Resource.RcvBuffer(i).colsPerFrame = Resource.Parameters.numRcvChannels;
        Resource.RcvBuffer(i).numFrames = P.HFRBufFrames;     %  frames for rcvDataLoop buffer.
    end
end

% check whether the chosen buffer sizes are suitable for Verasonics
for bufnum=1:length(Resource.RcvBuffer)
    if Resource.RcvBuffer(bufnum).rowsPerFrame > 524288
        disp(['Error. Too many samples for buffer ' num2str(bufnum)])
    else
        disp(['Buffer number ' num2str(bufnum) ' good to go'])
    end
end
Resource.VDAS.dmaTimeout = 3000;
%% Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [P.TransFrequency_true,0.67,2,1];

%% TX structure array.

TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'aperture', 33, ... % Changed this to 33 such that we use the middle 128 elements.
                   'Apod', PreM.apod_init.*ones(1,Resource.Parameters.numTransmit), ...
                   'focus', 0.0, ...
                   'Steer', [startAngle,0.0], ...
                   'Delay', zeros(1,Resource.Parameters.numTransmit)), 1, N_TX); % Last value was np*na for working with multiple apertures/angles

% - Set event specific TX attributes.
% If multiple apertures/angles have to be used, this has to be adapted. See
% previous script versions.
for TX_i=1:length(TX)
    TX(TX_i).Delay = computeTXDelays(TX(TX_i)); % non-inverted transmits 
end


%% Specify TGC Waveform structure.
% TGC.CntrlPts = [0,298,395,489,618,727,921,1023];
TGC.CntrlPts = [0,1023,1023,1023,1023,1023,1023,1023];

TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

%% Specify Receive structure arrays -

% Note: Apodization in the Kolo probe is not in the same fashion as regular
% linear probe line L11-4. Since we can not fire all 256 elements at the
% same time, we have to define aperture. By doing so, we can only use 128
% elements of Kolo at each time. Therefore, our
% "Resource.Parameters.numTransmit" is 128 (instead of ideally being 256).
% As a consequence of this CMUT probe coding, the available apodication
% will also be 128 (instead of ideally beng 256). Therefore, we have to
% define appodization per each apperture as follows:

% Jelle:
% Using middle aperture:
% aperture should be 33, and set apodization to 1 for all following 128
% elements rather than only 96 elements.
nr_HFRbufs = 2;


% One receive element used for each frame
%Live view + Premeasurement + HFR saving  (assuming na=1 and np=1) 

length_Receive = P.LiveRcvBufFrames + P.PreMBufFrames + nr_HFRbufs*P.HFRBufFrames*P.numAcqsSuperFrame;
if PreM.enablePartialCompensation
    length_Receive = P.LiveRcvBufFrames + P.PreMBufFrames + N_TX*P.HFRBufFrames*P.numAcqsSuperFrame;
end

clearvars nr_HFRbufs
Receive = repmat(struct('Apod', ones(1,Resource.Parameters.numTransmit), ... % all (used) elements apod 1
                        'aperture', 33, ...
                        'startDepth', P.startDepth, ...
                        'endDepth',P.VSXmaxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 1), 1,length_Receive);      % movepoints EVERY acquisition to illustrate superframe concept
                                                       % real-time images will look "jerky" but using the reconstructAll script,
                                                                                                                     % playback process all acquisitions and shows smooth displacement                                                                   
% - Set event specific Receive attributes --> Buffer 1, Live view
if np == 1
    track = 1;
    for i = 1:P.LiveRcvBufFrames
        for j=1:P.na % Not sure if angles will work, did not really focus on this yet. So assume P.na=1 for now.
            rcvNum_buffer1 = na*(i-1) + j;
            Receive(rcvNum_buffer1).bufnum = 1; %
            Receive(rcvNum_buffer1).framenum = i;
            Receive(rcvNum_buffer1).acqNum = j;
            track = track + 1;
        end
    end
    
    rcvNum_prem_start = track-1;
    
    
    % Premeasurement - goes to buffer 2
    for i = 1:P.PreMBufFrames   
        for j = 1:P.na % na acquisitions per frame
            rcvNum = rcvNum_prem_start+P.na*(i-1) + (j);
            Receive(rcvNum).framenum = i;
            Receive(rcvNum).acqNum = j;
            Receive(rcvNum).bufnum = 2;
            track=track+1;
        end  
    end
    % - Set event specific Receive attributes --> Buffer 3, Saving HFR
    rcvNum_HFR_start = track-1; %start value for receive values saved to buffer 3
    for frame = 1:P.HFRBufFrames
        for acqNum = 1:P.numAcqsSuperFrame
            for angNum = 1
                % -- Acquisitions for 'super' frame and angle.
                rcvNum = P.numAcqsSuperFrame*P.na*(frame-1) + P.na*(acqNum-1)+ angNum + rcvNum_HFR_start;
            
                Receive(rcvNum).framenum = frame;
                Receive(rcvNum).acqNum = na*(acqNum-1)+angNum;
                Receive(rcvNum).bufnum = 3; % HFR saves to buffer 3
                track=track+1;
            end
        end
    end
    % Buffer 4 - will capture data obtained with TX2 (equal apodisation) at
    % current voltage
    rcvNum_HFR_D2_TX2_Vend = track-1;
    for frame = 1:P.HFRBufFrames
        for acqNum = 1:P.numAcqsSuperFrame
            for angNum = 1
                % -- Acquisitions for 'super' frame and angle.
                rcvNum = P.numAcqsSuperFrame*P.na*(frame-1) + P.na*(acqNum-1)+ angNum + rcvNum_HFR_D2_TX2_Vend;
            
                Receive(rcvNum).framenum = frame;
                Receive(rcvNum).acqNum = na*(acqNum-1)+angNum;
                Receive(rcvNum).bufnum = 4; % HFR saves to buffer 4
                track=track+1;
            end
        end
    end
    
    if PreM.enablePartialCompensation
        % Buffers for partial compensation
        % Below should be arrays of same length as TX.
        
        PreM.PartialCompensation.RcvStartNums = [rcvNum_HFR_start, rcvNum_HFR_D2_TX2_Vend, nan(1,N_TX-2)];
        % 
        for TX_i=3:N_TX % First two are already assigned.
            rcvNumStart = track-1;
            PreM.PartialCompensation.RcvStartNums(TX_i) = rcvNumStart;
            for frame = 1:P.HFRBufFrames
                for acqNum = 1:P.numAcqsSuperFrame
                    for angNum = 1
                        % -- Acquisitions for 'super' frame and angle.
                        rcvNum = P.numAcqsSuperFrame*P.na*(frame-1) + P.na*(acqNum-1)+ angNum + rcvNumStart;
                    
                        Receive(rcvNum).framenum = frame;
                        Receive(rcvNum).acqNum = na*(acqNum-1)+angNum;
                        Receive(rcvNum).bufnum = PreM.PartialCompensation.BufferNumbers(TX_i); % HFR saves to this buffer 
                        track=track+1;
                    end
                end
            end
        end
    end
else    
   % Jelle: did not implement the multple aperture thing. See Ashkans
   % script in measurements of 12-12-2022 if I want to re-introduce this.
   % Note that that requires more rows in the Receive struct.
end


%% ReconInfo
% Here I added some new fields not originally present in the ReconInfo
% object, but which can be used in the external reconDAS function
% Basically, every parameter needed for DAS reconstruction that is not defined 
% in % the reconDAS function.
%         ReconInfo = struct('mode', 'replaceIntensity', ...
%                            'txnum', 1, ... 
%                            'rcvnum', 1, ...
%                            'regionnum', 1);
ReconInfo_DAS = struct('c', Resource.Parameters.speedOfSound, ...
                        'pas', P.ReconResolution, ...
                        'x_start', P.Recon_x_start, ...
                        'x_end', P.Recon_x_end, ...
                        'z_start', P.Recon_z_start, ...
                        'z_end', P.Recon_z_end,...
                        'first_frame', true, ...  % Will be changed to false after the first frame
                        'Nfrs_avg', 5,...
                        'demod_rf', 1,... % Whether or not hilbert transform will be applied for the live view.
                        'demod_rf_method', 'rf2iq'); % rf2iq or rf_hilb (recommended: rf2iq) 
displayInfo.dynRange = 55;

%% History and variables for premeasurement 
% (so we do not have to evaluate all variables seperately in premeasurement)
PreM.x = ReconInfo_DAS.x_start:ReconInfo_DAS.pas:ReconInfo_DAS.x_end;
PreM.z = ReconInfo_DAS.z_start:ReconInfo_DAS.pas:ReconInfo_DAS.z_end;
PreM.ReconInfo_DAS = ReconInfo_DAS;

% history initialisation
History(1).Voltage = NaN;
History(1).Apod = NaN;
History(1).Imax_ROI = NaN;
History(1).I_ROI_line = NaN;
History(1).I_L_used = NaN;
History(1).P = NaN;
History(1).I = NaN;
History(1).I_clip = NaN;
History(1).error_total = NaN;
History(1).el_clip = NaN;
History(1).Nit_stable = 0;
%% Specify Process structure array.
                     
% % timer
P_nr_timer_start = 1;
Process(P_nr_timer_start).classname = 'External';
Process(P_nr_timer_start).method = 'TimerInit';
Process(P_nr_timer_start).Parameters = {'srcbuffer','none',... 
                                         'srcbufnum','none',...
                                         'dstbuffer','none'};              
P_nr_timer_end = 2;                     
Process(P_nr_timer_end).classname = 'External';
Process(P_nr_timer_end).method = 'PrintTime';
Process(P_nr_timer_end).Parameters = {'srcbuffer','none',...  % name of buffer to process.
                                      'srcbufnum','none',...
                                      'dstbuffer','none'};  
% Premeasurement (=iterative procedure) - main function
P_nr_PreM = 3;
Process(P_nr_PreM).classname = 'External'; 
Process(P_nr_PreM).method = 'PreM_main';
Process(P_nr_PreM).Parameters = {'srcbuffer', 'image',...
                                     'srcbufnum', 2,... 
                                     'srcframenum', -1,... % Take the last reconstructed frame in buffer 2
                                     'dstbuffer', 'none'}; 
% Call DAS reconstruction during live view
P_nr_DASlive = 4;
Process(P_nr_DASlive).classname = 'External';
Process(P_nr_DASlive).method = 'reconDAS_L12_3v_fast';
Process(P_nr_DASlive).Parameters = {'srcbuffer', 'receive',... % We want to reconstruct the received data
                         'srcbufnum', 1,... % Buffer 1 is where our live view data is
                         'srcframenum',-1, ... % Last frame in the buffer
                         'dstbuffer', 'image',... % Reconstructed image should go to imgae buffer
                         'dstbufnum', 1,...
                         'dstframenum', -1}; % Not sure if this should be 1 or -1 or something else entirely
% Call DAS reconstruction during premeasurement
P_nr_DASpreM = 5;
Process(P_nr_DASpreM).classname = 'External';
Process(P_nr_DASpreM).method = 'reconDAS_L12_3v_fast';
Process(P_nr_DASpreM).Parameters = {'srcbuffer', 'receive',... % We want to reconstruct the received data
                         'srcbufnum', 2,... % Buffer 2 is where our premeasurement data is
                         'srcframenum',-1, ... % Last frame in the buffer
                         'dstbuffer', 'image',... % Reconstructed image should go to imgae buffer
                         'dstbufnum', 2,...
                         'dstframenum', -1}; % Not sure if this should be 1 or -1 or something else entirely

% Display image (including log compression) (only used for live view)
P_nr_display = 6;
Process(P_nr_display).classname = 'External';
Process(P_nr_display).method = 'displayLogComp';
Process(P_nr_display).Parameters = {'srcbuffer', 'image',... % We want to display the reconstructed data 
                         'srcbufnum', 1,... % Image Buffer 1 is where the live view data is stored
                         'srcframenum',-1,... % Last frame in the buffer
                         'dstbuffer','none'};
% save HFR data --> The "srcframenum" value was introduced to be zero as suggested by Verasonics team
% Accrording to them, this would solve the following problem:
% When some superframe's RF data are empty in the RCVData
P_nr_saveHFR = 7;
Process(P_nr_saveHFR).classname = 'External';
Process(P_nr_saveHFR).method = 'SaveRcvDataHFR_L12_3_TX1';
Process(P_nr_saveHFR).Parameters = {'srcbuffer','receive',...  
                         'srcbufnum',3,... % Buffer 3
                         'srcframenum',0,...        
                         'dstbuffer','none'};
% Saving for saving 2 mesurements subsequently.                     
P_nr_saveHFR_D1_TX1_Vend = 8;
Process(P_nr_saveHFR_D1_TX1_Vend).classname = 'External';
Process(P_nr_saveHFR_D1_TX1_Vend).method = 'SaveRcvDataHFR_L12_3_TX1_and_TX2';
Process(P_nr_saveHFR_D1_TX1_Vend).Parameters = {'srcbuffer','receive',...  
                         'srcbufnum',3,... % Buffer 3
                         'srcframenum',0,...        
                         'dstbuffer','none'};                    
P_nr_saveHFR_D2_TX2_Vend = 9;
Process(P_nr_saveHFR_D2_TX2_Vend).classname = 'External';
Process(P_nr_saveHFR_D2_TX2_Vend).method = 'SaveRcvDataHFR_L12_3_TX1_and_TX2';
Process(P_nr_saveHFR_D2_TX2_Vend).Parameters = {'srcbuffer','receive',...  
                         'srcbufnum',4,... % Buffer 4
                         'srcframenum',0,...        
                         'dstbuffer','none'};  
                     
% Saving of premeasurement image data
P_nr_PreMsave = 10;
Process(P_nr_PreMsave).classname = 'External'; 
Process(P_nr_PreMsave).method = 'SavePreM_L12_3v';
Process(P_nr_PreMsave).Parameters = {'srcbuffer', 'image',...
                                     'srcbufnum', 2,... % Image buffer 2 is where the image data of the premeasurement is stored.
                                     'srcframenum', 0,...
                                     'dstbuffer', 'none'}; 
% Reset startEvent
P_nr_resetStartEvent = 11;
Process(P_nr_resetStartEvent).classname = 'External'; 
Process(P_nr_resetStartEvent).method = 'resetStartEvent';
Process(P_nr_resetStartEvent).Parameters = {'srcbuffer', 'none',...
                                     'srcbufnum', 'none',...
                                     'dstbuffer', 'none'};
                                 
% Saving including partial transmits
if PreM.enablePartialCompensation
    % We can use processes 8 and 9 here. the rest we have to redefine.
    PreM.PartialCompensation.ProcessNumbers = nan(1,N_TX);
    for TX_i = 1:length(TX)
        PreM.PartialCompensation.ProcessNumbers(TX_i) = P_nr_resetStartEvent+TX_i;
        Process(P_nr_resetStartEvent+TX_i).classname = 'External';
        Process(P_nr_resetStartEvent+TX_i).method = 'SaveRcvDataHFR_L12_3_PartialCompensation';
        Process(P_nr_resetStartEvent+TX_i).Parameters = {'srcbuffer','receive',...  
                         'srcbufnum', PreM.PartialCompensation.BufferNumbers(TX_i) ,... 
                         'srcframenum',0,...        
                         'dstbuffer','none'}; 
    end
end
    
%% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start.
SeqControl(1).argument = 1;

SeqControl(2).command = 'timeToNextAcq';  % time between superframes in live
SeqControl(2).argument = P.PRP_Live;  % 10 msecs

SeqControl(3).command = 'returnToMatlab';

SeqControl(4).command = 'jump'; %  - Jump back to startEvent without the calcapod(=5).
SeqControl(4).argument = 1; % Note that this will be replaced by event_nr_live when creating the event struct below

% -- Change to TPC Profile 1 (The high voltage parameter will be adapted in calcApod)
SeqControl(5).command = 'setTPCProfile';
SeqControl(5).condition = 'immediate';
SeqControl(5).argument = 1;


% Used in saving HFR images
SeqControl(6).command = 'timeToNextAcq';  % time between acquisitions for HFR saving
SeqControl(6).argument = P.PRP_HFR;  % 1000 usecs --> capturing 1000 fps

SeqControl(7).command = 'sync';
SeqControl(7).argument = 1000000;
SeqControl(8).command = 'noop';
SeqControl(8).argument = 500000;
%%% ============ Not in use =============
SeqControl(9).command = 'timeToNextAcq';  % time between acquisitions for HFR and Live view
SeqControl(9).argument = P.PRP_AA;      % Note: not actually used when using 1 angle and 1 aperture.
%%% ======================================
SeqControl(10).command = 'noop';
SeqControl(10).argument = 10000; % 10000 musec = 10 ms. We use this between acquisitions in the premeasurement.

nsc = length(SeqControl)+1; % nsc is count of SeqControl objects
%% Events (loaded from external scripts)
n = 1; % n is count of Events
% Pre-measurement events:
Events_1_PreM
% Live-view events
% Per default, we start at the following loop
% This acquire all frames defined in RcvBuffer
% Update start event nr. for the live data loop
event_nr_live=n;
Resource.Parameters.startEvent = event_nr_live;
SeqControl(4).argument = event_nr_live;

Events_2_LiveView

% Events: Saving HFR (only using TX1 at the current voltage)
event_nr_HFR = n;

Events_3_SaveHFR_TX1

% Events: Saving HFR (first using TX1 at current voltage, then using TX 2
% at current voltage, then using TX2 at initial voltage)
event_nr_HFR_3times = n;

Events_4_SaveHFR_TX1_and_TX2
% Events: Saving acquired imageData of premeasurement
event_nr_PreM_saving = n;

Events_5_PreM_save_ImageData

% Events: saving HFR from TX1 to N_TX (for partial compensation)
if PreM.enablePartialCompensation
    event_nr_HFR_partial = n;

    Events_6_SaveHFR_withPartialCompensation
end

%% User specified UI Control Elements

% - Selecting current data
% should be a control of type Button
UI(3).Control = {'UserB5', 'Style', 'VsPushButton', 'Label', 'PreM 1/0'};
UI(3).Callback = text2cell('%CalcImageLevelCallback');
UI(4).Control = {'UserB4', 'Style', 'VsPushButton',...
                  'Label', 'Save PreM'};
UI(4).Callback = text2cell('-UI#4Callback');

UI(7).Control = {'UserB3', 'Style', 'VsPushButton',...
                  'Label', 'Save HFR'};
UI(7).Callback = text2cell('-UI#7Callback');

UI(8).Control = {'UserB2', 'Style', 'VsPushButton',...
                  'Label', 'Save HFR 2x'};

UI(8).Callback = text2cell('-UI#8Callback');
if PreM.enablePartialCompensation
    UI(9).Control = {'UserB1', 'Style', 'VsPushButton',...
                      'Label', 'HFR - Part'};
    UI(9).Callback = text2cell('-UI#9Callback');
end

UI(6).Control = {'UserC3','Style','VsSlider',...
    'Label','Dynamic range',... 
    'SliderMinMaxVal', [1.0,100,displayInfo.dynRange],...
    'SliderStep',[0.1 0.5]};
UI(6).Callback = text2cell('%DynamicRange');

% External functions (timer functions. (simple and therefore defined within this script)
EF(1).Function = text2cell('%-EF#1');
EF(2).Function = text2cell('%-EF#2');

%%
% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = P.numAcqsSuperFrame;

% Save all the structures to a .mat file.
save([matfilesavedir,filename]);
filename = [matfilesavedir,filename]; %such that VSX will find it.
disp([ mfilename ': NOTE -- Running VSX automatically!']), disp(' ')
VSX
commandwindow  % just makes the Command window active to show printout

return

%% **** Callback routines to be converted by text2cell function. ****
%-UI#4Callback - Save Premeasurement data
event_nr_PreM_saving = evalin('base','event_nr_PreM_saving');
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Parameters',1,'startEvent',event_nr_PreM_saving};
assignin('base','Control',Control);
return
%-UI#4Callback

%-UI#7Callback - Initiate HFR Acquisition
P = evalin('base','P');
acqtime = P.HFRBufFrames*P.numAcqsSuperFrame*(P.PRP_HFR*1e-6);
fprintf('Aquiring HFR Frames at %.2g fps for %.2g seconds\n',1/(P.PRP_HFR*1e-6), acqtime);
% assignin('base','timer',timer);
event_nr_HFR = evalin('base','event_nr_HFR');
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Parameters',1,'startEvent',event_nr_HFR};
assignin('base','Control',Control);
return
%-UI#7Callback

%-UI#8Callback - Initiate Full HFR Acquisition - TX1 at current voltage,
%TX2 at current voltage and initial voltage.
P = evalin('base','P');
acqtime = P.HFRBufFrames*P.numAcqsSuperFrame*(P.PRP_HFR*1e-6);
fprintf('Aquiring HFR Frames at %.2g fps for %.2g seconds\n',1/(P.PRP_HFR*1e-6), acqtime);
% assignin('base','timer',timer);
event_nr_HFR_3times = evalin('base','event_nr_HFR_3times');
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Parameters',1,'startEvent',event_nr_HFR_3times};
assignin('base','Control',Control);
return
%-UI#8Callback

%-UI#9Callback - Initiate Full HFR Acquisition with partial copmensation
P = evalin('base','P');
acqtime = P.HFRBufFrames*P.numAcqsSuperFrame*(P.PRP_HFR*1e-6);
fprintf('Aquiring HFR Frames at %.2g fps for %.2g seconds\n',1/(P.PRP_HFR*1e-6), acqtime);
% assignin('base','timer',timer);
event_nr_HFR_partial = evalin('base','event_nr_HFR_partial');
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Parameters',1,'startEvent',event_nr_HFR_partial};
assignin('base','Control',Control);
return
%-UI#9Callback

%-EF#1 - Initialize Timer
TimerInit()
tic
return
%-EF#1

%-EF#2 - Print Time Since Timer Initialized
PrintTime()
dt = toc;
fprintf('dt = %f \n', dt)
return
%-EF#2

%DynamicRange
% Adapted to own live view. 
displayInfo = evalin('base','displayInfo');
displayInfo.dynRange = UIValue;
assignin('base','displayInfo',displayInfo);
%DynamicRange

%CalcImageLevelCallback
Resource = evalin('base', 'Resource');
event_nr_live = evalin('base', 'event_nr_live');
current_startevent = Resource.Parameters.startEvent;
control_i = 1;
if current_startevent == 1 % Manual stop premeasurement
    Resource.Parameters.startEvent = event_nr_live;
    enablePartialCompensation = evalin('base', 'PreM.enablePartialCompensation')
    if enablePartialCompensation
        method = evalin('base','PreM.method')
        % Calculate parameters for partial compensation
        TX = evalin('base', 'TX');
        if strcmp(method, 'rayleigh')
                TX = calcTX_PartialCompensation_v_rayleigh(TX);
        end
   
        assignin('base', 'TX', TX)
        Control(control_i).Command = 'update&Run';
        Control(control_i).Parameters = {'TX'};
        control_i=control_i+1;
    end
elseif current_startevent == event_nr_live % Start Premeasurement
    Resource.Parameters.startEvent = 1;
else
    error("Something must have gone wrong with assigning startEvent" )
end
assignin('base', 'Resource', Resource);
Control = evalin('base', 'Control');
Control(control_i).Command = 'update&Run';
Control(control_i).Parameters = {'Parameters'};
assignin('base','Control', Control);
%CalcImageLevelCallback

%SaveImageDataPreM
event_nr_PreM_saving = evalin('base','event_nr_PreM_saving');
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Parameters',1,'startEvent',event_nr_PreM_saving};
assignin('base','Control',Control);
%SaveImageDataPreM