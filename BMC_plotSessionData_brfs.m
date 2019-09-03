%BMC_plotSessionData
%GOAL: plot LFP,CSD,fCSD, and aMUA for each session
%   Version 1.0
%   Brock Carlson -- created 8/27/19
%   Taken from BMC_rforiPlot.m but broken into smaller functions because
%   160510 rfori001 was crashing the system.


clear

%% EDITABLE VARIABLES
filename = '160510_E_brfs001';
if ispc
    directory = 'G:\LaCie\all BRFS\160510_E';
else
    directory = '/Volumes/Drobo/DATA/NEUROPHYS/carlsobm';
end
sinkAllocate = 'BMC_DfS';
pre = 50;
post = 250;
TM = -pre:1:post;

%% LOAD IN SESSION DATA
tic
%% Load session params
if ispc
    cd('G:\LaCie')
else
    cd(directory)
end
load('SessionParams.mat')
switch sinkAllocate
    case 'BMC_DfS'
        SessionParams.EvalSink = SessionParams.BMC_DfS;
    case 'Old_DfS'
        SessionParams.EvalSink = SessionParams.Old_DfS;
end

% get info needed for day I'm analyzing
dayID = strfind(SessionParams.Date',str2double(filename(1:6)));
PARAMS = SessionParams(dayID,:);

%Set Channel vector based on the number of electrode contacts
chans = [1:PARAMS.el];

%% make names to load later
% Neural data filenames
cd(directory)
readNEVfile = strcat(directory,filesep,filename,'.nev');
readNS2file = strcat(directory,filesep,filename,'.ns2');
readNS6file = strcat(directory,filesep,filename,'.ns6');

%% load grating from text file
[grating,readGRATINGfile] = formatAndOpenGratingTxt(directory,filename);

%% DIG INTO PREALLOCATED VARIABLES

%% Load event event codes and times
NEV = openNEV(readNEVfile,'noread','overwrite');
EventCodes      = NEV.Data.SerialDigitalIO.UnparsedData - 128;        % we don't know why we have to subtract 128 but we do
EventSamples    = NEV.Data.SerialDigitalIO.TimeStamp;                 % in samples 
EventTimes      = floor(NEV.Data.SerialDigitalIO.TimeStampSec.*1000); % convert to ms 
[pEvC, pEvT]    = parsEventCodesML(EventCodes,EventSamples);          % sorts codes, samps or times into trials
% [pEvT_photo,tf] = pEvtPhoto2_BMCbrfs(readGRATINGfile,pEvC,pEvT,mode(grating.ypos),[],'ainp1',0); % photo diode signal
%% SORT pEvC/pEvT TO MATCH GRAITN VARIABLE
% So far all of these data are from EVERY trial, including trials where
% animal breaks fixation. Lets get rid of those and make a new structure
% with the grating info and the stimulus onsets 
STIM            = sortStimandTimeData_bmc(grating,pEvC,pEvT,'stim'); % this is in nbanalysis. definitely double check it before you use it. 
% STIM_photo      =  sortStimandTimeData(grating,pEvC,pEvT_photo,'stim'); % uses the photodiode inputs

%%%%%%%%%%%% Note, STIM.onsets, is now an index of the column position for
%%%%%%%%%%%% the stimulus onset in the 30kHz data

STIM.onsetsdown         = floor(STIM.onsets./30);
% STIM_photo.onsetsdown   = floor(STIM_photo.onsets./30);

%% LOAD NEURAL DATA
% pin-by-pin function
 [ns2DAT,ns6DAT,elLabelsOut] = loadNeuralData(readNS2file,readNS6file,PARAMS);
 
%% FILTER AND DOWNSAMPLE
% filter  LFP and downsample the ns6LFP immediatly to align with STIM.onsets in 1kHz
    [ns2LFP,ns6LFPdown] = filterForLFP(ns2DAT,ns6DAT,readNS2file,readNS6file,PARAMS);
% filter aMUA and downsample immediatly to align with STIM.onsets in 1kHz
    [aMUAdown] = filterForaMUA(ns6DAT,readNS2file,readNS6file,PARAMS);

%% CLEAR UNNECESSARY VARIABLES
clear ns2DAT ns6DAT 
clear readGRATINGfile readNEVfile 
clear EventCodes EventSamples Event Times
%% CALCULATE CSD.
% NOTE: calcCSD operates along dimensions of (continuous data x channels).
ns2CSD   = padarray(calcCSD(ns2LFP'),[1 0],NaN,'replicate'); 
ns6CSD   = padarray(calcCSD(ns6LFPdown'),[1 0],NaN,'replicate'); 

%% Trigger data
[TRIG] = bmcTRIG(STIM,PARAMS,pre,post,readNS2file, ...
    readNS6file,ns2LFP,ns2CSD,ns6LFPdown,ns6CSD,aMUAdown);
% % % [TRIG_photo] = bmcTRIG_photo(STIM_photo,PARAMS,pre,post,readNS2file, ...
% % %     readNS6file,ns2LFP,ns2CSD,ns6LFPdown,ns6CSD,aMUAdown);

%% AVERAGE AND BASELINE-CORRECT TRIGGERED DATA

%% Average
% Avgerage TRIG
fields.TRIG = fieldnames(TRIG);
for avtr=1:numel(fields.TRIG)
    AVG.(fields.TRIG{avtr})  = mean(TRIG.(fields.TRIG{avtr}),3);
end


% % % % Avgerage TRIG_photo
% % % fields.TRIG_photo = fieldnames(TRIG_photo);
% % % for avtrp=1:numel(fields.TRIG_photo)
% % %     AVG_photo.(fields.TRIG_photo{avtrp})  = mean(TRIG_photo.(fields.TRIG_photo{avtrp}),3);
% % % end
%% Baseline corect
% bl AVG
fields.AVG = fieldnames(AVG);
for blavg=1:numel(fields.AVG)
    bl  = mean(AVG.(fields.AVG{blavg})(:,TM<0),2);
    BLavg.(fields.AVG{blavg})  = AVG.(fields.AVG{blavg}) - bl;
end

% % % % % bl AVG_photo
% % % % fields.AVG_photo = fieldnames(AVG_photo);
% % % % for blavgp=1:numel(fields.AVG_photo)
% % % %     bl_photo  = mean(AVG_photo.(fields.AVG_photo{blavgp})(:,TM<0),2);
% % % %     BLavg_photo.(fields.AVG_photo{blavgp})  = AVG_photo.(fields.AVG_photo{blavgp}) - bl_photo;
% % % % end


%% Filter and interpolate CSD, ns2 and ns5
% AVG
AVG.ns2fCSD = filterCSD(AVG.ns2CSD);
AVG.ns6fCSD = filterCSD(AVG.ns6CSD);

% % % % AVG_photo
% % % AVG_photo.ns2fCSD = filterCSD(AVG_photo.ns2CSD);
% % % AVG_photo.ns6fCSD = filterCSD(AVG_photo.ns6CSD);

% BLavg
BLavg.ns2fCSD       = filterCSD(BLavg.ns2CSD);
BLavg.ns6fCSD       = filterCSD(BLavg.ns6CSD);

% % % % BLavg_photo
% % % BLavg_photo.ns2fCSD       = filterCSD(BLavg_photo.ns2CSD);
% % % BLavg_photo.ns6fCSD       = filterCSD(BLavg_photo.ns6CSD);

%% POSTPROCESS
% downsample and get PSD

%Downsample LFP so I can easily calculate the psd

%PSD

%% Align Across Days
    %
    % create alignment matrices if necessary at a later date
    %
    
%% PLOT 1,2,3,4,5
% LFP, CSDshadedLine, fCSD, aMUA, PSD
figure
%LFPline
subplot(1,4,1)
f_ShadedLinePlotbyDepth(BLavg.ns6LFP,chans,TM,[],1)
plot([0 0], ylim,'k')
title({'LFP',filename}, 'Interpreter', 'none')
xlabel('time (ms)')
ylabel('Contacts indexed down from surface')

% CSD line
subplot(1,4,2)
f_ShadedLinePlotbyDepth(BLavg.ns6CSD,chans,TM,[],1)
plot([0 0], ylim,'k')
title({'CSD',filename}, 'Interpreter', 'none')
xlabel('time (ms)')

% fCSD
subplot(1,4,3)
imagesc(TM,chans,BLavg.ns6fCSD); 
colormap(flipud(jet));
climit = max(abs(get(gca,'CLim'))*.8);
set(gca,'CLim',[-climit climit],'Box','off','TickDir','out')
hold on;
plot([0 0], ylim,'k')
clrbar = colorbar;
title({'interpolate CSD',filename}, 'Interpreter', 'none')
xlabel('time (ms)')
clrbar.Label.String = 'nA/mm^3';


% aMUA line
subplot(1,4,4)
f_ShadedLinePlotbyDepth(BLavg.aMUA,chans,TM,[],1)
plot([0 0], ylim,'k')
title({'aMUA',filename}, 'Interpreter', 'none')
xlabel('time (ms)')

%% SCRIPT END
toc

load gong.mat;
sound(y);






