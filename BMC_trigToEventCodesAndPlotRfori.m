%BMC_trigToEventCodesAndPlotRfori.m
%GOAL: load the previously filtered and saved brfs data, trigger, and plot
%
%   Version 1.1
%   Brock Carlson -- created 9/11/19
%   
%   Current plotting goal: plot brfs monocular fCSDs for each ori and eye
%   under scaled conditions to see if an effect exists.



clear

%% EDITABLE VARIABLES
% % filename = {'160427_E_rfori002'}';
filename = {'160102_E_rfori002','160427_E_rfori002','160510_E_rfori002'}';
sinkAllocate = 'BMC_DfS';
pre = 50;
post = 250;
TM = -pre:1:post;
nameSaveType = 'LFPandCSDof';
figtype = 'rfori_ori1vsori2';


% Computer-specific editable variables 
if strcmp(getenv('USER'),'maierav')
    % @Alex -- fill in necessary information for your system here.
        % savefiledir = 'G:\LaCie\SfN 2019--figsAndMatVars\SfN 2019 figs\brfs conditions diagnostics';        
elseif strcmp(getenv('USERNAME'),'Brock Carlson')
    % variables for end of script
    savefiledir = 'G:\LaCie\SfN 2019--figsAndMatVars\SfN 2019 figs\Is there an effect in monocular CSD';
elseif ~ispc && strcmp(getenv('USERNAME'),'Brock')
    savefiledir = 'INSERT PATH HERE';
end



for a = 1:size(filename,1)
    clearvars -except a filename sinkAllocate pre post TM savefiledir nameSaveType figtype
    
    disp(filename{a})
    

%% Computer-specific directories 
if strcmp(getenv('USER'),'maierav')
    % @Alex -- fill in necessary information for your system here.
        % %     addpath(genpath('/Users/alex 1/Desktop/LAB/Brock'));
        % %     drname        = {'/Users/alex 1/Desktop/LAB/Brock/Data'};
        % load session params
    %session params
    sessionParamDir = 'G:/LaCie';
elseif strcmp(getenv('USERNAME'),'Brock Carlson')
    addpath(genpath('G:\LaCie\all BRFS'));
    dataDirectory = strcat('G:\LaCie\all BRFS\',filename{a}(1:8));
    %session params
    sessionParamDir = 'G:/LaCie';

elseif ~ispc && strcmp(getenv('USERNAME'),'Brock')
    dataDirectory = strcat('/Volumes/Drobo/DATA/NEUROPHYS/carlsobm/',filename{a}(1:8));
    %session params
    sessionParamDir = '/Volumes/Drobo/DATA/NEUROPHYS/carlsobm/';
end




%% LOAD IN SESSION DATA
tic

cd(sessionParamDir)
load('SessionParams.mat')
    
switch sinkAllocate
    case 'BMC_DfS'
        SessionParams.EvalSink = SessionParams.BMC_DfS;
    case 'Old_DfS'
        SessionParams.EvalSink = SessionParams.Old_DfS;
end

% get info needed for day I'm analyzing
dayID = strfind(SessionParams.Date',str2double(filename{a}(1:6)));
PARAMS = SessionParams(dayID,:);
if contains(filename{a},'160427_E_rfori002')
   PARAMS.PS = 125;
   PARAMS.NPS = 35;
end

%Set Channel vector based on the number of electrode contacts
chans = 1:PARAMS.el;

%% LOAD FILTERED CONTINUOUS DATA
% The saved .mat variables should contain both LFP and CSD
cd(savefiledir)
loadname = strcat(nameSaveType,filename{a},'.mat');
load(loadname)






%% DIG INTO PREALLOCATED VARIABLES
%
%% Load grating
[grating,readGRATINGfile] = formatAndOpenGratingTxt(dataDirectory,filename{a});

%% Load event event codes and times
readNEVfile = strcat(dataDirectory,filesep,filename{a},'.nev');
NEV = openNEV(readNEVfile,'noread','overwrite');
EventCodes      = NEV.Data.SerialDigitalIO.UnparsedData - 128;        % we don't know why we have to subtract 128 but we do
EventSamples    = NEV.Data.SerialDigitalIO.TimeStamp;                 % in samples 
EventTimes      = floor(NEV.Data.SerialDigitalIO.TimeStampSec.*1000); % convert to ms 
[pEvC, pEvT]    = parsEventCodesML(EventCodes,EventSamples);          % sorts codes, samps or times into trials
%% SORT pEvC/pEvT TO MATCH GRAITNG VARIABLE
% So far all of these data are from EVERY trial, including trials where
% animal breaks fixation. Lets get rid of those and make a new structure
% with the grating info and the stimulus onsets. 
    % NOTE: STIM.onsets, is now an index of the column position for the...
    % stimulus onset in the 30kHz data. This is why, for using it with the ns2
    % file, we need to downsample to STIM.onsetsdown.
STIM            = sortStimandTimeData(grating,pEvC,pEvT,'stim'); 
STIM.onsetsdown         = floor(STIM.onsets./30); % necessary even for ns2 file analysis and use of photo diode. pEvT is in 30kHz sampling time. 

% Get RFORI trials where the two orientations shown in brfs are displayed
% to the dominant eye
clear i count1 count2
count1 = 0; count2 = 0;
for i = 1:size(STIM.tilt,1)
    if STIM.tilt(i) == PARAMS.PS
        count1 = count1+1;
        STIM.ori1onsets(count1,:) = STIM.onsetsdown(i);
    elseif STIM.tilt(i) == PARAMS.NPS
        count2 = count2+1;
        STIM.ori2onsets(count2,:) = STIM.onsetsdown(i);
    end    
end

% get dimensions for preallocation
readNS2file = strcat(dataDirectory,filesep,filename{a},'.ns2');
NS2_Header      = openNSx(readNS2file,'noread');
contactLogical = strcmp({NS2_Header.ElectrodesInfo.ConnectorBank},PARAMS.V1bank{1,1}(2));
contactLabels = {NS2_Header.ElectrodesInfo(contactLogical).Label};
contactNum = size(contactLabels,2);
if contactNum ~= PARAMS.el
    error('Houston we have a problem') %I'm not totally sure this set-up will always work... BMC 8/26/19
end

%% TRIGGER SECTION
%%%
%%%
%%%
%% trig all onsets
numTriggers.STIM.onsetsdown = size(STIM.onsetsdown,1);
%preallocate
    TRIG.ns2LFP         = NaN( contactNum,(post + pre + 1),numTriggers.STIM.onsetsdown); 
    TRIG.ns2CSD         = NaN( contactNum,(post + pre + 1),numTriggers.STIM.onsetsdown);  
%%% DATA MUST BE DOWNSAMPLED TO 1KH IN ORDER TO TRIGGER PROPERLY
for singleCh = 1:contactNum 
    for singleTrigger = 1:numTriggers.STIM.onsetsdown
        timeOfTrigger = STIM.onsetsdown(singleTrigger);
        windowOfTrigger = timeOfTrigger-pre:timeOfTrigger+post;
        % output is (Ch x time x triggerNumber)
        TRIG.ns2LFP(singleCh,:,singleTrigger)       = ns2LFP(singleCh,windowOfTrigger); 
        TRIG.ns2CSD(singleCh,:,singleTrigger)       = ns2CSD(singleCh,windowOfTrigger);

    end
end


%% trig ori1 and ori2 onsets
numTriggers.STIM.ori1onsets = size(STIM.ori1onsets,1);
numTriggers.STIM.ori2onsets = size(STIM.ori2onsets,1);
%preallocate
TRIG.ori1LFP         = NaN( contactNum,(post + pre + 1),numTriggers.STIM.ori1onsets); 
TRIG.ori1CSD         = NaN( contactNum,(post + pre + 1),numTriggers.STIM.ori1onsets); 
TRIG.ori2LFP         = NaN( contactNum,(post + pre + 1),numTriggers.STIM.ori2onsets); 
TRIG.ori2CSD         = NaN( contactNum,(post + pre + 1),numTriggers.STIM.ori2onsets);   
%%% DATA MUST BE DOWNSAMPLED TO 1KH IN ORDER TO TRIGGER PROPERLY
for singleCh = 1:contactNum 
    for singleTrigger = 1:numTriggers.STIM.ori1onsets
        timeOfTrigger = STIM.ori1onsets(singleTrigger);
        windowOfTrigger = timeOfTrigger-pre:timeOfTrigger+post;
        % output is (Ch x time x triggerNumber)
        TRIG.ori1LFP(singleCh,:,singleTrigger)       = ns2LFP(singleCh,windowOfTrigger); 
        TRIG.ori1CSD(singleCh,:,singleTrigger)       = ns2CSD(singleCh,windowOfTrigger);

    end
end
for singleCh = 1:contactNum 
    for singleTrigger = 1:numTriggers.STIM.ori2onsets
        timeOfTrigger = STIM.ori2onsets(singleTrigger);
        windowOfTrigger = timeOfTrigger-pre:timeOfTrigger+post;
        % output is (Ch x time x triggerNumber)
        TRIG.ori2LFP(singleCh,:,singleTrigger)       = ns2LFP(singleCh,windowOfTrigger); 
        TRIG.ori2CSD(singleCh,:,singleTrigger)       = ns2CSD(singleCh,windowOfTrigger);

    end
end



%% AVERAGE AND BASELINE-CORRECT TRIGGERED DATA 

%% Average
% Avgerage TRIG
fields.TRIG = fieldnames(TRIG);
for avtr=1:numel(fields.TRIG)
    AVG.(fields.TRIG{avtr})  = mean(TRIG.(fields.TRIG{avtr}),3);
end


%% Baseline corect
% bl AVG
fields.AVG = fieldnames(AVG);
for blavg=1:numel(fields.AVG)
    bl  = mean(AVG.(fields.AVG{blavg})(:,TM<0),2);
    BLavg.(fields.AVG{blavg})  = AVG.(fields.AVG{blavg}) - bl;
end


%% Filter and interpolate CSD
% filter AVG struct
AVG.ns2fCSD_allonsets   = filterCSD(AVG.ns2CSD);
AVG.ns2fCSD_ori1        = filterCSD(AVG.ori1CSD);
AVG.ns2fCSD_ori2        = filterCSD(AVG.ori2CSD);

% filter BLavg struct
BLavg.ns2fCSD_allonsets   = filterCSD(BLavg.ns2CSD);
BLavg.ns2fCSD_ori1        = filterCSD(BLavg.ori1CSD);
BLavg.ns2fCSD_ori2        = filterCSD(BLavg.ori2CSD);

%% POSTPROCESS
% downsample and get PSD

%Downsample LFP so I can easily calculate the psd

%PSD

%% Align Across Days
    %
    % create alignment matrices if necessary at a later date
    %
    
%% Plot
climit = {600,1000,1000};

figure;
set(gcf, 'Position', [680 338 711 760]);
subplot(1,2,1)
imagesc(TM,chans,BLavg.ns2fCSD_ori1); 
colormap(flipud(jet));
set(gca,'CLim',[-climit{a} climit{a}],'Box','off','TickDir','out')
% % climit = max(abs(get(gca,'CLim'))*.8);
% % set(gca,'CLim',[-climit climit],'Box','off','TickDir','out')
hold on;
plot([0 0], ylim,'k')
plot([800 800], ylim,'k')
clrbar = colorbar;
title({'rfori_ori1',filename{a}}, 'Interpreter', 'none')
xlabel('time (ms)')
clrbar.Label.String = 'nA/mm^3';


subplot(1,2,2)
imagesc(TM,chans,BLavg.ns2fCSD_ori2); 
colormap(flipud(jet));
set(gca,'CLim',[-climit{a} climit{a}],'Box','off','TickDir','out')
% % climit = max(abs(get(gca,'CLim'))*.8);
% % set(gca,'CLim',[-climit climit],'Box','off','TickDir','out')
hold on;
plot([0 0], ylim,'k')
plot([800 800], ylim,'k')
clrbar = colorbar;
title({'rfori_ori2',filename{a}}, 'Interpreter', 'none')
xlabel('time (ms)')
clrbar.Label.String = 'nA/mm^3';

set(gcf, 'Position', [75 501 517 569]);


%% SAVE figs
cd(savefiledir)
figsavename = strcat(figtype,'_',filename{a});
saveas(gcf, figsavename, 'fig')
saveas(gcf, figsavename, 'pdf')
saveas(gcf, figsavename, 'png')
    



end




% % % % %% PLOT 1,2,3,4 "Session Data diagnostic"
% % % % % LFP, CSDshadedLine, fCSD, aMUA, 
% % % % figure
% % % % %LFPline
% % % % subplot(1,3,1)
% % % % f_ShadedLinePlotbyDepth(BLavg_BRFS.diop_simult_NPS.ns2LFP,chans,TM,[],1)
% % % % plot([0 0], ylim,'k')
% % % % plot([800 800], ylim,'k')
% % % % title({'LFP',filename{a}}, 'Interpreter', 'none')
% % % % xlabel('time (ms)')
% % % % ylabel('Contacts indexed down from surface')
% % % % 
% % % % % CSD line
% % % % subplot(1,3,2)
% % % % f_ShadedLinePlotbyDepth(BLavg_BRFS.diop_simult_NPS.ns2CSD,chans,TM,[],1)
% % % % plot([0 0], ylim,'k')
% % % % plot([800 800], ylim,'k')
% % % % title({'CSD',filename{a}}, 'Interpreter', 'none')
% % % % xlabel('time (ms)')
% % % % 
% % % % % fCSD
% % % % subplot(1,3,3)
% % % % imagesc(TM,chans,BLavg_BRFS.diop_simult_NPS.ns2fCSD); 
% % % % colormap(flipud(jet));
% % % % climit = max(abs(get(gca,'CLim'))*.8);
% % % % set(gca,'CLim',[-climit climit],'Box','off','TickDir','out')
% % % % hold on;
% % % % plot([0 0], ylim,'k')
% % % % plot([800 800], ylim,'k')
% % % % clrbar = colorbar;
% % % % title({'interpolate CSD',filename{a}}, 'Interpreter', 'none')
% % % % xlabel('time (ms)')
% % % % clrbar.Label.String = 'nA/mm^3';
% % % % 





