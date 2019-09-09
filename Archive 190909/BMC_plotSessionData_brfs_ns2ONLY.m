%BMC_plotSessionData_brfs_ns2ONLY
%GOAL: plot LFP,CSD,and fCSD,for each session
%   Version 2.0
%   Brock Carlson -- created 9/5/19
%   Taken from BMC_plotSessionDat_rforia.m but commented out all ns6 file stuff
%   and allowed for new sortBRFSStimandTimeData.m function
%   Also removed all photodiode triggers
%
%   This code marks the end of the August Update branch in github to be
%   merged onto master. This code is heavily based on the BMC_rforiPlot.m,
%   which is the current gold-standard for the BMC analysis path. However,
%   it is broken into functions that work as of 9/9/19.
%
%   Future iterations of this code will simply process the data and save
%   matlab variables. Plotting will be done in a seperate matlab script so
%   that the pre-processing and does not have to be re-run every
%   single time.
%
%  % used as basis for BMC_processBrfsDat.m and BMC_plotBrfs.m
%
%   Moved to archive on 09/09/19

clear

%% EDITABLE VARIABLES
filename = {'160102_E_brfs001'}';
% % filename = {'160102_E_rfori002','160102_E_brfs001','160427_E_rfori002',...
% %     '160427_E_brfs001','160510_E_rfori001','160510_E_brfs001'}';
sinkAllocate = 'BMC_DfS';
pre = 50;
post = 1600;
TM = -pre:1:post;

% variables for end of script
savefiledir = 'G:\LaCie\SfN 2019--figsAndMatVars\SfN 2019 figs\brfs conditions diagnostics';



for a = 1:size(filename,1)
    clearvars -except a filename sinkAllocate pre post TM savefiledir
    
    disp(filename{a})
    
if ispc
    directory = strcat('G:\LaCie\all BRFS\',filename{a}(1:8));
else
    directory = '/Volumes/Drobo/DATA/NEUROPHYS/carlsobm';
end

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
dayID = strfind(SessionParams.Date',str2double(filename{a}(1:6)));
PARAMS = SessionParams(dayID,:);

%Set Channel vector based on the number of electrode contacts
chans = [1:PARAMS.el];

%% make names to load later
% Neural data filenames
cd(directory)
readNEVfile = strcat(directory,filesep,filename{a},'.nev');
readNS2file = strcat(directory,filesep,filename{a},'.ns2');
% % % % readNS6file = strcat(directory,filesep,filename{a},'.ns6');

[grating,readGRATINGfile] = formatAndOpenGratingTxt(directory,filename{a});







%% DIG INTO PREALLOCATED VARIABLES

%% Load event event codes and times
NEV = openNEV(readNEVfile,'noread','overwrite');
EventCodes      = NEV.Data.SerialDigitalIO.UnparsedData - 128;        % we don't know why we have to subtract 128 but we do
EventSamples    = NEV.Data.SerialDigitalIO.TimeStamp;                 % in samples 
EventTimes      = floor(NEV.Data.SerialDigitalIO.TimeStampSec.*1000); % convert to ms 
[pEvC, pEvT]    = parsEventCodesML(EventCodes,EventSamples);          % sorts codes, samps or times into trials
% % % % [pEvT_photo,tf] = pEvtPhoto2(readGRATINGfile,pEvC,pEvT,mode(grating.ypos),[],'ainp1',0); % photo diode signal
%% SORT pEvC/pEvT TO MATCH GRAITN VARIABLE
% So far all of these data are from EVERY trial, including trials where
% animal breaks fixation. Lets get rid of those and make a new structure
% with the grating info and the stimulus onsets 
STIM            = sortStimandTimeData(grating,pEvC,pEvT,'stim'); % this is in nbanalysis. definitely double check it before you use it. 
% % % % STIM_photo      =  sortStimandTimeData(grating,pEvC,pEvT_photo,'stim'); % uses the photodiode inputs

%%%%%%%%%%%% Note, STIM.onsets, is now an index of the column position for
%%%%%%%%%%%% the stimulus onset in the 30kHz data

STIM.onsetsdown         = floor(STIM.onsets./30); % necessary even for ns2 file analysis and use of photo diode. pEvT is in 30kHz sampling time. 
% % % % STIM_photo.onsetsdown   = floor(STIM_photo.onsets./30);

STIM_BRFS = sortBrfsStimandTimeData(grating,pEvC,pEvT,PARAMS);


%% LOAD NEURAL DATA

%% Load pin-by-pin
% Load pin-by-pin and order
%'noread' for Header
NS2_Header      = openNSx(readNS2file,'noread');
% % % % NS6_Header      = openNSx(readNS6file,'noread');
% create a logical array indexing the position of the contacts for the
% electrode penetrated into V1 for the day
contactLogical = strcmp({NS2_Header.ElectrodesInfo.ConnectorBank},PARAMS.V1bank{1,1}(2));
contactLabels = {NS2_Header.ElectrodesInfo(contactLogical).Label};
contactNum = size(contactLabels,2);
if contactNum ~= PARAMS.el
    error('Houston we have a problem') %I'm not totally sure this set-up will always work... BMC 8/26/19
end

count = 0;
for i = 1:size(contactLogical,2)
    if contactLogical(i) == 1
        count = count+1;
        disp(strcat('i=',num2str(count)))

        %sort electrode contacts for Plexon DfS
        contactName = strcat(PARAMS.V1bank{1,1},sprintf('%02d',i));
        idx = contains(contactLabels,contactName);
        contactPosition = find(idx);
        pinNum = sprintf('c:%u',contactPosition);
        NS2         = openNSx(readNS2file,pinNum,'read','uV');
% % % %         NS6         = openNSx(readNS6file,pinNum,'read','uV');
        %preallocate
            if count == 1
               sampleNumNS2 = length(NS2.Data); 
               ns2DAT_predivide = zeros(contactNum,sampleNumNS2);
% % % %                sampleNumNS6 = length(NS6.Data);
% % % %                ns6DAT_predivide = zeros(contactNum,sampleNumNS6);
% % % %                elLabelsOut = strings(contactNum,1);
            end
        ns2DAT_predivide(count,:) = NS2.Data;
% % % %         ns6DAT_predivide(count,:) = NS6.Data;
        elLabelsOut(count,:) = contactName;
        clear pinNum NS2 NS6 contactName idx
    else
        continue
    end
end
ns2DAT = ns2DAT_predivide./4;%convert units to  uV
% % % % ns6DAT = ns6DAT_predivide./4;
%flip if NN,  most superficial channel on top, regardless of number
if strcmp(string(PARAMS.SortDirection), 'descending')
    disp('flipped for NN array')
    ns2DAT = flipud(ns2DAT);
% % % %     ns6DAT = flipud(ns6DAT);
    elLabelsOut = flipud(elLabelsOut);
end
 





    
%% FILTER AND DOWNSAMPLE

%% filter  LFP and downsample the ns6LFP immediatly to align with STIM.onsets in 1kHz
    fc = 50;
    fs = 1000;
    [butter_b,butter_a] = butter(4,fc/(fs/2));
    ns2LFP = nan(size(ns2DAT));
% % % %     ns6LFP = nan(size(ns6DAT));
    for j = 1:contactNum
        disp(strcat('j=',num2str(j)))
        clear ns2FiltVec ns6FiltVec
        ns2FiltVec = ns2DAT(j,:)';
% % % %         ns6FiltVecLFP = ns6DAT(j,:)';
        ns2LFP(j,:) = filtfilt(butter_b,butter_a,ns2FiltVec);
% % % %         ns6LFP(j,:) = filtfilt(butter_b,butter_a,ns6FiltVecLFP);
    end
    

    
    
    
    
%% CALCULATE CSD

ns2CSD   = padarray(calcCSD(ns2LFP'),[1 0],NaN,'replicate'); 

%% TRIGGER TO pEvC/pEvT at 1 kHz 
% Structure of section.
    % A) Trigger to STIM.onsetsdown
        % i) trigger ns2
                % output is TRIG.ns2LFP T TRIG.ns2CSD



        
% A) Trigger to STIM.onsetsdown 
% remove TPs that are too close to start or end
STIM.onsetsdown(STIM.onsetsdown - pre  < 0) = [];
STIM.onsetsdown(STIM.onsetsdown + post > size(ns2LFP,2)) = [];

% get dimensions for preallocation
contactLogical = strcmp({NS2_Header.ElectrodesInfo.ConnectorBank},PARAMS.V1bank{1,1}(2));
contactLabels = {NS2_Header.ElectrodesInfo(contactLogical).Label};
contactNum = size(contactLabels,2);
if contactNum ~= PARAMS.el
    error('Houston we have a problem') %I'm not totally sure this set-up will always work... BMC 8/26/19
end

numTriggers.STIM = size(STIM.onsetsdown,1);

%preallocate
TRIG.ns2LFP         = NaN( contactNum,(post + pre + 1),numTriggers.STIM); 
TRIG.ns2CSD         = NaN( contactNum,(post + pre + 1),numTriggers.STIM);


    
%%% DATA MUST BE DOWNSAMPLED TO 1KH IN ORDER TO TRIGGER PROPERLY
for singleCh = 1:contactNum 
    for singleTrigger = 1:numTriggers.STIM
        timeOfTrigger = STIM.onsetsdown(singleTrigger);
        windowOfTrigger = timeOfTrigger-pre:timeOfTrigger+post;
        % output is (Ch x time x triggerNumber)
        TRIG.ns2LFP(singleCh,:,singleTrigger)       = ns2LFP(singleCh,windowOfTrigger); 
        TRIG.ns2CSD(singleCh,:,singleTrigger)       = ns2CSD(singleCh,windowOfTrigger);

    end
end

%%%%
%%%%%%% TRIGGER BRFS 
%%%%
STIM_BRFS.diop_simult_PS.start_noSoaDown(STIM_BRFS.diop_simult_PS.start_noSoaDown - pre  < 0) = [];
STIM_BRFS.diop_simult_PS.start_noSoaDown(STIM_BRFS.diop_simult_PS.start_noSoaDown + post > size(ns2LFP,2)) = [];

numTriggers.diop_simult_PS = size(STIM_BRFS.diop_simult_PS.start_noSoaDown,1);

%preallocate
TRIG_BRFS.diop_simult_PS.ns2LFP         = NaN( contactNum,(post + pre + 1),numTriggers.diop_simult_PS); 
TRIG_BRFS.diop_simult_PS.ns2CSD         = NaN( contactNum,(post + pre + 1),numTriggers.diop_simult_PS);


    
%%% DATA MUST BE DOWNSAMPLED TO 1KH IN ORDER TO TRIGGER PROPERLY
for singleCh = 1:contactNum 
    for singleTrigger = 1:numTriggers.diop_simult_PS
        timeOfTrigger = STIM_BRFS.diop_simult_PS.start_noSoaDown(singleTrigger);
        windowOfTrigger = timeOfTrigger-pre:timeOfTrigger+post;
        % output is (Ch x time x triggerNumber)
        TRIG_BRFS.diop_simult_PS.ns2LFP(singleCh,:,singleTrigger)       = ns2LFP(singleCh,windowOfTrigger); 
        TRIG_BRFS.diop_simult_PS.ns2CSD(singleCh,:,singleTrigger)       = ns2CSD(singleCh,windowOfTrigger);

    end
end


%% AVERAGE AND BASELINE-CORRECT TRIGGERED DATA

%% Average
% Avgerage TRIG
fields.TRIG = fieldnames(TRIG);
for avtr=1:numel(fields.TRIG)
    AVG.(fields.TRIG{avtr})  = mean(TRIG.(fields.TRIG{avtr}),3);
end

firstfields.TRIG_BRFS = fieldnames(TRIG_BRFS);
for ffav = 1:numel(firstfields.TRIG_BRFS)
    subfields.TRIG_BRFS = fieldnames(TRIG_BRFS.(firstfields.TRIG_BRFS{ffav}));
    for avtr=1:numel(subfields.TRIG_BRFS)
        AVG_BRFS.(firstfields.TRIG_BRFS{ffav}).(subfields.TRIG_BRFS{avtr})  = mean(TRIG_BRFS.(firstfields.TRIG_BRFS{ffav}).(subfields.TRIG_BRFS{avtr}),3);
    end
end

%% Baseline corect
% bl AVG
fields.AVG = fieldnames(AVG);
for blavg=1:numel(fields.AVG)
    bl  = mean(AVG.(fields.AVG{blavg})(:,TM<0),2);
    BLavg.(fields.AVG{blavg})  = AVG.(fields.AVG{blavg}) - bl;
end

firstfields.AVG_BRFS = fieldnames(AVG_BRFS);
for ffbl = 1:numel(firstfields.AVG_BRFS)
    subfields.AVG_BRFS = fieldnames(AVG_BRFS.(firstfields.AVG_BRFS{ffbl}));
    for blavg=1:numel(subfields.AVG_BRFS)
        bl  = mean(AVG_BRFS.(firstfields.AVG_BRFS{ffbl}).(subfields.AVG_BRFS{blavg})(:,TM<0),2);
        BLavg_BRFS.(firstfields.AVG_BRFS{ffbl}).(subfields.AVG_BRFS{blavg})  = AVG_BRFS.(firstfields.AVG_BRFS{ffbl}).(subfields.AVG_BRFS{blavg}) - bl;
    end
end
%% Filter and interpolate CSD, ns2 and ns5
% AVG
AVG.ns2fCSD = filterCSD(AVG.ns2CSD);
AVG_BRFS.diop_simult_PS.ns2fCSD = filterCSD(AVG_BRFS.diop_simult_PS.ns2CSD);

% BLavg
BLavg.ns2fCSD       = filterCSD(BLavg.ns2CSD);
BLavg_BRFS.diop_simult_PS.ns2fCSD = filterCSD(BLavg_BRFS.diop_simult_PS.ns2CSD);



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
subplot(1,3,1)
f_ShadedLinePlotbyDepth(BLavg_BRFS.diop_simult_PS.ns2LFP,chans,TM,[],1)
plot([0 0], ylim,'k')
% plot([800 800], ylim,'k')
title({'LFP',filename{a}}, 'Interpreter', 'none')
xlabel('time (ms)')
ylabel('Contacts indexed down from surface')

% CSD line
subplot(1,3,2)
f_ShadedLinePlotbyDepth(BLavg_BRFS.diop_simult_PS.ns2CSD,chans,TM,[],1)
plot([0 0], ylim,'k')
% plot([800 800], ylim,'k')
title({'CSD',filename{a}}, 'Interpreter', 'none')
xlabel('time (ms)')

% fCSD
subplot(1,3,3)
imagesc(TM,chans,BLavg_BRFS.diop_simult_PS.ns2fCSD); 
colormap(flipud(jet));
climit = max(abs(get(gca,'CLim'))*.8);
set(gca,'CLim',[-climit climit],'Box','off','TickDir','out')
hold on;
plot([0 0], ylim,'k')
% plot([800 800], ylim,'k')
clrbar = colorbar;
title({'interpolate CSD',filename{a}}, 'Interpreter', 'none')
xlabel('time (ms)')
clrbar.Label.String = 'nA/mm^3';

cd(savefiledir)

%% SAVE
% % % cd(savefiledir)
% % % saveas(gcf, filename{a}, 'fig')
% % % saveas(gcf, filename{a}, 'pdf')
% % % saveas(gcf, filename{a}, 'png')

%% SCRIPT END
toc

end

load gong.mat;
sound(y);





