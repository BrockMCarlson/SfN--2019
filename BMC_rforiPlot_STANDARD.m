% BMC_rforiPlot.m
%   This is the gold-standard for the BMC analysis path. It was approved at
%   the 9/2/2019 lab meeting. This code does not work for large files as it
%   does not use space effeciently. However, all further analysis should be
%   based on this code because it proper plots LFP, CSD, fCSD, and aMUA
%   with correct indexing, using both ns6 files (with capability of using
%   the ns2 files if necessary), and all triggered to the photodiode.
%   Future iterations of this code broke it down into functions but this
%   code is ideal because everything is written out line-by-line. 
%
%   Note: this code simply plots all stimulus onsets. It does not work for
%   specific triggering in BRFS.
%
%   As of 9/9/19 this code runs.
%
%



%% Make CSD of rfori
clear

%% EDITABLE VARIABLES
filename = '160102_E_rfori002';
if ispc
    directory = 'G:\LaCie\all BRFS\160102_E';
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

% grating file name
patterns   = {'rforidrft','rfsfdrft','posdisparitydrft','disparitydrft','cinterocdrft','coneinterocdrft','conedrft', ...
                'colorflicker','bwflicker','rfori','rfsize','cinteroc','color','rfsf','mcosinteroc','dotmapping'}; 
for p = 1:length(patterns)
   pattern      = patterns{p}; 
   if any(strfind(filename,pattern))
       startlog = strfind(filename,pattern); 
       if ~isequal(filename(startlog:end-3),pattern)
            continue
       else
            match    = patterns{p}; 
       end
   end   
end
if isequal(match,'dotmapping')
    gratingext  = '.gDotsXY_di';
else
    gratingext  = ['.g' upper(match) 'Grating_di']; 
end
readGRATINGfile = strcat(directory,filesep,filename,gratingext);


%% Load text file
if contains(gratingext,'DRFT')
      grating     = readgDRFTGrating([filename gratingext]); % from nbanalysis 
elseif contains(gratingext,'Dots')
      grating     = readgDotsXY([filename gratingext]);
else
      grating     = readgGrating([filename gratingext]);
end






%% DIG INTO PREALLOCATED VARIABLES

%% Load event event codes and times
NEV = openNEV(readNEVfile,'noread','overwrite');
EventCodes      = NEV.Data.SerialDigitalIO.UnparsedData - 128;        % we don't know why we have to subtract 128 but we do
EventSamples    = NEV.Data.SerialDigitalIO.TimeStamp;                 % in samples 
EventTimes      = floor(NEV.Data.SerialDigitalIO.TimeStampSec.*1000); % convert to ms 
[pEvC, pEvT]    = parsEventCodesML(EventCodes,EventSamples);          % sorts codes, samps or times into trials
[pEvT_photo,tf] = pEvtPhoto2(readGRATINGfile,pEvC,pEvT,mode(grating.ypos),[],'ainp1',0); % photo diode signal
%% SORT pEvC/pEvT TO MATCH GRAITN VARIABLE
% So far all of these data are from EVERY trial, including trials where
% animal breaks fixation. Lets get rid of those and make a new structure
% with the grating info and the stimulus onsets 
STIM            = sortStimandTimeData(grating,pEvC,pEvT,'stim'); % this is in nbanalysis. definitely double check it before you use it. 
STIM_photo      =  sortStimandTimeData(grating,pEvC,pEvT_photo,'stim'); % uses the photodiode inputs

%%%%%%%%%%%% Note, STIM.onsets, is now an index of the column position for
%%%%%%%%%%%% the stimulus onset in the 30kHz data

STIM.onsetsdown         = floor(STIM.onsets./30);
STIM_photo.onsetsdown   = floor(STIM_photo.onsets./30);



%% LOAD NEURAL DATA

%% Load pin-by-pin
% Load pin-by-pin and order
%'noread' for Header
NS2_Header      = openNSx(readNS2file,'noread');
NS6_Header      = openNSx(readNS6file,'noread');
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
        NS6         = openNSx(readNS6file,pinNum,'read','uV');
        %preallocate
            if count == 1
               sampleNumNS2 = length(NS2.Data); 
               ns2DAT_predivide = zeros(contactNum,sampleNumNS2);
               sampleNumNS6 = length(NS6.Data);
               ns6DAT_predivide = zeros(contactNum,sampleNumNS6);
               elLabelsOut = strings(contactNum,1);
            end
        ns2DAT_predivide(count,:) = NS2.Data;
        ns6DAT_predivide(count,:) = NS6.Data;
        elLabelsOut(count,:) = contactName;
        clear pinNum NS2 NS6 contactName idx
    else
        continue
    end
end
ns2DAT = ns2DAT_predivide./4;%convert units to  uV
ns6DAT = ns6DAT_predivide./4;
%flip if NN,  most superficial channel on top, regardless of number
if strcmp(string(PARAMS.SortDirection), 'descending')
    disp('flipped for NN array')
    ns2DAT = flipud(ns2DAT);
    ns6DAT = flipud(ns6DAT);
    elLabelsOut = flipud(elLabelsOut);
end
 





    
%% FILTER AND DOWNSAMPLE

%% filter  LFP and downsample the ns6LFP immediatly to align with STIM.onsets in 1kHz
    fc = 50;
    fs = 1000;
    [butter_b,butter_a] = butter(4,fc/(fs/2));
    ns2LFP = nan(size(ns2DAT));
    ns6LFP = nan(size(ns6DAT));
    for j = 1:contactNum
        disp(strcat('j=',num2str(j)))
        clear ns2FiltVec ns6FiltVec
        ns2FiltVec = ns2DAT(j,:)';
        ns6FiltVecLFP = ns6DAT(j,:)';
        ns2LFP(j,:) = filtfilt(butter_b,butter_a,ns2FiltVec);
        ns6LFP(j,:) = filtfilt(butter_b,butter_a,ns6FiltVecLFP);
    end
    
    %Downsample ns6LFP to 1kHz to align with triggerpoints below
    ns6LFPtrans = ns6LFP';
    ns6LFPdownTrans = downsample(ns6LFPtrans,30);
    ns6LFPdown = ns6LFPdownTrans';
    
%% filter aMUA and downsample immediatly to align with STIM.onsets in 1kHz
    lpc1 = 500;
    hpc  = 5000;
    nyq  = NS6_Header.MetaTags.SamplingFreq/2;
    lpc2 = lpc1 / 2;
    hWn = hpc / nyq;
    [ bwb1, bwa1 ] = butter( 4, hWn, 'high' );
    lWn1 = lpc1 / nyq;
    [ bwb2, bwa2 ] = butter( 4, lWn1, 'low' );
    lWn2 = lpc2 / nyq;
    [ bwb3, bwa3 ] = butter( 4, lWn2, 'low' );
    aMUA = zeros(size(ns6DAT));
    for k = 1: contactNum
        disp(strcat('k=',num2str(k)))
        clear ns6FiltVecaMUA hpMUA lpMUA
        ns6FiltVecaMUA = ns6DAT(k,:)';
        hpMUA = filtfilt(bwb1,bwa1,ns6FiltVecaMUA);       
        lpMUA = abs( filtfilt( bwb2, bwa2, hpMUA ) ); 
        aMUA(k,:) = filtfilt( bwb3, bwa3, lpMUA );     
    end
    
    %Downsample aMUA to 1kHz to align with triggerpoints below
    aMUAtrans = aMUA';
    aMUAdownTrans = downsample(aMUAtrans,30);
    aMUAdown = aMUAdownTrans';

    
    
    
    
    
    
%% CALCULATE CSD

%% CSD fro double derivative. Simply used to compare the method. NOT FOR POSTERS OR PAPERS!!!
    %derivitive 
    ns2CSD_diff         = diff(ns2LFP,2,1);
    ns2CSD_diffpad      = padarray(ns2CSD_diff,[1 0],NaN,'replicate');    
    ns6CSD_diff         = diff(ns6LFPdown,2,1);
    ns6CSD_diffpad      = padarray(ns6CSD_diff,[1 0],NaN,'replicate');
    
%% calcCSD_BMC    
    elSpaces = 0.1:0.1:contactNum/10;
    distance = mean(diff(elSpaces));
    for m = 1:(contactNum-2) %creates CSD calculation performance matrix
        for n = 1:contactNum
            if (m == n-1)
                calced(m,n) = -2/distance^2;
            elseif (abs(m-n+1) == 1)
                calced(m,n) = 1/distance^2;
            else
                calced(m,n) = 0;
            end
        end            
    end
    ns2CSD(:,:) = (-1*calced*ns2LFP);
    ns2CSDpad   = padarray(ns2CSD,[1 0],NaN,'replicate');
    ns6CSD(:,:) = (-1*calced*ns6LFPdown);
    ns6CSDpad   = padarray(ns6CSD,[1 0],NaN,'replicate');

%% TRIGGER TO pEvC/pEvT at 1 kHz 
% Structure of section.
    % A) Trigger to STIM.onsetsdown
        % i) trigger ns2
        % i) trigger ns6
                % output is TRIG.ns2LFP TRIG.ns2CSD_diff TRIG.ns2CSD
                % TRIG.ns6LFP TRIG.ns6CSD_diff TRIG.ns6CSD and TRIG.aMUA
    % b) Trigger to STIM_photo.onsetsdown
        % i) trigger ns2
        % i) trigger ns6

        
% A) Trigger to STIM.onsetsdown 
% remove TPs that are too close to start or end
STIM.onsetsdown(STIM.onsetsdown - pre  < 0) = [];
STIM.onsetsdown(STIM.onsetsdown + post > size(ns6LFPdown,2)) = [];

% get dimensions for preallocation
contactLogical = strcmp({NS2_Header.ElectrodesInfo.ConnectorBank},PARAMS.V1bank{1,1}(2));
contactLabels = {NS2_Header.ElectrodesInfo(contactLogical).Label};
contactNum = size(contactLabels,2);
if contactNum ~= PARAMS.el
    error('Houston we have a problem') %I'm not totally sure this set-up will always work... BMC 8/26/19
end
numTriggers = size(STIM.onsetsdown,1);

%preallocate
TRIG.ns2LFP         = NaN( contactNum,(post + pre + 1),numTriggers); 
TRIG.ns2CSD_diff    = NaN( contactNum,(post + pre + 1),numTriggers);
TRIG.ns2CSD         = NaN( contactNum,(post + pre + 1),numTriggers);
TRIG.ns6LFP         = NaN( contactNum,(post + pre + 1),numTriggers);
TRIG.ns6CSD_diff    = NaN( contactNum,(post + pre + 1),numTriggers);
TRIG.ns6CSD         = NaN( contactNum,(post + pre + 1),numTriggers);
TRIG.aMUA           = NaN( contactNum,(post + pre + 1),numTriggers);

    
%%% DATA MUST BE DOWNSAMPLED TO 1KH IN ORDER TO TRIGGER PROPERLY
for singleCh = 1:contactNum 
    for singleTrigger = 1:numTriggers
        timeOfTrigger = STIM.onsetsdown(singleTrigger);
        windowOfTrigger = timeOfTrigger-pre:timeOfTrigger+post;
        % output is (Ch x time x triggerNumber)
        TRIG.ns2LFP(singleCh,:,singleTrigger)       = ns2LFP(singleCh,windowOfTrigger); 
        TRIG.ns2CSD_diff(singleCh,:,singleTrigger)	= ns2CSD_diffpad(singleCh,windowOfTrigger);
        TRIG.ns2CSD(singleCh,:,singleTrigger)       = ns2CSDpad(singleCh,windowOfTrigger);
        TRIG.ns6LFP(singleCh,:,singleTrigger)       = ns6LFPdown(singleCh,windowOfTrigger);
        TRIG.ns6CSD_diff(singleCh,:,singleTrigger)	= ns6CSD_diffpad(singleCh,windowOfTrigger);
        TRIG.ns6CSD(singleCh,:,singleTrigger)       = ns6CSDpad(singleCh,windowOfTrigger);
        TRIG.aMUA(singleCh,:,singleTrigger)         = aMUAdown(singleCh,windowOfTrigger);
    end
end

% B) Trigger to STIM_photo.onsetsdown 
% remove TPs that are too close to start or end
STIM_photo.onsetsdown(STIM_photo.onsetsdown - pre  < 0) = [];
STIM_photo.onsetsdown(STIM_photo.onsetsdown + post > size(ns6LFPdown,2)) = [];

% get dimensions for preallocation
contactLogical = strcmp({NS2_Header.ElectrodesInfo.ConnectorBank},PARAMS.V1bank{1,1}(2));
contactLabels = {NS2_Header.ElectrodesInfo(contactLogical).Label};
contactNum = size(contactLabels,2);
if contactNum ~= PARAMS.el
    error('Houston we have a problem') %I'm not totally sure this set-up will always work... BMC 8/26/19
end
numTriggers_photo = size(STIM_photo.onsetsdown,1);

%preallocate
TRIG_photo.ns2LFP         = NaN( contactNum,(post + pre + 1),numTriggers_photo); 
TRIG_photo.ns2CSD_diff    = NaN( contactNum,(post + pre + 1),numTriggers_photo);
TRIG_photo.ns2CSD         = NaN( contactNum,(post + pre + 1),numTriggers_photo);
TRIG_photo.ns6LFP         = NaN( contactNum,(post + pre + 1),numTriggers_photo);
TRIG_photo.ns6CSD_diff    = NaN( contactNum,(post + pre + 1),numTriggers_photo);
TRIG_photo.ns6CSD         = NaN( contactNum,(post + pre + 1),numTriggers_photo);
TRIG_photo.aMUA           = NaN( contactNum,(post + pre + 1),numTriggers_photo);

    
%%% DATA MUST BE DOWNSAMPLED TO 1KH IN ORDER TO TRIGGER PROPERLY
for singleCh = 1:contactNum 
    for singleTrigger = 1:numTriggers_photo
        timeOfTrigger_photo = STIM_photo.onsetsdown(singleTrigger);
        windowOfTrigger_photo = timeOfTrigger_photo-pre:timeOfTrigger_photo+post;
        % output is (Ch x time x triggerNumber)
        TRIG_photo.ns2LFP(singleCh,:,singleTrigger)       = ns2LFP(singleCh,windowOfTrigger_photo); 
        TRIG_photo.ns2CSD_diff(singleCh,:,singleTrigger)	= ns2CSD_diffpad(singleCh,windowOfTrigger_photo);
        TRIG_photo.ns2CSD(singleCh,:,singleTrigger)       = ns2CSDpad(singleCh,windowOfTrigger_photo);
        TRIG_photo.ns6LFP(singleCh,:,singleTrigger)       = ns6LFPdown(singleCh,windowOfTrigger_photo);
        TRIG_photo.ns6CSD_diff(singleCh,:,singleTrigger)	= ns6CSD_diffpad(singleCh,windowOfTrigger_photo);
        TRIG_photo.ns6CSD(singleCh,:,singleTrigger)       = ns6CSDpad(singleCh,windowOfTrigger_photo);
        TRIG_photo.aMUA(singleCh,:,singleTrigger)         = aMUAdown(singleCh,windowOfTrigger_photo);
    end
end




%% AVERAGE AND BASELINE-CORRECT TRIGGERED DATA

%% Average
% Avgerage TRIG
fields.TRIG = fieldnames(TRIG);
for avtr=1:numel(fields.TRIG)
    AVG.(fields.TRIG{avtr})  = mean(TRIG.(fields.TRIG{avtr}),3);
end


% Avgerage TRIG_photo
fields.TRIG_photo = fieldnames(TRIG_photo);
for avtrp=1:numel(fields.TRIG_photo)
    AVG_photo.(fields.TRIG_photo{avtrp})  = mean(TRIG_photo.(fields.TRIG_photo{avtrp}),3);
end
%% Baseline corect
% bl AVG
fields.AVG = fieldnames(AVG);
for blavg=1:numel(fields.AVG)
    bl  = mean(AVG.(fields.AVG{blavg})(:,TM<0),2);
    BLavg.(fields.AVG{blavg})  = AVG.(fields.AVG{blavg}) - bl;
end

% bl AVG_photo
fields.AVG_photo = fieldnames(AVG_photo);
for blavgp=1:numel(fields.AVG_photo)
    bl_photo  = mean(AVG_photo.(fields.AVG_photo{blavgp})(:,TM<0),2);
    BLavg_photo.(fields.AVG_photo{blavgp})  = AVG_photo.(fields.AVG_photo{blavgp}) - bl_photo;
end


%% Filter and interpolate CSD, ns2 and ns5
% AVG
AVG.ns2fCSD = filterCSD(AVG.ns2CSD);
AVG.ns6fCSD = filterCSD(AVG.ns6CSD);

% AVG_photo
AVG_photo.ns2fCSD = filterCSD(AVG_photo.ns2CSD);
AVG_photo.ns6fCSD = filterCSD(AVG_photo.ns6CSD);

% BLavg
BLavg.ns2fCSD       = filterCSD(BLavg.ns2CSD);
BLavg.ns2fCSD_diff  = filterCSD(BLavg.ns2CSD_diff);
BLavg.ns6fCSD       = filterCSD(BLavg.ns6CSD);
BLavg.ns6fCSD_diff  = filterCSD(BLavg.ns6CSD_diff);

% BLavg_photo
BLavg_photo.ns2fCSD       = filterCSD(BLavg_photo.ns2CSD);
BLavg_photo.ns2fCSD_diff  = filterCSD(BLavg_photo.ns2CSD_diff);
BLavg_photo.ns6fCSD       = filterCSD(BLavg_photo.ns6CSD);
BLavg_photo.ns6fCSD_diff  = filterCSD(BLavg_photo.ns6CSD_diff);


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
f_ShadedLinePlotbyDepth(BLavg_photo.ns6LFP,chans,TM,[],1)
plot([0 0], ylim,'k')
title({'LFP',filename}, 'Interpreter', 'none')
xlabel('time (ms)')
ylabel('Contacts indexed down from surface')

% CSD line
subplot(1,4,2)
f_ShadedLinePlotbyDepth(BLavg_photo.ns6CSD,chans,TM,[],1)
plot([0 0], ylim,'k')
title({'CSD',filename}, 'Interpreter', 'none')
xlabel('time (ms)')

% fCSD
subplot(1,4,3)
imagesc(TM,chans,BLavg_photo.ns6fCSD); 
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
f_ShadedLinePlotbyDepth(BLavg_photo.aMUA,chans,TM,[],1)
plot([0 0], ylim,'k')
title({'aMUA',filename}, 'Interpreter', 'none')
xlabel('time (ms)')

%% SCRIPT END
toc

load gong.mat;
sound(y);




