% Make CSD of rfori
clear

%% EDITABLE VARIABLES
filename = '160102_E_rfori002';
directory = 'G:\LaCie\all BRFS\160102_E';
sinkAllocate = 'BMC_DfS';
pre = 50;
post = 250;
subBaseline = true;
manualUv = true;

%% LOAD AND ORDER PINS
%Load session params
cd('G:\LaCie')
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

%make names to load later
cd(directory)
readNEVfile = strcat(directory,filesep,filename,'.nev');
readNS2file = strcat(directory,filesep,filename,'.ns2');
readNS6file = strcat(directory,filesep,filename,'.ns6');
readGRATINGfile = strcat(directory,filesep,filename,'.gRFORIGrating_di');

%nev loading and set-up
NEV = openNEV(readNEVfile); 

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
    error('Houston we have a problem')
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
    
    
%% PREPROCESS, 
% filter data and calculate CSD
    %filter  LFP
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
    
    %filter aMUA
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
    
    %Downsample aMUA to 1kHz to aling with triggerpoints below
    aMUAtrans = aMUA';
    aMUAdownTrans = downsample(aMUAtrans,30);
    aMUAdown = aMUAdownTrans';

%CSD
    %derivitive 
    ns2CSD_diff = diff(ns2LFP,2,1);
    ns6CSD_diff = diff(ns6LFP,2,1);
    
    %calcCSD_BMC    
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
    ns6CSD(:,:) = (-1*calced*ns6LFP);

    
%% TRIGGER TO DIODE ====> Problem Solved perhaps?
% trigger to BNC and to pEvC/pEvT to double-check against photodiode
% get triggered LFP, CSD, and aMUA out.

% Trigger to BNC
EventCodes = NEV.Data.SerialDigitalIO.UnparsedData - 128;
EventTimes = double(NEV.Data.SerialDigitalIO.TimeStamp); %30kHz samples
grating = readgGrating(readGRATINGfile);
    [pEvC,pEvT] = parsEventCodesML(EventCodes,EventTimes);
    [pEvT_photo,tf] = pEvtPhoto2(readGRATINGfile,pEvC,pEvT,mode(grating.ypos),[],'ainp1',0);

%% TRIGGER TO pEvC/pEvT at 1 kHz 
EventCodes = NEV.Data.SerialDigitalIO.UnparsedData - 128;
EventTimes = floor(NEV.Data.SerialDigitalIO.TimeStampSec .* 1000); %ms, to match 1kHz
triggerpoints = EventTimes(EventCodes == 23 | EventCodes == 25 | EventCodes == 27 | EventCodes == 29| EventCodes == 31);
    % remove TPs that are too close to start or end
    triggerpoints(triggerpoints - pre  < 0) = [];
    triggerpoints(triggerpoints + post > size(ns2CSD,2)) = [];
    
 
    %32 channel loop.
    numChanaMUA = size(aMUAdown, 1);
    numTriggers = size(triggerpoints,2);
    aMUAtrig    = NaN( numChanaMUA,[post + pre + 1],numTriggers);

    %%% DATA MUST BE DOWNSAMPLED TO 1KH IN ORDER TO TRIGGER PROPERLY
    for singleCh = 1:numChanaMUA 
        for singleTrigger = 1:numTriggers
            timeOfTrigger = triggerpoints(singleTrigger);
            windowOfTrigger = timeOfTrigger-pre:timeOfTrigger+post;
            aMUAtrig(singleCh,:,singleTrigger)= aMUAdown(singleCh,windowOfTrigger);   
            % output is Ch x time x trial/triggerNumber
        end
    end
    aMUAavg = mean(aMUAtrig,3); 
    if subBaseline == true
        aMUAbl  = mean(aMUAavg(:,csdTM<0),2);
        aMUAavg = aMUAavg - aMUAbl;        
    end

    
    
    
    
    %CSD 30 channel loop
    numChan = size(ns2CSD, 1);
    numTriggers   = length(triggerpoints);
    csdTrig     = NaN( numChan,(post + pre + 1),numTriggers); %changed to be (ch x sample x tr)
    aMUAtrig    = NaN( numChan,(post + pre + 1),numTriggers);
    csdTM     = (0:(size(csdTrig,2)-1)) - pre;
    %loop through each channel and trigger each trial
    for singleCh = 1:numChan
        for singleTrigger = 1:numTriggers
            timeOfTrigger = triggerpoints(singleTrigger);
            timeOfTrigger = round(timeOfTrigger);
            windowOfTrigger = timeOfTrigger-pre:timeOfTrigger+post;
            csdTrig(singleCh,:,singleTrigger) = ns2CSD(singleCh,windowOfTrigger); %flipped the dimension from trigData.m. Now Channel is first and time is second
            % output is Ch x time x trial/triggerNumber
        end
    end
    csdAvg  = mean(csdTrig,3);
    if subBaseline == true
        csdBL  = mean(csdAvg(:,csdTM<0),2);
        csdAvg = csdAvg - csdBL;       
    end
    csdPad = padarray(csdAvg,[1 0],NaN,'replicate');
    csdFilter = filterCSD(csdPad);
    


   

    
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
subplot(1,2,1)
imagesc(csdTM,chans,csdFilter); colormap(flipud(jet));
climit = max(abs(get(gca,'CLim'))*.8);
if PARAMS.SortDirection == 'descending'
    ydir = 'reverse';
elseif PARAMS.SortDirection == 'ascending'
    ydir = 'normal';
end
set(gca,'CLim',[-climit climit],'Ydir',ydir,'Box','off','TickDir','out')
hold on;
plot([0 0], ylim,'k')
clrbar = colorbar;
title({'CSD code rewrite',filename}, 'Interpreter', 'none')
xlabel('time (ms)')
ylabel('Contacts indexed down from surface')
clrbar.Label.String = 'nA/mm^3';
set(gcf,'Position',[1 40 331 662]);
%test

subplot(1,2,2)
f_ShadedLinePlotbyDepth(aMUAavg,chans,csdTM,[],1)
plot([0 0], ylim,'k')
title({'aMUA code rewrite',filename}, 'Interpreter', 'none')
xlabel('time (ms)')
ylabel('Contacts indexed down from surface')
set(gcf,'Position',[1 40 331 662]);

