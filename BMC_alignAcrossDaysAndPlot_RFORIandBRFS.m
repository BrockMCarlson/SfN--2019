%BMC_alignAcrossDaysAndPlot_RFORIandBRFS.m
%GOAL: align and average across days the PEvsNPEvsPSvsNPS. The 2x2plot. The
%RFORI sessions contain the DE response and the BRFS contain the NDE
%response under monocular conditions.
%
%   Version 1.0
%   Brock Carlson -- created 9/16/19
%   
%   Current plotting goal: plot rfori for each orientation presented in
%   brfs and the two orientations of brfs under monocular conditions.
%
%   DOES NOT TRIGGER TO PHOTO DIODE



clear

%% EDITABLE VARIABLES
% % filename = {'160102_E_brfs001'};
filename = {'160102_E_rfori002','160102_E_brfs001','160427_E_rfori002','160427_E_brfs001','160510_E_rfori002','160510_E_brfs001'}';
sinkAllocate = 'BMC_DfS';
pre = 50;
post = 250;
TM = -pre:1:post;
nameSaveType = 'LFPandCSDof';
figtype = 'AlignedSessionEyeandOri';



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
    clearvars -except a filename sinkAllocate pre post TM savefiledir nameSaveType figtype ALIGNED
    
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

if contains(filename{a},'brfs')
    STIM_BRFS = sortBrfsStimandTimeData(grating,pEvC,pEvT,PARAMS);
elseif contains(filename{a},'rfori')
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


%% trig brfs onsets or rfori onsets
if contains(filename{a},'brfs')
    fields.STIM_BRFS = fieldnames(STIM_BRFS);
    for ffstim=1:numel(fields.STIM_BRFS)

        if sum(strcmp(fieldnames(STIM_BRFS.(fields.STIM_BRFS{ffstim})),'start_noSoaDown')) == 1
            numTriggers.STIM_BRFS.(fields.STIM_BRFS{ffstim}) = size(STIM_BRFS.(fields.STIM_BRFS{ffstim}).start_noSoaDown,1);
            %preallocate
            TRIG_BRFS.(fields.STIM_BRFS{ffstim}).ns2LFP         = NaN( contactNum,(post + pre + 1),numTriggers.STIM_BRFS.(fields.STIM_BRFS{ffstim})); 
            TRIG_BRFS.(fields.STIM_BRFS{ffstim}).ns2CSD         = NaN( contactNum,(post + pre + 1),numTriggers.STIM_BRFS.(fields.STIM_BRFS{ffstim}));
            %%% DATA MUST BE DOWNSAMPLED TO 1KH IN ORDER TO TRIGGER PROPERLY
            for singleCh = 1:contactNum 
                for singleTrigger = 1: numTriggers.STIM_BRFS.(fields.STIM_BRFS{ffstim})
                    timeOfTrigger = STIM_BRFS.(fields.STIM_BRFS{ffstim}).start_noSoaDown(singleTrigger);
                    windowOfTrigger = timeOfTrigger-pre:timeOfTrigger+post;
                    % output is (Ch x time x triggerNumber)
                    TRIG_BRFS.(fields.STIM_BRFS{ffstim}).ns2LFP(singleCh,:,singleTrigger)       = ns2LFP(singleCh,windowOfTrigger); 
                    TRIG_BRFS.(fields.STIM_BRFS{ffstim}).ns2CSD(singleCh,:,singleTrigger)       = ns2CSD(singleCh,windowOfTrigger);

                end
            end

        elseif sum(strcmp(fieldnames(STIM_BRFS.(fields.STIM_BRFS{ffstim})),'start2Down')) == 1
            numTriggers.STIM_BRFS.(fields.STIM_BRFS{ffstim}) = size(STIM_BRFS.(fields.STIM_BRFS{ffstim}).start2Down,1);
            %preallocate
            TRIG_BRFS.(fields.STIM_BRFS{ffstim}).ns2LFP         = NaN( contactNum,(post + pre + 1),numTriggers.STIM_BRFS.(fields.STIM_BRFS{ffstim})); 
            TRIG_BRFS.(fields.STIM_BRFS{ffstim}).ns2CSD         = NaN( contactNum,(post + pre + 1),numTriggers.STIM_BRFS.(fields.STIM_BRFS{ffstim}));
            %%% DATA MUST BE DOWNSAMPLED TO 1KH IN ORDER TO TRIGGER PROPERLY
            for singleCh = 1:contactNum 
                for singleTrigger = 1: numTriggers.STIM_BRFS.(fields.STIM_BRFS{ffstim})
                    timeOfTrigger = STIM_BRFS.(fields.STIM_BRFS{ffstim}).start2Down(singleTrigger);
                    windowOfTrigger = timeOfTrigger-pre:timeOfTrigger+post;
                    % output is (Ch x time x triggerNumber)
                    TRIG_BRFS.(fields.STIM_BRFS{ffstim}).ns2LFP(singleCh,:,singleTrigger)       = ns2LFP(singleCh,windowOfTrigger); 
                    TRIG_BRFS.(fields.STIM_BRFS{ffstim}).ns2CSD(singleCh,:,singleTrigger)       = ns2CSD(singleCh,windowOfTrigger);

                end
            end


        end       
    end
    % Concatonate all monocular conditions (i.e. average monocular PS,
    % diop_800soa_PS, and dichop_800soa_NPS as the PS monocular condition. Ten
    % average monocular NPS, diop_800soa_NPS, and dichop_800soa_PS as the NPS
    % monocualr condition)     
        TRIG_BRFS.allmonocPS.ns2CSD =  cat(3,TRIG_BRFS.monoc_PS.ns2CSD,TRIG_BRFS.diop_800soa_PS.ns2CSD);
        TRIG_BRFS.allmonocPS.ns2CSD =  cat(3,TRIG_BRFS.allmonocPS.ns2CSD,TRIG_BRFS.dichop_800soa_brfsNPSflash.ns2CSD);
        TRIG_BRFS.allmonocNPS.ns2CSD =  cat(3,TRIG_BRFS.monoc_NPS.ns2CSD,TRIG_BRFS.diop_800soa_NPS.ns2CSD);
        TRIG_BRFS.allmonocNPS.ns2CSD =  cat(3,TRIG_BRFS.allmonocNPS.ns2CSD,TRIG_BRFS.dichop_800soa_brfsPSflash.ns2CSD);
        %lfp needed for proper indexing later
        TRIG_BRFS.allmonocPS.ns2LFP =  cat(3,TRIG_BRFS.monoc_PS.ns2LFP,TRIG_BRFS.diop_800soa_PS.ns2LFP);
        TRIG_BRFS.allmonocPS.ns2LFP =  cat(3,TRIG_BRFS.allmonocPS.ns2LFP,TRIG_BRFS.dichop_800soa_brfsNPSflash.ns2LFP);
        TRIG_BRFS.allmonocNPS.ns2LFP =  cat(3,TRIG_BRFS.monoc_NPS.ns2LFP,TRIG_BRFS.diop_800soa_NPS.ns2LFP);
        TRIG_BRFS.allmonocNPS.ns2LFP =  cat(3,TRIG_BRFS.allmonocNPS.ns2LFP,TRIG_BRFS.dichop_800soa_brfsPSflash.ns2LFP);

elseif contains(filename{a},'rfori')
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

end


%% AVERAGE AND BASELINE-CORRECT TRIGGERED DATA 

%% Average
% Avgerage TRIG
fields.TRIG = fieldnames(TRIG);
for avtr=1:numel(fields.TRIG)
    AVG.(fields.TRIG{avtr})  = mean(TRIG.(fields.TRIG{avtr}),3);
end

if contains(filename{a},'brfs')
    firstfields.TRIG_BRFS = fieldnames(TRIG_BRFS);
    for ffav = 1:numel(firstfields.TRIG_BRFS)
        subfields.TRIG_BRFS = fieldnames(TRIG_BRFS.(firstfields.TRIG_BRFS{ffav}));
        for avtr=1:numel(subfields.TRIG_BRFS)
            AVG_BRFS.(firstfields.TRIG_BRFS{ffav}).(subfields.TRIG_BRFS{avtr})  = mean(TRIG_BRFS.(firstfields.TRIG_BRFS{ffav}).(subfields.TRIG_BRFS{avtr}),3);
        end

    end
end



%% Baseline corect
% bl AVG
fields.AVG = fieldnames(AVG);
for blavg=1:numel(fields.AVG)
    bl  = mean(AVG.(fields.AVG{blavg})(:,TM<0),2);
    BLavg.(fields.AVG{blavg})  = AVG.(fields.AVG{blavg}) - bl;
end

if contains(filename{a},'brfs')
    firstfields.AVG_BRFS = fieldnames(AVG_BRFS);
    for ffbl = 1:numel(firstfields.AVG_BRFS)
        subfields.AVG_BRFS = fieldnames(AVG_BRFS.(firstfields.AVG_BRFS{ffbl}));
        for blavg=1:numel(subfields.AVG_BRFS)
            bl  = mean(AVG_BRFS.(firstfields.AVG_BRFS{ffbl}).(subfields.AVG_BRFS{blavg})(:,TM<0),2);
            BLavg_BRFS.(firstfields.AVG_BRFS{ffbl}).(subfields.AVG_BRFS{blavg})  = AVG_BRFS.(firstfields.AVG_BRFS{ffbl}).(subfields.AVG_BRFS{blavg}) - bl;
        end
    end
end

%% Align Across Days -- sink on ch 25


% cut the contacts' index based on the cortex limits (rftop and rfbtm)

    if contains(filename{a},'brfs')
        firstfields.TRIG_BRFS = fieldnames(TRIG_BRFS);
        for ffav = 1:numel(firstfields.TRIG_BRFS)
            subfields.TRIG_BRFS = fieldnames(TRIG_BRFS.(firstfields.TRIG_BRFS{ffav}));
            for avtr=1:numel(subfields.TRIG_BRFS)
                cutmatrix.(firstfields.TRIG_BRFS{ffav}).(fields.AVG{avtr})  = AVG.(fields.AVG{avtr})(PARAMS.rftop:PARAMS.rfbtm,:);
                ALIGNED(a).(firstfields.TRIG_BRFS{ffav}).(fields.AVG{avtr}) = nan(50,size(TM,2));
                % creat row assignments for the ALIGNED matrix
                clear x y z
                x = PARAMS.EvalSink-PARAMS.rftop+1; % the sink's row# in cutmatrix
                y = 25-x+1; %rftop's location in ALIGNED
                z = y+size(cutmatrix.(firstfields.TRIG_BRFS{ffav}).(fields.AVG{avtr}),1)-1; %rfbtm's locaiton in ALIGNED.    
                ALIGNED(a).(firstfields.TRIG_BRFS{ffav}).(fields.AVG{avtr})(y:z,:) = cutmatrix.(firstfields.TRIG_BRFS{ffav}).(fields.AVG{avtr}) ;                      
            end
        end
    elseif contains(filename{a},'rfori')
        fields.AVG = fieldnames(AVG);
        for avtr=1:numel(fields.AVG)
            cutmatrix.rfori.(fields.AVG{avtr})  = AVG.(fields.AVG{avtr})(PARAMS.rftop:PARAMS.rfbtm,:);
            ALIGNED(a).rfori.(fields.AVG{avtr}) = nan(50,size(TM,2));
            % creat row assignments for the ALIGNED matrix
            x = PARAMS.EvalSink-PARAMS.rftop+1; % the sink's row# in cutmatrix
            y = 25-x+1; %rftop's location in ALIGNED
            z = y+size(cutmatrix.rfori.(fields.AVG{avtr}),1)-1; %rfbtm's locaiton in ALIGNED.    
            ALIGNED(a).rfori.(fields.AVG{avtr})(y:z,:) = cutmatrix.rfori.(fields.AVG{avtr}) ;
        end
    end




end
clearvars -except  filename sinkAllocate pre post TM savefiledir nameSaveType figtype ALIGNED


%% Average across alignment
%
% %%%%%% I'm having a problem here. How do I average across sessions? If I
% still split rfor and brfs it might be easier. But may also be a headache
% down the road. Should I only put in one file type at a time in the
% filenames?
% %%%%%% For now, I will simply focus on rfori and the allmonoc fields
% %%%%%% I'm realizing that this was a bit of sloppy coding. May have to
%        re-do in the future.
%

% PExPS -- rfori ori1
PExPS = nan(50,size(TM,2),3);
PExPS(:,:,1) = ALIGNED(1).rfori.ori1CSD;
PExPS(:,:,2) = ALIGNED(3).rfori.ori1CSD;
PExPS(:,:,3) = ALIGNED(5).rfori.ori1CSD;
AlAvg.PExPS = nanmean(PExPS,3);

% PExNPS -- rfori ori2
PExNPS = nan(50,size(TM,2),3);
PExNPS(:,:,1) = ALIGNED(1).rfori.ori2CSD;
PExNPS(:,:,2) = ALIGNED(3).rfori.ori2CSD;
PExNPS(:,:,3) = ALIGNED(5).rfori.ori2CSD;
AlAvg.PExNPS = nanmean(PExNPS,3);

% NPExPS -- brfs allmonocPS
NPExPS = nan(50,size(TM,2),3);
NPExPS(:,:,1) = ALIGNED(2).allmonocPS.ns2CSD;
NPExPS(:,:,2) = ALIGNED(4).allmonocPS.ns2CSD;
NPExPS(:,:,3) = ALIGNED(6).allmonocPS.ns2CSD;
AlAvg.NPExPS = nanmean(NPExPS,3);

% NPExNPS -- brfs allmonocNPS
NPExNPS = nan(50,size(TM,2),3);
NPExNPS(:,:,1) = ALIGNED(2).allmonocNPS.ns2CSD;
NPExNPS(:,:,2) = ALIGNED(4).allmonocNPS.ns2CSD;
NPExNPS(:,:,3) = ALIGNED(6).allmonocNPS.ns2CSD;
AlAvg.NPExNPS = nanmean(NPExNPS,3);

%% Baseline average the alignment

clear bl
bl.PExPS  = mean(AlAvg.PExPS(:,TM<0),2);
bl.PExNPS  = mean(AlAvg.PExNPS(:,TM<0),2);
bl.NPExPS  = mean(AlAvg.NPExPS(:,TM<0),2);
bl.NPExNPS  = mean(AlAvg.NPExNPS(:,TM<0),2);
AlBLavg.PExPS  = AlAvg.PExPS - bl.PExPS;
AlBLavg.PExNPS  = AlAvg.PExNPS - bl.PExNPS;
AlBLavg.NPExPS  = AlAvg.NPExPS - bl.NPExPS;
AlBLavg.NPExNPS  = AlAvg.NPExNPS - bl.NPExNPS;

%% Cut matrix for cortical depth
AlCut.PExPS = AlBLavg.PExPS(14:31,:);
AlCut.PExNPS = AlBLavg.PExNPS(14:31,:);
AlCut.NPExPS = AlBLavg.NPExPS(14:31,:);
AlCut.NPExNPS = AlBLavg.NPExNPS(14:31,:);

%% Filter and interpolate CSD across averaged alignment
% filter AVG struct
AlFilt.PExPS = filterCSD(AlCut.PExPS);
AlFilt.PExNPS = filterCSD(AlCut.PExNPS);
AlFilt.NPExPS = filterCSD(AlCut.NPExPS);
AlFilt.NPExNPS = filterCSD(AlCut.NPExNPS);






%% Plot
corticaldepth = (1.1:-0.1:-0.5);
climit = 1000;

figure;
set(gcf, 'Position', [680 338 711 760]);
subplot(2,2,1)
imagesc(TM,corticaldepth,AlFilt.PExPS); 
colormap(flipud(jet));
% % % set(gca,'CLim',[-climit{a} climit{a}],'Box','off','TickDir','out')
% climit = max(abs(get(gca,'CLim'))*.8);
set(gca,'CLim',[-climit climit],'YDir','normal','Box','off','TickDir','out')
hold on;
plot([0 0], ylim,'k')
plot([800 800], ylim,'k')
clrbar = colorbar;
title({'Preferred eye. Preferred stim'}, 'Interpreter', 'none')
ylabel('cortical depth')
xlabel('time (ms)')
clrbar.Label.String = 'nA/mm^3';


set(gcf, 'Position', [680 338 711 760]);
subplot(2,2,2)
imagesc(TM,corticaldepth,AlFilt.PExNPS); 
colormap(flipud(jet));
% % % set(gca,'CLim',[-climit{a} climit{a}],'Box','off','TickDir','out')
% climit = max(abs(get(gca,'CLim'))*.8);
set(gca,'CLim',[-climit climit],'YDir','normal','Box','off','TickDir','out')
hold on;
plot([0 0], ylim,'k')
plot([800 800], ylim,'k')
clrbar = colorbar;
title({'Preferred eye. Non-Pref stim'}, 'Interpreter', 'none')
xlabel('time (ms)')
clrbar.Label.String = 'nA/mm^3';

set(gcf, 'Position', [680 338 711 760]);
subplot(2,2,3)
imagesc(TM,corticaldepth,AlFilt.NPExPS); 
colormap(flipud(jet));
% % % set(gca,'CLim',[-climit{a} climit{a}],'Box','off','TickDir','out')
% climit = max(abs(get(gca,'CLim'))*.8);
set(gca,'CLim',[-climit climit],'YDir','normal','Box','off','TickDir','out')
hold on;
plot([0 0], ylim,'k')
plot([800 800], ylim,'k')
clrbar = colorbar;
title({'Non-pref eye. Preferred stim'}, 'Interpreter', 'none')
ylabel('cortical depth')
xlabel('time (ms)')
clrbar.Label.String = 'nA/mm^3';


set(gcf, 'Position', [680 338 711 760]);
subplot(2,2,4)
imagesc(TM,corticaldepth,AlFilt.NPExNPS); 
colormap(flipud(jet));
% % % set(gca,'CLim',[-climit{a} climit{a}],'Box','off','TickDir','out')
% climit = max(abs(get(gca,'CLim'))*.8);
set(gca,'CLim',[-climit climit],'YDir','normal','Box','off','TickDir','out')
hold on;
plot([0 0], ylim,'k')
plot([800 800], ylim,'k')
clrbar = colorbar;
title({'Non-pref eye. Non-pref stim'}, 'Interpreter', 'none')
xlabel('time (ms)')
clrbar.Label.String = 'nA/mm^3';

set(gcf, 'Position',[680 54 711 1044]);


%% SAVE figs
cd(savefiledir)
figsavename = strcat(figtype);
saveas(gcf, figsavename, 'fig')
saveas(gcf, figsavename, 'pdf')
saveas(gcf, figsavename, 'png')
    








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





