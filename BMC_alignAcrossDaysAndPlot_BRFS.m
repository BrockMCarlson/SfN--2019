%BMC_alignAcrossDaysAndPlot_BRFS-ONLY.m
%GOAL: align and average across days the PEvsNPEvsPSvsNPS. The 2x2plot. The
%RFORI sessions contain the DE response and the BRFS contain the NDE
%response under monocular conditions.
%
%   Version 1.0
%   Brock Carlson -- created 9/17/19
%   
%   Current plotting goal: plot brfs PS after flash vs null. This code was
%   taken from the RFORIandBRFS alignment code and removed all rfori
%   elements. Hopefully this will clean it up.
%
%   DOES NOT TRIGGER TO PHOTO DIODE
%
%   Version 1.2 - SfN fig2 draft



clear

%% EDITABLE VARIABLES
% % filename = {'160102_E_brfs001'};
filename = {'160102_E_brfs001','160427_E_brfs001','160510_E_brfs001'}';
sinkAllocate = 'BMC_DfS';
pre = 850;
post = 800;
TM = -pre:1:post;
nameLoadType = 'LFPandCSDof';
figtype = 'SfN_fig2';



% Computer-specific editable variables 
[ret, hostname] = system('hostname');
if strcmp(getenv('USER'),'maierav')
    % @Alex -- fill in necessary information for your system here.
        % loadfiledir = 'G:\LaCie\SfN 2019--figsAndMatVars\SfN 2019 figs\brfs conditions diagnostics'; 
         sessionParamDir = 'G:\LaCie';
elseif strcmp(getenv('USERNAME'),'Brock Carlson')
    % variables for end of script
    loadfiledir =  'G:\LaCie\SfN 2019--figsAndMatVars\SfN 2019 figs\SfN_fig2';
    sessionParamDir = 'G:/LaCie';
elseif  ~ispc && contains(hostname,'Brocks-MacBook-Air')
    loadfiledir = '/Volumes/SfN_2019/SfN 2019 MacBook Figs/SfN_fig2 -- MacBook';
    sessionParamDir = '/Volumes/SfN_2019/';
end





for a = 1:size(filename,1)
    clearvars -except a filename sinkAllocate pre post TM loadfiledir nameLoadType figtype ALIGNED sessionParamDir
    
    disp(filename{a})
    

%% Computer-specific directories 
[ret, hostname] = system('hostname');
if strcmp(getenv('USER'),'maierav')
    % @Alex -- fill in necessary information for your system here.
        addpath(genpath('/Users/alex 1/Desktop/LAB/Brock'));
    dataDirectory = strcat('/Users/alex 1/Desktop/LAB/Brock',filesep,filename{a}(1:8));    
elseif strcmp(getenv('USERNAME'),'Brock Carlson')
    addpath(genpath('G:\LaCie\all BRFS'));
    dataDirectory = strcat('G:\LaCie\all BRFS\',filename{a}(1:8));    
elseif ~ispc && contains(hostname,'Brocks-MacBook-Air')
    dataDirectory = strcat('/Volumes/SfN_2019/all BRFS/',filename{a}(1:8));   
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
cd(loadfiledir)
loadname = strcat(nameLoadType,filename{a},'.mat');
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

    STIM_BRFS = sortBrfsStimandTimeData(grating,pEvC,pEvT,PARAMS);


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






%% Align Across Days -- sink on ch 25


% cut the contacts' index based on the cortex limits (rftop and rfbtm)


        firstfields.AVG_BRFS = fieldnames(AVG_BRFS);
        for ffal = 1:numel(firstfields.AVG_BRFS)
            subfields.AVG_BRFS = fieldnames(AVG_BRFS.(firstfields.AVG_BRFS{ffal}));
            for sfal=1:numel(subfields.AVG_BRFS)
                cutmatrix.(firstfields.AVG_BRFS{ffal}).(subfields.AVG_BRFS{sfal})  = AVG_BRFS.(firstfields.AVG_BRFS{ffal}).(subfields.AVG_BRFS{sfal})(PARAMS.rftop:PARAMS.rfbtm,:);
                ALIGNED(a).(firstfields.AVG_BRFS{ffal}).(subfields.AVG_BRFS{sfal}) = nan(50,size(TM,2));
                % creat row assignments for the ALIGNED matrix
                clear x y z
                x = PARAMS.EvalSink-PARAMS.rftop+1; % the sink's row# in cutmatrix
                y = 25-x+1; %rftop's location in ALIGNED
                z = y+size(cutmatrix.(firstfields.AVG_BRFS{ffal}).(subfields.AVG_BRFS{sfal}),1)-1; %rfbtm's locaiton in ALIGNED.    
                ALIGNED(a).(firstfields.AVG_BRFS{ffal}).(subfields.AVG_BRFS{sfal})(y:z,:) = cutmatrix.(firstfields.AVG_BRFS{ffal}).(subfields.AVG_BRFS{sfal}) ;                      
            end
        end




end


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


% PSflash -- brfs dichop_800soa_brfsPSflash
PSflash = nan(50,size(TM,2),3);
PSflash(:,:,1) = ALIGNED(1).dichop_800soa_brfsPSflash.ns2CSD;
PSflash(:,:,2) = ALIGNED(2).dichop_800soa_brfsPSflash.ns2CSD;
PSflash(:,:,3) = ALIGNED(3).dichop_800soa_brfsPSflash.ns2CSD;
AlAvg.PSflash = nanmean(PSflash,3);

% nullFlash -- brfs dichop_800soa_brfsNPSflash
nullFlash = nan(50,size(TM,2),3);
nullFlash(:,:,1) = ALIGNED(1).dichop_800soa_brfsNPSflash.ns2CSD;
nullFlash(:,:,2) = ALIGNED(2).dichop_800soa_brfsNPSflash.ns2CSD;
nullFlash(:,:,3) = ALIGNED(3).dichop_800soa_brfsNPSflash.ns2CSD;
AlAvg.nullFlash = nanmean(nullFlash,3);

%% Baseline average the alignment
searchvector = (-50:1:0);
[~,blidx] = ismember(searchvector,TM);


clear bl
bl.PSflash  = mean(AlAvg.PSflash(:,blidx),2);
bl.nullFlash  = mean(AlAvg.nullFlash(:,blidx),2);

AlBLavg.PSflash  = AlAvg.PSflash - bl.PSflash;
AlBLavg.nullFlash  = AlAvg.nullFlash - bl.nullFlash;

%% Cut matrix for cortical depth

AlCut.PSflash = AlBLavg.PSflash(14:31,:);
AlCut.nullFlash = AlBLavg.nullFlash(14:31,:);

%% Filter and interpolate CSD across averaged alignment
% filter AVG struct

AlFilt.PSflash = filterCSD(AlCut.PSflash);
AlFilt.nullFlash = filterCSD(AlCut.nullFlash);






%% Plot
corticaldepth = (1.1:-0.1:-0.5);
% % climit = 1000;

figure;
subplot(1,2,1)
imagesc(TM,corticaldepth,AlFilt.PSflash); 
colormap(flipud(jet));
% % % set(gca,'CLim',[-climit{a} climit{a}],'Box','off','TickDir','out')
climit = max(abs(get(gca,'CLim'))*.8);
set(gca,'CLim',[-climit climit],'YDir','normal','Box','off','TickDir','out')
hold on;
plot([0 0], ylim,'k')
plot([-800 -800], ylim,'k')
clrbar = colorbar;
title({'BRFS. PS flashed'}, 'Interpreter', 'none')
ylabel('cortical depth')
xlabel('time (ms)')
clrbar.Label.String = 'nA/mm^3';


subplot(1,2,2)
imagesc(TM,corticaldepth,AlFilt.nullFlash); 
colormap(flipud(jet));
% % % set(gca,'CLim',[-climit{a} climit{a}],'Box','off','TickDir','out')
climit = max(abs(get(gca,'CLim'))*.8);
set(gca,'CLim',[-climit climit],'YDir','normal','Box','off','TickDir','out')
hold on;
plot([0 0], ylim,'k')
plot([-800 -800], ylim,'k')
clrbar = colorbar;
title({'BRFS. Null Flashed'}, 'Interpreter', 'none')
xlabel('time (ms)')
clrbar.Label.String = 'nA/mm^3';


set(gcf, 'Position',[151 397 1694 701]);


%% SAVE figs
cd(loadfiledir)
figsavename = strcat(figtype);
saveas(gcf, figsavename, 'fig')
saveas(gcf, figsavename, 'pdf')
saveas(gcf, figsavename, 'png')
    










