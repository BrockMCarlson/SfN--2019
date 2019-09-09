%% BMC - find Preffered Stim in brfs from discretized MUA spiking data
% combination of BMC conditions (looking at monocular) and KD ppnev
% analysis.

% Update 6/12/2019

clear 
close all

brdrname      = 'G:\LaCie\BRFSdays ppnev rfori and brfs';
cd(brdrname)
BRdatafile    = '160102_E_rfori003'; 
filename      = [brdrname filesep BRdatafile]; 

% get SessionParams
cd('G:\LaCie')
load('SessionParams.mat')
SinkAllocate = 'BMC';
switch SinkAllocate
    case 'BMC'
        SessionParams.EvalSink = SessionParams.BMC_DfS;
    case 'Old'
        SessionParams.EvalSink = SessionParams.Old_DfS;
end
        for j = 1:size(SessionParams.Date,1)
      
        pattern = string(SessionParams.Date(j));
        dateFound = strfind(BRdatafile,pattern);
            if isempty(dateFound) == 0
                idxSessionParams = j;
            end
        end
 sink = SessionParams.EvalSink(idxSessionParams);
%  sink = sink-2:1:sink+1;

% load ppNEV
load(strcat(filename,'.ppnev'),'-MAT');
allchanIDs  = {ppNEV.ElectrodesInfo.ElectrodeLabel};  % every pin on entire 128 channel system 
spikechs    = nanunique(ppNEV.Data.Spikes.Electrode); % every pin with spike data 

% load grating info
cd('G:\LaCie\all BRFS\160102_E')
fileForGrating = strcat(BRdatafile,'.gRFORIGrating_di');
Grating = readgGrating(fileForGrating);
% fileForGrating = strcat(BRdatafile,'.gBrfsGratings');
% Grating = readBRFS(fileForGrating);


            

% the below code will sort the pins in order from channel label 01 to 24 or
% 32. ADJUST FOR NN ARRAYS 
for e = 1:length(spikechs)
    chans(:,e) = allchanIDs{spikechs(e)}(1:4)';  %#ok<SAGROW>
end
els  = nanunique(chans(2,:)); 
nums = sort(nanunique(str2num(chans(3:4,:)'))); %#ok<ST2NM>

pinorder = []; 
for e = 1:length(els)
    for n = 1:length(nums)
    elname    = sprintf('e%s%02u',els(e),nums(n));
    pinorder  = [pinorder find(ismember(chans',elname,'rows'))];   %#ok<AGROW>
    end
end
orderedpins = spikechs(pinorder); 

% get event codes from NEV
clear chans; 
EventCodes    = ppNEV.Data.SerialDigitalIO.UnparsedData - 128;
EventTimes    = floor(ppNEV.Data.SerialDigitalIO.TimeStampSec .* 1000); %ms, to match 1kHz
EventSamples  = ppNEV.Data.SerialDigitalIO.TimeStamp;
Fs            = ppNEV.MetaTags.TimeRes; 
[pEvC,pEvT] = parsEventCodesML(EventCodes,EventTimes);
stim1 = 0;
stim2 = 0;
%%
for t = 1:length(pEvC)

if ~any(pEvC{t} == 96) % This is not necessary on the evp trails
     % skip if trial was aborted and animal was not rewarded (event code 96)
     continue
 end     

stimon   =  pEvC{t} == 23;
stimoff  =  pEvC{t} == 24;

%only pull and save the monocular data
stimQuant = find(stimon);
if	numel(stimQuant) == 1     % there is no soa.
    start_noSoa	=  pEvT{t}(stimon);
    stop        =  pEvT{t}(stimoff);

    % Assign monocular timepoint == without soa.
    % create MonocTime
    if strcmp('Monocular',Grating.stim(t))
        if stim1 == 0
            stim1 = stim1+1;
            stim1_ori = Grating.tilt(t);
            MonocTime1(stim1,:) = [start_noSoa stop];
        elseif stim1 > 0 && Grating.tilt(t) == stim1_ori
            stim1 = stim1+1;
            MonocTime1(stim1,:) = [start_noSoa stop];
        elseif stim1 > 0 && Grating.tilt(t) ~= stim1_ori
            stim2 = stim2+1;
            stim2_ori = Grating.tilt(t);
            MonocTime2(stim2,:) = [start_noSoa stop];            
        end
    elseif ~strcmp('Monocular', Grating.stim(t)) 
        continue % skips and dismisses all non-Monocular trials.
    else
       disp('error, please check numel(idx)==1 if statement') 
       disp(t)
    end
end

 end

onsets1 = MonocTime1(:,1)';
onsets2 = MonocTime2(:,1)';

% Load matching neural data
pre           = -50; 
post          = 200; 
chans         = []; 
spkbin1        = zeros(length(pre:post),length(onsets1),length(spikechs)); 
spksdf1        = zeros(length(pre:post),length(onsets1),length(spikechs)); 
spkbin2        = zeros(length(pre:post),length(onsets2),length(spikechs)); 
spksdf2        = zeros(length(pre:post),length(onsets2),length(spikechs)); 

%only interested in the sink. Deleted the all electrodes loop (e loop)
    clear IDX SPK 
    e = sink;
    IDX        = ppNEV.Data.Spikes.Electrode == orderedpins(e); 
    SPK        = ppNEV.Data.Spikes.TimeStamp(IDX); 
    chans(e,:) = allchanIDs{orderedpins(e)}(1:4)';  
    
    % convolve spikes 
    sdf        = spk2sdf(SPK,Fs); 
    
    % trigger spikes to events
    spkbin1(:,:,e) = trigBinaryData(floor(SPK./30),pre,post,onsets1); % binary spikes. 1== spike. 0 == no event
    spksdf1(:,:,e) = trigData(sdf',onsets1,-pre,post); 
    spkbin2(:,:,e) = trigBinaryData(floor(SPK./30),pre,post,onsets2); % binary spikes. 2== spike. 0 == no event
    spksdf2(:,:,e) = trigData(sdf',onsets2,-pre,post); 
    



%% plot
% chans -- channel label
% spkbin -- 0s and 1s with binary spike data// samples x trials x channels 
% spksdf -- convolved spike data // samples x trials x channels 

tvec = pre:post; 
% only interested in the sink
ch = sink;
e = sink;

    clear sem mn IDX 
    figure, set(gcf,'color','w','position',[1 1 660 300]); 
    
    %stim1
    sem = nanstd(spksdf1(:,:,ch),0,2)./sqrt(size(spksdf1,2)); 
    mn1  = nanmean(spksdf1(:,:,ch),2); 
    L(1) = plot(tvec,mn1,'linewidth',2,'color','b'); hold on; 
    plot(tvec,mn1 - sem,'linewidth',1,'color','b'); hold on; 
    plot(tvec,mn1 + sem,'linewidth',1,'color','b'); hold on; 


    
    % stim2
    sem = nanstd(spksdf2(:,:,ch),0,2)./sqrt(size(spksdf2,2)); 
    mn2  = nanmean(spksdf2(:,:,ch),2); 
    L(2) = plot(tvec,mn2,'linewidth',2,'color','r'); hold on; 
    plot(tvec,mn2 - sem,'linewidth',1,'color','r'); hold on; 
    plot(tvec,mn2 + sem,'linewidth',1,'color','r'); hold on; 
    

    ylabel('spks/s'); xlabel('t(ms)');
    xlim([pre post]);
    ylim([0 50])
    v = vline(0); set(v,'color','k','linestyle','-','linewidth',2); 
    set(gca,'box','off','tickdir','out','linewidth',2); 
    title(gca,BRdatafile); 
    
    legend(L, {num2str(stim1_ori), num2str(stim2_ori)})

%% ttest
stimavg(1) = mean(mn1);
stimavg(2) = mean(mn2);
[M,I] = max(stimavg);
if I == 1
    disp(strcat('pref stim equals...',num2str(stim1_ori)))
elseif I == 2
    disp(strcat('pref stim equals...',num2str(stim2_ori)))
end


[h,p] = ttest(mn1,mn2);
disp(strcat('p value is...',num2str(p)))

