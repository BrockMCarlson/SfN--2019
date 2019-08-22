%% BMC - find Preffered Stim in brfs from discretized MUA spiking data
% combination of BMC conditions (looking at monocular) and KD ppnev
% analysis.

% Update 6/20/2019

clear 
close all

pre           = 50; 
post          = 150; 

brdrname      = 'G:\LaCie\all BRFS_recentB-52\190410_B';
cd(brdrname)
BRdatafile    = '190410_E_brfs001'; 
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


% load ppNEV
load(strcat(filename,'.ppnev'),'-MAT');
allchanIDs  = {ppNEV.ElectrodesInfo.ElectrodeLabel};  % every pin on entire 128 channel system 
spikechs    = nanunique(ppNEV.Data.Spikes.Electrode); % every pin with spike data 

% load grating info
cd('G:\LaCie\all BRFS\160102_E')
fileForGrating = strcat(BRdatafile,'.gBrfsGratings');
Grating = readBRFS(fileForGrating);


            

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


%% AVERAGE AROUND SINK LOOP
% Load matching neural data
sinkRng = sink-2:1:sink+1;
chans         = []; 
spksdf.A.dat   = zeros(length(-pre:post),length(onsets1),length(sinkRng)); 
spksdf.B.dat   = zeros(length(-pre:post),length(onsets2),length(sinkRng)); 

for e = 1:length(sinkRng)
    clear IDX SPK 
   
    IDX        = ppNEV.Data.Spikes.Electrode == orderedpins(sinkRng(e)); 
    SPK        = ppNEV.Data.Spikes.TimeStamp(IDX); 
    chans(e,:) = allchanIDs{orderedpins(sinkRng(e))}(1:4)';  
    
    % convolve spikes 
    sdf        = spk2sdf(SPK,Fs); 
    
    % trigger spikes to events
    spksdf.A.dat(:,:,e) = trigData(sdf',onsets1,pre,post); 
    spksdf.B.dat(:,:,e) = trigData(sdf',onsets2,pre,post); 
    
end

spksdf.A.avg = squeeze(mean(spksdf.A.dat,3));
spksdf.B.avg = squeeze(mean(spksdf.B.dat,3));

%plot
tvec = -pre:post; 
clear sem mn IDX 
    figure, set(gcf,'color','w','position',[1 1 660 300]); 
    
    %stim1
    sem = nanstd(spksdf.A.avg(:,:),0,2)./sqrt(size(spksdf.A.avg,2)); 
    mn1  = nanmean(spksdf.A.avg(:,:),2); %samples x trials, avgbyTrials
    L(1) = plot(tvec,mn1,'linewidth',2,'color','b'); hold on; 
    plot(tvec,mn1 - sem,'linewidth',1,'color','b'); hold on; 
    plot(tvec,mn1 + sem,'linewidth',1,'color','b'); hold on; 


    
    % stim2
    sem = nanstd(spksdf.B.avg(:,:),0,2)./sqrt(size(spksdf.B.avg,2)); 
    mn2  = nanmean(spksdf.B.avg(:,:),2); 
    L(2) = plot(tvec,mn2,'linewidth',2,'color','r'); hold on; 
    plot(tvec,mn2 - sem,'linewidth',1,'color','r'); hold on; 
    plot(tvec,mn2 + sem,'linewidth',1,'color','r'); hold on; 
    

    ylabel('spks/s'); xlabel('t(ms)');
    xlim([-pre post]);
%     ylim([0 50])
    v = vline(0); set(v,'color','k','linestyle','-','linewidth',2); 
    set(gca,'box','off','tickdir','out','linewidth',2); 
    
    legend(L, {num2str(stim1_ori), num2str(stim2_ori)})

    % ttest
    stimavg(1) = mean(mn1); %avg across all samples from avg from all trials
    stimavg(2) = mean(mn2);
    [M,I] = max(stimavg);
    if I == 1
        pforitxt = strcat('Pref stim=',num2str(stim1_ori));
    elseif I == 2
        pforitxt = strcat('Pref stim=',num2str(stim2_ori));
    end

    [h,p] = ttest(mn1,mn2);

    % plot title
    titletext = {strcat('Avg for chans--',num2str(sinkRng)),...
        strcat(BRdatafile,'.',pforitxt,'.','p value=',num2str(p))};
    title(gca,titletext,'interpreter','none'); 
    
    %export granular section
    cd('G:\LaCie\SfN 2019\SfN 2019 figs\MUA_PS')
    export_fig(strcat(BRdatafile,'_Granular_v2'),'-pdf','-nocrop') 


%% FULL ELECTRODE LOOP

% Load matching neural data

chans         = []; 
spksdf1        = zeros(length(-pre:post),length(onsets1),length(spikechs)); 
spksdf2        = zeros(length(-pre:post),length(onsets2),length(spikechs)); 



for e = 1:length(orderedpins)
    clear IDX SPK 
   
    IDX        = ppNEV.Data.Spikes.Electrode == orderedpins(e); 
    SPK        = ppNEV.Data.Spikes.TimeStamp(IDX); 
    chans(e,:) = allchanIDs{orderedpins(e)}(1:4)';  
    
    % convolve spikes 
    sdf        = spk2sdf(SPK,Fs); 
    
    % trigger spikes to events
    spksdf1(:,:,e) = trigData(sdf',onsets1,pre,post); 
    spksdf2(:,:,e) = trigData(sdf',onsets2,pre,post); 
    
end


% plot
% chans -- channel label
% spkbin -- 0s and 1s with binary spike data// samples x trials x channels 
% spksdf -- convolved spike data // samples x trials x channels 

tvec = -pre:post; 
PSRESULT = {'Chan','Sig?','PVal','PS'};


for ch = 1:size(spksdf1,3)
%%PLOT
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
    xlim([-pre post]);
%     ylim([0 50])
    v = vline(0); set(v,'color','k','linestyle','-','linewidth',2); 
    set(gca,'box','off','tickdir','out','linewidth',2); 
    
    legend(L, {num2str(stim1_ori), num2str(stim2_ori)})

%% ttest
stimavg(1) = mean(mn1);
stimavg(2) = mean(mn2);
[M,I] = max(stimavg);
if I == 1
    pforitxt = strcat('Pref stim=',num2str(stim1_ori));
elseif I == 2
    pforitxt = strcat('pref stim=',num2str(stim2_ori));
end


[h,p] = ttest(mn1,mn2);


%% plot title
titletext = {strcat('Channel --',num2str(ch)),...
        strcat(BRdatafile,'.',pforitxt,'.','p value=',num2str(p))};
title(gca,titletext,'interpreter','none'); 

%% Append a pdf
cd('G:\LaCie\SfN 2019\SfN 2019 figs\MUA_PS')
if ch ==1
    export_fig(strcat(BRdatafile,'_allCh_v2'),'-pdf','-nocrop') 
else
    export_fig(strcat(BRdatafile,'_allCh_v2'),'-pdf','-nocrop','-append')
end

%% Make variable for full electrode for PS vs NPS for day

PSRESULT{ch+1,1} = ch;

if p <.05
    PSRESULT{ch+1,2} = 1;
else
    PSRESULT{ch+1,2} = 0;
end

PSRESULT{ch+1,3} = p;

if I == 1
    PSRESULT{ch+1,4} = stim1_ori;
elseif I == 2
     PSRESULT{ch+1,4} = stim2_ori;
end





end

save(strcat(BRdatafile,'_PSRESULT'),'PSRESULT')