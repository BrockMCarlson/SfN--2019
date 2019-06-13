% example code for ppNEV files
clear 

brdrname      = 'E:\LaCie\BRFSdays ppnev rfori and brfs'; 
BRdatafile    = '160104_E_rfori002'; 
filename      = [brdrname '/' BRdatafile]; 

% load ppNEV
load(strcat(filename,'.ppnev'),'-MAT');
allchanIDs  = {ppNEV.ElectrodesInfo.ElectrodeLabel};  % every pin on entire 128 channel system 
spikechs    = nanunique(ppNEV.Data.Spikes.Electrode); % every pin with spike data 

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
onsets        = EventTimes(EventCodes == 23 | EventCodes == 25 | EventCodes == 27 | EventCodes == 29 | EventCodes == 31); 
pre           = -50; 
post          = 200; 
chans         = []; 
spkbin        = zeros(length(pre:post),length(onsets),length(spikechs)); 
spksdf        = zeros(length(pre:post),length(onsets),length(spikechs)); 
for e = 1:length(orderedpins)
    clear IDX SPK 
    
    IDX        = ppNEV.Data.Spikes.Electrode == orderedpins(e); 
    SPK        = ppNEV.Data.Spikes.TimeStamp(IDX); 
    chans(e,:) = allchanIDs{orderedpins(e)}(1:4)';  %#ok<SAGROW>
    
    % convolve spikes 
    sdf        = spk2sdf(SPK,Fs); 
    
    % trigger spikes to events
    spkbin(:,:,e) = trigBinaryData(floor(SPK./30),pre,post,onsets); % binary spikes. 1== spike. 0 == no event
    spksdf(:,:,e) = trigData(sdf',onsets,-pre,post); 
    

end


% chans -- channel label
% spkbin -- 0s and 1s with binary spike data// samples x trials x channels 
% spksdf -- convolved spike data // samples x trials x channels 

tvec = pre:post; 
for ch = 1:size(spksdf,3)
    clear sem mn IDX 
    figure, set(gcf,'color','w','position',[1 1 660 300]); 
    subplot(1,2,1)
    sem = nanstd(spksdf(:,:,ch),0,2)./sqrt(size(spksdf,2)); 
    mn  = nanmean(spksdf(:,:,ch),2); 
    plot(tvec,mn,'linewidth',2,'color','b'); hold on; 
    plot(tvec,mn - sem,'linewidth',1,'color','b'); hold on; 
    plot(tvec,mn + sem,'linewidth',1,'color','b'); hold on; 
    
    ylabel('spks/s'); xlabel('t(ms)'); 
    v = vline(0); set(v,'color','k','linestyle','-','linewidth',2); 
    set(gca,'box','off','tickdir','out','linewidth',2); 
    title(gca,sprintf('%s',chans(ch,:))); 
    xlim([pre post]);
    
    subplot(1,2,2)
    IDX = spkbin(:,:,e); 
    imagesc(tvec,[],IDX'); 
    colormap(flipud(colormap('gray'))); 
    set(gca,'tickdir','out','box','off','linewidth',2);
    v = vline(0); set(v,'color','k','linestyle','-','linewidth',2); 
    ylabel('trial'); xlabel('t(ms)'); 
    xlim([pre post]); 

end
