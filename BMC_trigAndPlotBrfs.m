%BMC_trigAndPlotBrfs.m
%GOAL: load the previously filtered and saved brfs data, trigger, and plot
%
%   Version 1.0
%   Brock Carlson -- created 9/9/19
%   
%   Current plotting goal: plot brfs monocular fCSDs for each ori and eye
%   under scaled conditions to see if an effect exists.

%% LOAD PROCESSED DATA



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
STIM_BRFS.diop_simult_NPS.start_noSoaDown(STIM_BRFS.diop_simult_NPS.start_noSoaDown - pre  < 0) = [];
STIM_BRFS.diop_simult_NPS.start_noSoaDown(STIM_BRFS.diop_simult_NPS.start_noSoaDown + post > size(ns2LFP,2)) = [];

numTriggers.diop_simult_NPS = size(STIM_BRFS.diop_simult_NPS.start_noSoaDown,1);

%preallocate
TRIG_BRFS.diop_simult_NPS.ns2LFP         = NaN( contactNum,(post + pre + 1),numTriggers.diop_simult_NPS); 
TRIG_BRFS.diop_simult_NPS.ns2CSD         = NaN( contactNum,(post + pre + 1),numTriggers.diop_simult_NPS);


    
%%% DATA MUST BE DOWNSAMPLED TO 1KH IN ORDER TO TRIGGER PROPERLY
for singleCh = 1:contactNum 
    for singleTrigger = 1:numTriggers.diop_simult_NPS
        timeOfTrigger = STIM_BRFS.diop_simult_NPS.start_noSoaDown(singleTrigger);
        windowOfTrigger = timeOfTrigger-pre:timeOfTrigger+post;
        % output is (Ch x time x triggerNumber)
        TRIG_BRFS.diop_simult_NPS.ns2LFP(singleCh,:,singleTrigger)       = ns2LFP(singleCh,windowOfTrigger); 
        TRIG_BRFS.diop_simult_NPS.ns2CSD(singleCh,:,singleTrigger)       = ns2CSD(singleCh,windowOfTrigger);

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
AVG_BRFS.diop_simult_NPS.ns2fCSD = filterCSD(AVG_BRFS.diop_simult_NPS.ns2CSD);

% BLavg
BLavg.ns2fCSD       = filterCSD(BLavg.ns2CSD);
BLavg_BRFS.diop_simult_NPS.ns2fCSD = filterCSD(BLavg_BRFS.diop_simult_NPS.ns2CSD);



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
f_ShadedLinePlotbyDepth(BLavg_BRFS.diop_simult_NPS.ns2LFP,chans,TM,[],1)
plot([0 0], ylim,'k')
plot([800 800], ylim,'k')
title({'LFP',filename{a}}, 'Interpreter', 'none')
xlabel('time (ms)')
ylabel('Contacts indexed down from surface')

% CSD line
subplot(1,3,2)
f_ShadedLinePlotbyDepth(BLavg_BRFS.diop_simult_NPS.ns2CSD,chans,TM,[],1)
plot([0 0], ylim,'k')
plot([800 800], ylim,'k')
title({'CSD',filename{a}}, 'Interpreter', 'none')
xlabel('time (ms)')

% fCSD
subplot(1,3,3)
imagesc(TM,chans,BLavg_BRFS.diop_simult_NPS.ns2fCSD); 
colormap(flipud(jet));
climit = max(abs(get(gca,'CLim'))*.8);
set(gca,'CLim',[-climit climit],'Box','off','TickDir','out')
hold on;
plot([0 0], ylim,'k')
plot([800 800], ylim,'k')
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










