function [TRIG_photo] = bmcTRIG_photo(STIM_photo,PARAMS,pre,post,readNS2file, ...
    readNS6file,ns2LFP,ns2CSD,ns6LFPdown,ns6CSD,aMUAdown)
%bmcTRIG.m
% BMC bersion 1.0 -- 8/28/19
%
% TRIGGER TO pEvC/pEvT at 1 kHz 
% Structure of function.
    % A) Trigger to STIM.onsetsdown
        % i) trigger ns2
        % i) trigger ns6            
    % b) Trigger to STIM_photo.onsetsdown
        % i) trigger ns2
        % i) trigger ns6
%
% outputs are TRIG and TRIG_photo with field .ns2LFP .ns2CSD .ns6LFP
% .ns6CSD and .aMUA
%
%This function can triggers data to BOTH the event codes and the
%photodiode. The photo diode is ovbiously the better signal, so that is the
%output that should be used if possible. However, often the photodiode
%didnt work or was forgotten. In this case, the pEvT codes should be used.
%This function may have to be modified to ignore the photodiode inputs in
%that case.


%% Load info I usually need
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
%% Trigger function
% % % A) Trigger to STIM.onsetsdown 
% % % remove TPs that are too close to start or end
% % STIM.onsetsdown(STIM.onsetsdown - pre  < 0) = [];
% % STIM.onsetsdown(STIM.onsetsdown + post > size(ns6LFPdown,2)) = [];
% % 
% % numTriggers = size(STIM.onsetsdown,1);
% % 
% % %preallocate
% % TRIG.ns2LFP         = NaN( contactNum,(post + pre + 1),numTriggers); 
% % TRIG.ns2CSD         = NaN( contactNum,(post + pre + 1),numTriggers);
% % TRIG.ns6LFP         = NaN( contactNum,(post + pre + 1),numTriggers);
% % TRIG.ns6CSD         = NaN( contactNum,(post + pre + 1),numTriggers);
% % TRIG.aMUA           = NaN( contactNum,(post + pre + 1),numTriggers);
% % 
% %     
% % %%% DATA MUST BE DOWNSAMPLED TO 1KH IN ORDER TO TRIGGER PROPERLY
% % for singleCh = 1:contactNum 
% %     for singleTrigger = 1:numTriggers
% %         timeOfTrigger = STIM.onsetsdown(singleTrigger);
% %         windowOfTrigger = timeOfTrigger-pre:timeOfTrigger+post;
% %         % output is (Ch x time x triggerNumber)
% %         TRIG.ns2LFP(singleCh,:,singleTrigger)       = ns2LFP(singleCh,windowOfTrigger); 
% %         TRIG.ns2CSD(singleCh,:,singleTrigger)       = ns2CSD(singleCh,windowOfTrigger);
% %         TRIG.ns6LFP(singleCh,:,singleTrigger)       = ns6LFPdown(singleCh,windowOfTrigger);
% %         TRIG.ns6CSD(singleCh,:,singleTrigger)       = ns6CSD(singleCh,windowOfTrigger);
% %         TRIG.aMUA(singleCh,:,singleTrigger)         = aMUAdown(singleCh,windowOfTrigger);
% %     end
% % end

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
TRIG_photo.ns2CSD         = NaN( contactNum,(post + pre + 1),numTriggers_photo);
TRIG_photo.ns6LFP         = NaN( contactNum,(post + pre + 1),numTriggers_photo);
TRIG_photo.ns6CSD         = NaN( contactNum,(post + pre + 1),numTriggers_photo);
TRIG_photo.aMUA           = NaN( contactNum,(post + pre + 1),numTriggers_photo);

    
%%% DATA MUST BE DOWNSAMPLED TO 1KH IN ORDER TO TRIGGER PROPERLY
for singleCh = 1:contactNum 
    for singleTrigger = 1:numTriggers_photo
        timeOfTrigger_photo = STIM_photo.onsetsdown(singleTrigger);
        windowOfTrigger_photo = timeOfTrigger_photo-pre:timeOfTrigger_photo+post;
        % output is (Ch x time x triggerNumber)
        TRIG_photo.ns2LFP(singleCh,:,singleTrigger)       = ns2LFP(singleCh,windowOfTrigger_photo); 
        TRIG_photo.ns2CSD(singleCh,:,singleTrigger)       = ns2CSD(singleCh,windowOfTrigger_photo);
        TRIG_photo.ns6LFP(singleCh,:,singleTrigger)       = ns6LFPdown(singleCh,windowOfTrigger_photo);
        TRIG_photo.ns6CSD(singleCh,:,singleTrigger)       = ns6CSD(singleCh,windowOfTrigger_photo);
        TRIG_photo.aMUA(singleCh,:,singleTrigger)         = aMUAdown(singleCh,windowOfTrigger_photo);
    end
end
end

