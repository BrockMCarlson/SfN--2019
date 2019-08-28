function [ns2LFP,ns6LFPdown] = filterForLFP(ns2DAT,ns6DAT,readNS2file,readNS6file,PARAMS)
%BMC filterForLFP.m
%   Version 1.0
%   Brock Carlson -- created 8/27/19
%   filter  LFP and downsample the ns6LFP immediatly to align with 
%   STIM.onsets in 1kHz.

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

%% Custom filtering ns2 and ns65 for lfp, contact by contact.
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
    clear ns2FiltVec ns6FiltVec  butter_b butter_a
    %Downsample ns6LFP to 1kHz to align with triggerpoints below
    
    ns6LFPdown = downsample(ns6LFP',30)';
    
end

