function [aMUAdown] = filterForaMUA(ns6DAT,readNS2file,readNS6file,PARAMS)
%BMC filterForaMUA.m
%   Version 1.0
%   Brock Carlson -- created 8/27/19
%   filter  ns6 continuous data and downsample immediatly to align with the
%   STIM.onsets (or STIM.onsetsdown ??) in 1kHz.
%
%   readNS2file only needed for NS2_header which is used to establish
%   contactLogicals etc.

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
    clear ns6FiltVecaMUA hpMUA lpMUA bwb1 bwa1 bwb2 bwa2 bwb3 bwa3
    %Downsample aMUA to 1kHz to align with triggerpoints below
    aMUAdown = downsample(aMUA',30)';
    

end

