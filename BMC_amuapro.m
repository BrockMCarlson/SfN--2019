function [ jnmSig ] = jnm_amuapro( jnmMUA, deci, hpc, lpc1, nyq )
% BMC_amuapro.m editeded from jnm_amuapro.m 
% 190807 edited

%jnm_amuapro.m
%Jake Westerberg
%Vanderbilt University
%August 4, 2016
%v1.0.0
%   Generate analog mua signal from raw signal
lpc2 = l0pc1 / 2;
hWn = hpc / nyq;
[ bwb, bwa ] = butter( 4, hWn, 'high' );
hpMUA = filtfilt( bwb, bwa, jnmMUA );
lWn = lpc1 / nyq;
[ bwb, bwa ] = butter( 4, lWn, 'low' );
hpMUA = abs( filtfilt( bwb, bwa, hpMUA ) );
if ~(deci == false)
   hpMUA = decimate( hpMUA, deci );
   nyq = nyq/(deci);
end
lWn = lpc2 / nyq;
[ bwb, bwa ] = butter( 4, lWn, 'low' );
jnmSig = filtfilt( bwb, bwa, hpMUA );
end