% editeded from Jake's analyPSDdepth_redo
%Jake Westerberg
%Vanderbilt University
%August 19, 2016
%v1.0.0
%   Detailed explanation goes here


tic
clear

addpath(genpath('G:\LaCie\all BRFS'));
drname        = {'G:\LaCie\all BRFS\160102_E\'};%%,'G:\LaCie\all BRFS\160427_E\',...
BRdatafile    = {'160102_E_rfori001'};%%'160427_E_brfs001','160510_E_brfs001'};
exportfigtext = {'rfori_160102'};%%'brfs_60427','brfs_160510'};

a = 1;
cd(drname{a})
jnmfile = [BRdatafile{a} '.ns2']; 

if ~exist( 'val', 'var' )    
    val = 15; 
end

t = openNSx( jnmfile, 'noread' );

tn = ~strcmp( 'E', { t.ElectrodesInfo.ConnectorBank } );
nsx.E = length( tn );
nsx.neural = sum( tn );
nsx.anlg = sum( ~tn );

if exist('electrode')
error('what goes here?')   
end

neuralChan = find( tn );
bncChan = find( ~tn );

nsx.neuralL = { t.ElectrodesInfo( tn ).Label };
nsx.neuralI = t.ElectrodesInfo( tn );

nsx.bncL = { t.ElectrodesInfo( ~tn ).Label };
nsx.bncI = t.ElectrodesInfo( ~tn );

nsx.rnge = ...
    ( ( length( t.ElectrodesInfo( 1 ).MinAnalogValue : ...
    t.ElectrodesInfo( 1 ).MaxAnalogValue ) ) ./ ...
    ( length( t.ElectrodesInfo( 1 ).MinDigiValue : ...
    t.ElectrodesInfo( 1 ).MaxDigiValue ) ) );

nsx.fs = t.MetaTags.SamplingFreq;
nsx.nyq = nsx.fs / 2;
nsx.deci = nsx.fs / 1000;

electD = openNSx( jnmfile, 'c:1', 'read' );
tData = double( electD.Data );

samples = length( tData );

bnc = zeros( nsx.anlg, ceil( samples / nsx.deci )  );
lfp = zeros( nsx.neural, ceil( samples / nsx.deci ) );
    
electD = openNSx( jnmfile );

tData = ( double( electD.Data ) ).' ;
t2Data = tData( :, bncChan );
tData = tData( :, neuralChan );

electD.Data = [];

tData = tData ./ 4;

lfp = tData.';
bnc = t2Data.';

clear tData t2Data electD 

if size( lfp, 1 ) < 33
    
    pNum = 1;
    
else
    
    pNum = 2;
    
end


if pNum == 2
    pChan(1) = 24;
    idx(1,:) = [1 24];
    pChan(2) = 24;
    idx(2,:) = [25 48];
else
    
pChan = size( lfp, 1 ) / pNum;
idx(1,:) = [1 pChan];
end


for k = 1 : pNum
    
lfp2 = jnm_reorder( lfp(idx(k,1):idx(k,2),:), nsx.neuralL(idx(k,1):idx(k,2)), 'BR', pNum, pChan(k) );

valn = 2^val;

pdrplot = figure( 'Position', [ 0 0 750 750 ] );
    
    jnm = zeros( pChan(k), 257 );
    
    for i = 1:pChan(k)
        
        [ jnm( i, : ) ] = jnm_psd( lfp2( i, end - valn : end ), ...
            512, 1000, 512, 0);
        
    end
    
    jnm2 = zeros( pChan(k), 52 );
    
    for i = 1 : 52
        
        for j = 1 : pChan(k)
            
            jnm2( j, i ) = ( jnm( j, i ) - mean( jnm( :, i ) ) ) ...
                / mean( jnm( :, i ) ) * 100;
            
        end
    end
    
    image( jnm2 );
    colormap('hot');
    set(gca,  'XTickLabel', [] );
    
end

toc
[~,fname] = fileparts(BRdatafile);
title(fname,'interpreter','none')
