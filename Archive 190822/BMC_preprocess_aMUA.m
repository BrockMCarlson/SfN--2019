%function jnm_preprocess( jnmList, GBL, GPS, overwrite, freqs )
% BMC_preprocess.m editeded from jnm_preprocess.m 
% 190807 edited

%jnm_preprocess.m
%Jake Westerberg
%Vanderbilt University
%August 3, 2016
%v1.0.0
%   take raw continuous signal data file and parse out the individual
%   components and put them into the .jnm format.

jnmListL = length( jnmList );
jnmraw = jnm_rawpath();

for i = 1 : jnmListL
    
    jnmIdent = jnmList{ i }( 1 : ( end - 4 ) );
    jnmFile = jnmList{ i };
    jnmOut = jnm_pathout( 'preprocess' );
    jnmOutN = [ jnmOut jnmIdent '\' ];
    jnmOut = jnmOutN;
    
    if strcmp( 'nev', jnmList{ i }( end - 2 : end ) )
        
        completed = jnm_checkexist( jnmOut, [ jnmIdent '_nev.jnm' ] );
        
    elseif ~strcmp( jnmList{ i }( end - 2 : end ), 'nev' ) && ...
            ~strcmp( jnmList{ i }( end - 2 : end ), 'ns6' )
        
        completed = true;
        
    else
        
        completed = jnm_checkexist( jnmOut, [ jnmIdent '_amua_p1.jnm' ] );
        
    end
    
    if completed == true && overwrite == false
        
        continue
        
    else
        
        mkdir( jnmOut );
        
        jnmi = jnm_sessioninfo( jnmIdent( 1 : 6 ) );
        
        pNum = jnmi.pNum;
        pChan = jnmi.pCnl( 1 ); %fix at a later date
        pSpc = jnmi.pSpc( 1 ); %fix at a later date
        
        fprintf( '\nFile %u of %u: Initializing Preprocessing Sequence\n', ...
            i, jnmListL );
        
        switch jnmList{ i }( (end - 2) : end )
            
            case 'nev'
                
                if strcmp( 'rsd', jnmIdent( 10 : 12 ) )
                    
                    continue
                    
                end
                
                tempNEV = openNEV( [ jnmraw jnmFile], 'noread', 'nomat', 'nosave' );
                
                if ~isstruct( tempNEV )
                    
                    continue
                    
                end
                
                if ~isempty( tempNEV.Data.SerialDigitalIO.UnparsedData )
                    
                    nevOI( 1, : ) = ...
                        tempNEV.Data.SerialDigitalIO.UnparsedData - 128;
                    nevOI( 2, : ) = ...
                        tempNEV.Data.SerialDigitalIO.TimeStampSec * 1000;
                    
                    save( [ jnmOut jnmIdent '_nev.jnm' ], 'nevOI', ...
                        '-mat', '-v7.3');
                    
                    clear nevOI
                    
                end
                
                clear tempNEV
                
            case 'ns6'
                
                t = openNSx( [jnmraw jnmFile], 'noread' );
                
                tn = ~strcmp( 'E', { t.ElectrodesInfo.ConnectorBank } );
                nsx.E = length( tn );
                nsx.neural = sum( tn );
                nsx.anlg = sum( ~tn );
                
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
                
                if abs( nsx.rnge - 0.25 ) > .001
                    
                    warning( 'Check D2A range on File %d\n', i );
                    
                end
                
                nsx.fs = t.MetaTags.SamplingFreq;
                nsx.nyq = nsx.fs / 2;
                nsx.deci = nsx.fs / 1000;
                
                electD= openNSx( [ jnmraw jnmFile ], 'c:1', 'read' );
                ttData = double( electD.Data );
                
                samples = length( ttData );
                
                bnc = zeros( nsx.anlg, ceil( samples / nsx.deci )  );
                lfp = zeros( nsx.neural, ceil( samples / nsx.deci ) );
                amua = zeros( nsx.neural, ceil( samples / nsx.deci ) );
                dmua = zeros( nsx.neural, ceil( samples / nsx.deci ) );
                
                for nct = min( neuralChan ) : max( neuralChan )
                    
                    elect = sprintf( 'c:%u', nct );
                    electD= openNSx( [ jnmraw jnmFile ], elect, 'read' );
                    
                    tData = ( double( electD.Data ) ).' ;
                    electD.Data = [];
                    
                    tData = tData ./ 4;
                    
                    lfp( nct, : ) = jnm_lfppro( tData, ...
                        nsx.deci, 500, nsx.nyq );
                    
                    amua( nct, : ) = jnm_amuapro( tData, ...
                        nsx.deci, 500, 5000, 200, nsx.nyq );
                    
                    dmua( nct, : ) = jnm_dmuapro( tData, 2, ...
                        nsx.deci, false, false, 500, 5000, 200, nsx.nyq );
                    
                end

end


