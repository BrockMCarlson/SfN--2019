function [ jnmNew ] = jnm_reorder( jnm, jnmOrder, oType, nElect, nChan )
%jnm_reorder.m
%Jake Westerberg
%Vanderbilt University
%August 2, 2016
%v1.0.0
%   Reorder jnm with the vector jnmOrder

if strcmp( oType, 'BR' )
    
    for i = 1:length( jnmOrder )
        
        if strcmp( jnmOrder{ i }( 1 ), 'b' )
            
            continue
            
        end
        
        tNum = str2num( jnmOrder{ i }( 3 : 4 ) );
        jnmNewO( i ) = tNum;
        
    end
    
    if nElect == 2
        
        jnmNewO( ( nChan + 1 ) : end ) = ...
            jnmNewO( ( nChan + 1 ) : end ) + nChan;
        
    end
    
    for i = 1 : length( jnmNewO )
        
        jnmNew( i, : ) = jnm( find( jnmNewO == i ), : );
        
    end
    
else
    
    disp('This function does not support reordering of this data type.')
        
end

end

