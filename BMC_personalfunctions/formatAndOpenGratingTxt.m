function [grating,readGRATINGfile] = formatAndOpenGratingTxt(directory,filename)
%BMC formatAndOpenGratingTxt.m
%   Version 1.0
%   Brock Carlson -- created 8/27/19
%   This is largely taken from Kacie's bootcamp that she held before she
%   left the lab. While certainly not necessary to use, it can be
%   frusturating to format the inputs into the various readGrating.m
%   functions


% grating file name
patterns   = {'rforidrft','rfsfdrft','posdisparitydrft','disparitydrft','cinterocdrft','coneinterocdrft','conedrft', ...
                'colorflicker','bwflicker','rfori','rfsize','cinteroc','color','rfsf','mcosinteroc','dotmapping','brfs'}; 
for p = 1:length(patterns)
   pattern      = patterns{p}; 
   if any(strfind(filename,pattern))
       startlog = strfind(filename,pattern); 
       if ~isequal(filename(startlog:end-3),pattern)
            continue
       else
            match    = patterns{p}; 
       end
   end   
end
if isequal(match,'dotmapping')
    gratingext  = '.gDotsXY_di';
elseif isequal(match,'brfs')
    gratingext = '.gBrfsGratings';
else
    gratingext  = ['.g' upper(match) 'Grating_di']; 
end
readGRATINGfile = strcat(directory,filesep,filename,gratingext);


%% Load text file
if contains(gratingext,'DRFT')
      grating     = readgDRFTGrating([filename gratingext]); % from nbanalysis 
elseif contains(gratingext,'Dots')
      grating     = readgDotsXY([filename gratingext]);
elseif contains(gratingext,'Brfs')
      grating = readBRFS([filename gratingext]);
else
      grating     = readgGrating([filename gratingext]);
end
end

