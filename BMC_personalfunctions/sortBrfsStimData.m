function [STIM_BRFS] = sortBrfsStimData(STIM)
%UNTITLED2 
%   takes STIM and gets the brfs info for each of the desired
%conditions

conditioncount = 0; 
clear gratingOnsets
switch condition{c}
    
 case 'biPSNOsoa'
    for  e = 1:length(pEvC)
     if strcmp('Binocular',grating.stim(e))           && ... 
            grating.soa(e)          ==  0       && ...
            grating.s1_contrast(e)  >= .5     && ...
            grating.s2_contrast(e)  >= .5     && ...
            grating.s1_tilt(e)      ==  SessionParamsForCondition.PS(a)  && ...         
            grating.s2_tilt(e)      ==  SessionParamsForCondition.PS(a)        
        if ~any(pEvC{e} == 96) 
             continue
        end     
        stimon   =  pEvC{e} == 23;
        stimoff  =  pEvC{e} == 24;
        idx = find(stimon);
        if      numel(idx) == 1     % there is no soa.
            start_noSoa  =  pEvT{e}(idx); 
        else
            disp('error, please check idx loop')
        end
        stop    =  pEvT{e}(stimoff);
        conditioncount = conditioncount +1;
        gratingOnsets(conditioncount,:) = [start_noSoa stop];
     end
    column = 1;  
    end
     
 case 'biNPSNOsoa'
    for  e = 1:length(pEvC)
     if strcmp('Binocular',grating.stim(e))           && ... 
            grating.soa(e)          ==  0       && ...
            grating.s1_contrast(e)  >= .5     && ...
            grating.s2_contrast(e)  >= .5     && ...
            grating.s1_tilt(e)      ==  SessionParamsForCondition.NPS(a)  && ...         
            grating.s2_tilt(e)      ==  SessionParamsForCondition.NPS(a)        
        if ~any(pEvC{e} == 96) 
             continue
        end     
        stimon   =  pEvC{e} == 23;
        stimoff  =  pEvC{e} == 24;
        idx = find(stimon);
        if      numel(idx) == 1     % there is no soa.
            start_noSoa  =  pEvT{e}(idx); 
        else
            disp('error, please check idx loop')
        end
        stop    =  pEvT{e}(stimoff);
        conditioncount = conditioncount +1;
        gratingOnsets(conditioncount,:) = [start_noSoa stop];
     end
    column = 1;  
    end
   if ~exist('gratingOnsets')
      disp(strcat('Binocular NPS does not exist for a=',string(a)))
       continue
   end
   
 case 'dicopNOsoa'
   for  e = 1:length(pEvC)
     if strcmp('dCOS',grating.stim(e))           && ... 
            grating.soa(e)          ==  0       && ...
            grating.s1_contrast(e)  >= .5     && ...
            grating.s2_contrast(e)  >=  .5
        if ~any(pEvC{e} == 96) 
             continue
        end     
        stimon   =  pEvC{e} == 23;
        stimoff  =  pEvC{e} == 24;
        idx = find(stimon);
        if      numel(idx) == 1     % there is no soa.
            start_noSoa  =  pEvT{e}(idx); 
        else
            disp('error, please check idx loop')
        end
        stop    =  pEvT{e}(stimoff);
        conditioncount = conditioncount +1;
        gratingOnsets(conditioncount,:) = [start_noSoa stop];
     end
    column = 1;  
   end 
   
 case 'biPSWsoa' 
    for  e = 1:length(pEvC)
     if strcmp('Binocular',grating.stim(e))           && ... 
            grating.soa(e)          ==  800       && ...
            grating.s1_contrast(e)  >= .5     && ...
            grating.s2_contrast(e)  >= .5     && ...
            grating.s1_tilt(e)      ==  SessionParamsForCondition.PS(a)  && ...         
            grating.s2_tilt(e)      ==  SessionParamsForCondition.PS(a)        
        if ~any(pEvC{e} == 96) 
             continue
        end     
        stimon   =  pEvC{e} == 23;
        stimoff  =  pEvC{e} == 24;
        idx = find(stimon);
        if	numel(idx) == 2     %there is indeed soa
            start1  =  pEvT{e}(idx(1));
            start2  =  pEvT{e}(idx(2));
        else
            disp('error, please check idx loop')
        end
        stop    =  pEvT{e}(stimoff);
        conditioncount = conditioncount +1;
        gratingOnsets(conditioncount,:) = [start1 start2 stop];
    end
    column = 1;  
    end
   
     case 'biNPSWsoa' 
    for  e = 1:length(pEvC)
     if strcmp('Binocular',grating.stim(e))           && ... 
            grating.soa(e)          ==  800       && ...
            grating.s1_contrast(e)  >= .5     && ...
            grating.s2_contrast(e)  >= .5     && ...
            grating.s1_tilt(e)      ==  SessionParamsForCondition.NPS(a)  && ...         
            grating.s2_tilt(e)      ==  SessionParamsForCondition.NPS(a)        
        if ~any(pEvC{e} == 96) 
             continue
        end     
        stimon   =  pEvC{e} == 23;
        stimoff  =  pEvC{e} == 24;
        idx = find(stimon);
        if	numel(idx) == 2     %there is indeed soa
            start1  =  pEvT{e}(idx(1));
            start2  =  pEvT{e}(idx(2));
        else
            disp('error, please check idx loop')
        end
        stop    =  pEvT{e}(stimoff);
        conditioncount = conditioncount +1;
        gratingOnsets(conditioncount,:) = [start1 start2 stop];
    end
    column = 1;  
   end
    
 case 'dichopWsoa_fullTrialPS'
   for  e = 1:length(pEvC)
    if strcmp('dCOS',grating.stim(e))           && ... 
            grating.soa(e)          ==  800     && ...
            grating.s1_contrast(e)  >= .5     && ...
            grating.s2_contrast(e)  >= .5     && ...
            grating.s2_tilt(e)      ==  SessionParamsForCondition.PS(a)
        if ~any(pEvC{e} == 96) 
             continue
        end     
        stimon   =  pEvC{e} == 23;
        stimoff  =  pEvC{e} == 24;
        idx = find(stimon);
        if	numel(idx) == 2     %there is indeed soa
            start1  =  pEvT{e}(idx(1));
            start2  =  pEvT{e}(idx(2));
        else
            disp('error, please check idx loop')
        end
        stop    =  pEvT{e}(stimoff);
        conditioncount = conditioncount +1;
        gratingOnsets(conditioncount,:) = [start1 start2 stop];
    end
    column = 1;  
  
   end
    
 case 'dicopWsoa_PSflash'
   for  e = 1:length(pEvC)
    if strcmp('dCOS',grating.stim(e))           && ... 
            grating.soa(e)          ==  800     && ...
            grating.s1_contrast(e)  >= .5     && ...
            grating.s2_contrast(e)  >= .5     && ...
            grating.s2_tilt(e)      ==  SessionParamsForCondition.PS(a)
        if ~any(pEvC{e} == 96) 
             continue
        end     
        stimon   =  pEvC{e} == 23;
        stimoff  =  pEvC{e} == 24;
        idx = find(stimon);
        if	numel(idx) == 2     %there is indeed soa
            start1  =  pEvT{e}(idx(1));
            start2  =  pEvT{e}(idx(2));
        else
            disp('error, please check idx loop')
        end
        stop    =  pEvT{e}(stimoff);
        conditioncount = conditioncount +1;
        gratingOnsets(conditioncount,:) = [start1 start2 stop];
    end
    column = 2;  
   end
   
    case 'dicopWsoa_NPSflash'
   for  e = 1:length(pEvC)
    if strcmp('dCOS',grating.stim(e))           && ... 
            grating.soa(e)          ==  800     && ...
            grating.s1_contrast(e)  >= .5     && ...
            grating.s2_contrast(e)  >= .5     && ...
            grating.s2_tilt(e)      ==  SessionParamsForCondition.NPS(a)
        if ~any(pEvC{e} == 96) 
             continue
        end     
        stimon   =  pEvC{e} == 23;
        stimoff  =  pEvC{e} == 24;
        idx = find(stimon);
        if	numel(idx) == 2     %there is indeed soa
            start1  =  pEvT{e}(idx(1));
            start2  =  pEvT{e}(idx(2));
        else
            disp('error, please check idx loop')
        end
        stop    =  pEvT{e}(stimoff);
        conditioncount = conditioncount +1;
        gratingOnsets(conditioncount,:) = [start1 start2 stop];
    end
    column = 2;  
   end
   
end



end

