function [STIM_BRFS] = sortBrfsStimandTimeData(grating,pEvC,pEvT,PARAMS)
%sortBrfsStimData
%   from KD's sortStimandTimeData
%   allows for dual stim onsets
%   ONLY works for brfs files


if ~isfield(grating,'pathw')
    if isfield(grating,'path')
        grating.pathw = grating.path;
    elseif isfield(grating,'s1_contrast')
        grating.pathw = grating.s1_contrast;
    else
        grating.pathw = grating.contrast;
    end
end


clear spkTPs STIM
if any(strfind(grating.filename,'brfs'))
    stimfeatures = {...
        'trial'...
        'xpos'...
        'ypos'...
        'tilt'...
        'sf'...
        'diameter'...
        'eye'...
        'eye_contrast'...
        'other_contrast'...
        'oridist'...
        'gabor'...
        'gabor_std'...
        'total_stim_dur'...
        'isi_dur'...
        's1_eye'...
        's2_eye'...
        's1_tilt'       ...
        's2_tilt'       ...
        's1_contrast'   ...
        's2_contrast'   ...
        'soa'           ...
        'stim'          ...
        'pres'          ...
        };
elseif ~any(strfind(grating.filename,'brfs'))
    error('this function is for brfs use only')
end

if isfield(grating,'fix_x')
    stimfeatures{end+1} = 'fix_x';
    stimfeatures{end+1} = 'fix_y';
end





clear STIM_BRFS
%% Monocular PS
%%%%% ADD
% Monocular NPS
    %%%%% ADD

 %'diop_simult_PS'
     obs = 0; clear e codes stim times
    for  e = 1:length(pEvC)
    codes  = pEvC{e};
    stim   = find(grating.trial == e);
    times  = pEvT{e}; 
         if strcmp('Binocular',grating.stim(e))           && ... 
                grating.soa(e)          ==  0       && ...
                grating.s1_contrast(e)  >= .5     && ...
                grating.s2_contrast(e)  >= .5     && ...
                grating.s1_tilt(e)      ==  PARAMS.PS  && ...         
                grating.s2_tilt(e)      ==  PARAMS.PS       
            if ~any(pEvC{e} == 96) % I could add something here to pull out other onsets if needed.
                 continue
            end     
            stimon   =  pEvC{e} == 23;
            stimoff  =  pEvC{e} == 24;
            idx = find(stimon);
            if      numel(idx) == 1     % there is no soa.
                start_noSoa  =  pEvT{e}(idx); 
            else
                disp('error, please check diop_simult_PS loop')
            end
            stop    =  pEvT{e}(stimoff);
            obs = obs +1;
            STIM_BRFS.diop_simult_PS.start_noSoa(obs,:) = start_noSoa;
            STIM_BRFS.diop_simult_PS.stop(obs,:) = stop;
            STIM_BRFS.diop_simult_PS.start_noSoaDown(obs,:) = start_noSoa./30;
            STIM_BRFS.diop_simult_PS.stopDown(obs,:) = stop./30;            
            STIM_BRFS.diop_simult_PS.trstart(obs,:)  = double(times(find(codes == 9,1,'first')));
            STIM_BRFS.diop_simult_PS.trend(obs,:)    = double(times(find(codes == 18,1,'first')));
            STIM_BRFS.diop_simult_PS.obs(obs,:)  = obs; 
            for f = 1:length(stimfeatures)
                 STIM_BRFS.diop_simult_NPS.(stimfeatures{f})(obs,:) = grating.(stimfeatures{f})(stim);
            end
         end

    end
     
 % 'diop_simult_NPS'
     obs = 0; clear e codes stim times
    for  e = 1:length(pEvC)
    codes  = pEvC{e};
    stim   = find(grating.trial == e);
    times  = pEvT{e};     
         if strcmp('Binocular',grating.stim(e))           && ... 
                grating.soa(e)          ==  0       && ...
                grating.s1_contrast(e)  >= .5     && ...
                grating.s2_contrast(e)  >= .5     && ...
                grating.s1_tilt(e)      ==  PARAMS.NPS  && ...         
                grating.s2_tilt(e)      ==  PARAMS.NPS        
            if ~any(pEvC{e} == 96) 
                 continue
            end     
            stimon   =  pEvC{e} == 23;
            stimoff  =  pEvC{e} == 24;
            idx = find(stimon);
            if      numel(idx) == 1     % there is no soa.
                start_noSoa  =  pEvT{e}(idx); 
            else
                disp('error, please check diop_simult_NPS loop')
            end
            stop    =  pEvT{e}(stimoff);
            obs = obs +1;
            STIM_BRFS.diop_simult_NPS.start_noSoa(obs,:) = start_noSoa;
            STIM_BRFS.diop_simult_NPS.stop(obs,:) = stop;
            STIM_BRFS.diop_simult_NPS.start_noSoaDown(obs,:) = start_noSoa./30;
            STIM_BRFS.diop_simult_NPS.stopDown(obs,:) = stop./30; 
            STIM_BRFS.diop_simult_NPS.trstart(obs,:)  = double(times(find(codes == 9,1,'first')));
            STIM_BRFS.diop_simult_NPS.trend(obs,:)    = double(times(find(codes == 18,1,'first')));
            STIM_BRFS.diop_simult_NPS.obs(obs,:)  = obs; 
            for f = 1:length(stimfeatures)
                 STIM_BRFS.diop_simult_NPS.(stimfeatures{f})(obs,:) = grating.(stimfeatures{f})(stim);
            end
        end
    end
    
    
   
 % 'dichop_simult'
     obs = 0; clear e codes stim times
   for  e = 1:length(pEvC)
    codes  = pEvC{e};
    stim   = find(grating.trial == e);
    times  = pEvT{e};      
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
                disp('error, please check dichop_simult loop')
            end
            stop    =  pEvT{e}(stimoff);
            obs = obs +1;
            STIM_BRFS.dichop_simult.start_noSoa(obs,:) = start_noSoa;
            STIM_BRFS.dichop_simult.stop(obs,:) = stop;
            STIM_BRFS.dichop_simult.start_noSoaDown(obs,:) = start_noSoa./30;
            STIM_BRFS.dichop_simult.stopDown(obs,:) = stop./30; 
            STIM_BRFS.dichop_simult.trstart(obs,:)  = double(times(find(codes == 9,1,'first')));
            STIM_BRFS.dichop_simult.trend(obs,:)    = double(times(find(codes == 18,1,'first')));
            STIM_BRFS.dichop_simult.obs(obs,:)  = obs; 
            for f = 1:length(stimfeatures)
                 STIM_BRFS.dichop_simult.(stimfeatures{f})(obs,:) = grating.(stimfeatures{f})(stim);
            end
         end
   end 
   
 % 'diop_800soa_PS'
     obs = 0; clear e codes stim times
    for  e = 1:length(pEvC)
    codes  = pEvC{e};
    stim   = find(grating.trial == e);
    times  = pEvT{e};
        if strcmp('Binocular',grating.stim(e))           && ... 
                grating.soa(e)          ==  800       && ...
                grating.s1_contrast(e)  >= .5     && ...
                grating.s2_contrast(e)  >= .5     && ...
                grating.s1_tilt(e)      ==  PARAMS.PS  && ...         
                grating.s2_tilt(e)      ==  PARAMS.PS        
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
                disp('error, please check diop_800soa_PS loop')
            end
            stop    =  pEvT{e}(stimoff);
            obs = obs +1;
            STIM_BRFS.diop_800soa_PS.start1(obs,:) = start1;
            STIM_BRFS.diop_800soa_PS.start2(obs,:) = start2;
            STIM_BRFS.diop_800soa_PS.stop(obs,:) = stop;
            STIM_BRFS.diop_800soa_PS.start1Down(obs,:) = start1./30;
            STIM_BRFS.diop_800soa_PS.start2Down(obs,:) = start2./30;
            STIM_BRFS.diop_800soa_PS.stopDown(obs,:) = stop./30;
            STIM_BRFS.diop_800soa_PS.trstart(obs,:)  = double(times(find(codes == 9,1,'first')));
            STIM_BRFS.diop_800soa_PS.trend(obs,:)    = double(times(find(codes == 18,1,'first')));
            STIM_BRFS.diop_800soa_PS.obs(obs,:)  = obs; 
            for f = 1:length(stimfeatures)
                 STIM_BRFS.diop_800soa_PS.(stimfeatures{f})(obs,:) = grating.(stimfeatures{f})(stim);
            end
        end
    end
   
     % 'diop_800soa_NPS' 
         obs = 0; clear e codes stim times
    for  e = 1:length(pEvC)
    codes  = pEvC{e};
    stim   = find(grating.trial == e);
    times  = pEvT{e};   
         if strcmp('Binocular',grating.stim(e))           && ... 
                grating.soa(e)          ==  800       && ...
                grating.s1_contrast(e)  >= .5     && ...
                grating.s2_contrast(e)  >= .5     && ...
                grating.s1_tilt(e)      ==  PARAMS.NPS  && ...         
                grating.s2_tilt(e)      ==  PARAMS.NPS        
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
                disp('error, please check diop_800soa_NPS loop')
            end
            stop    =  pEvT{e}(stimoff);
            obs = obs +1;
            STIM_BRFS.diop_800soa_NPS.start1(obs,:) = start1;
            STIM_BRFS.diop_800soa_NPS.start2(obs,:) = start2;
            STIM_BRFS.diop_800soa_NPS.stop(obs,:) = stop;
            STIM_BRFS.diop_800soa_NPS.start1Down(obs,:) = start1./30;
            STIM_BRFS.diop_800soa_NPS.start2Down(obs,:) = start2./30;
            STIM_BRFS.diop_800soa_NPS.stopDown(obs,:) = stop./30;
            STIM_BRFS.diop_800soa_NPS.trstart(obs,:)  = double(times(find(codes == 9,1,'first')));
            STIM_BRFS.diop_800soa_NPS.trend(obs,:)    = double(times(find(codes == 18,1,'first')));
            STIM_BRFS.diop_800soa_NPS.obs(obs,:)  = obs; 
            for f = 1:length(stimfeatures)
                 STIM_BRFS.diop_800soa_NPS.(stimfeatures{f})(obs,:) = grating.(stimfeatures{f})(stim);
            end
         end
   end
    
 %  'dichop_800soa_brfsPSflash'
     obs = 0; clear e codes stim times
   for  e = 1:length(pEvC)
    codes  = pEvC{e};
    stim   = find(grating.trial == e);
    times  = pEvT{e};
        if strcmp('dCOS',grating.stim(e))           && ... 
                grating.soa(e)          ==  800     && ...
                grating.s1_contrast(e)  >= .5     && ...
                grating.s2_contrast(e)  >= .5     && ...
                grating.s2_tilt(e)      ==  PARAMS.PS
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
                disp('error, please check dichop_800soa_brfsPSflash loop')
            end
            stop    =  pEvT{e}(stimoff);
            obs = obs +1;
            STIM_BRFS.dichop_800soa_brfsPSflash.start1(obs,:) = start1;
            STIM_BRFS.dichop_800soa_brfsPSflash.start2(obs,:) = start2;
            STIM_BRFS.dichop_800soa_brfsPSflash.stop(obs,:) = stop;
            STIM_BRFS.dichop_800soa_brfsPSflash.start1Down(obs,:) = start1./30;
            STIM_BRFS.dichop_800soa_brfsPSflash.start2Down(obs,:) = start2./30;
            STIM_BRFS.dichop_800soa_brfsPSflash.stopDown(obs,:) = stop./30;
            STIM_BRFS.dichop_800soa_brfsPSflash.trstart(obs,:)  = double(times(find(codes == 9,1,'first')));
            STIM_BRFS.dichop_800soa_brfsPSflash.trend(obs,:)    = double(times(find(codes == 18,1,'first')));
            STIM_BRFS.dichop_800soa_brfsPSflash.obs(obs,:)  = obs; 
            for f = 1:length(stimfeatures)
                 STIM_BRFS.dichop_800soa_brfsPSflash.(stimfeatures{f})(obs,:) = grating.(stimfeatures{f})(stim);
            end
        end
  
   end
    

   
% 'dichop_800soa_brfsNPSflash'
obs = 0; clear e codes stim times
   for  e = 1:length(pEvC)
    codes  = pEvC{e};
    stim   = find(grating.trial == e);
    times  = pEvT{e};
    if strcmp('dCOS',grating.stim(e))           && ... 
            grating.soa(e)          ==  800     && ...
            grating.s1_contrast(e)  >= .5     && ...
            grating.s2_contrast(e)  >= .5     && ...
            grating.s2_tilt(e)      ==  PARAMS.NPS
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
            disp('error, please check -- dichop_800soa_brfsNPSflash -- loop')
        end
        stop    =  pEvT{e}(stimoff);
        
        obs = obs +1;
        STIM_BRFS.dichop_800soa_brfsNPSflash.start1(obs,:) = start1;
        STIM_BRFS.dichop_800soa_brfsNPSflash.start2(obs,:) = start2;
        STIM_BRFS.dichop_800soa_brfsNPSflash.stop(obs,:) = stop;
        STIM_BRFS.dichop_800soa_brfsNPSflash.start1Down(obs,:) = start1./30;
        STIM_BRFS.dichop_800soa_brfsNPSflash.start2Down(obs,:) = start2./30;
        STIM_BRFS.dichop_800soa_brfsNPSflash.stopDown(obs,:) = stop./30;
        STIM_BRFS.dichop_800soa_brfsNPSflash.trstart(obs,:)  = double(times(find(codes == 9,1,'first')));
        STIM_BRFS.dichop_800soa_brfsNPSflash.trend(obs,:)    = double(times(find(codes == 18,1,'first')));
        STIM_BRFS.dichop_800soa_brfsNPSflash.obs(obs,:)  = obs; 
        for f = 1:length(stimfeatures)
             STIM_BRFS.dichop_800soa_brfsNPSflash.(stimfeatures{f})(obs,:) = grating.(stimfeatures{f})(stim);
        end
        
    end

   end
   



end



