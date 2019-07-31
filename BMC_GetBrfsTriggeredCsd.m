%% BMC_GetBrfsTriggeredCsd
% edit from BMC_AverageCSDEffect_matchedSessions in VSS2019 repository

% Concerned with 800ms post stimulus onset for...
%   Monoc PS
%   Monoc NPS
%   diPostSoa PS
%   diPostSoa NPS
% This will be done on the 8 matched sessions analyzed for VSS2019 in
% preperation for SfN 2019. Several session will not have NPS stimulus, and
% that is ok. A first step should still be plotting all of the conditions I
% have for each day before averaging across days.

% End varible goal: SessAvgCSD (day x channels [row] x CSD values at all 800 ms [column] )

clear
close all

%% 1. Editable Variables
savetitle = 'SessionAverageCSD';
pre = 100;
post = 800;
sinkAllocate = 'BMC_DfS'; %'BMC_DfS','BMC_ChanNum','Old_DfS','Old_ChanNum'
baseDirectory = 'E:\LaCie';
stepDirectory = 'E:\LaCie\all BRFS';


%% 2. Load SessionParams AND establish EvalSink. (or anything else in SessionParams)
cd(baseDirectory)
load('SessionParams.mat')
switch sinkAllocate
    case 'BMC_DfS'
        SessionParams.EvalSink = SessionParams.BMC_DfS;
    case 'Old_DfS'
        SessionParams.EvalSink = SessionParams.Old_DfS;
end

%% 3. Loop through available sessions for selected condition
cd(stepDirectory)
allfolders = dir(stepDirectory);
allfolders = allfolders(3:size(allfolders,1)); % cuts out the first two Dir outputs. Unnecessary '.' and '..'

% 3.a. establish recursive loop
    % 3.a.i. get the SessionParams data down to the conditions for matchedDays 
    count = 0;
 
for x = 1:size(SessionParams.MatchedExists,1)
    if SessionParams.MatchedExists(x)
        count = count+1;
        SessionParamsforMatched.Date(count,1) = SessionParams.Date(x);
        SessionParamsforMatched.el(count,1)   = SessionParams.el(x);
        SessionParamsforMatched.sortdirection(count,1) = SessionParams.SortDirection(x);
        SessionParamsforMatched.EvalSink(count,1) = SessionParams.EvalSink(x);
        SessionParamsforMatched.PS(count,1) = SessionParams.PS(x);
        SessionParamsforMatched.NPS(count,1) = SessionParams.NPS(x);
    end
end   
         
    % 3.a.ii. Full Loop for each condition
timerange = [-pre:post];
AllCSDaligned(:,:,:) = nan(100,[size(timerange,2)],[size(SessionParamsforMatched.Date,1)]);

for a = 1:size(SessionParamsforMatched.Date,1) %big loop
    clearvars -except savetitle AllCSDaligned a allfolders pre post SessionParamsforMatched 
    disp(a);
% 3.b. enter folder
for b = 1:size(allfolders,1)
    folderFound = strfind(allfolders(b).name,string(SessionParamsforMatched.Date(a)));
    if ~isempty(folderFound)
       sessionDay = strcat(allfolders(b).folder,filesep,allfolders(b).name,filesep);
       filenameGrating = strcat(allfolders(b).name,'_brfs001.gBrfsGratings');
       filenameNEV = strcat(allfolders(b).name,'_brfs001.nev');
       filenameNs2 = strcat(allfolders(b).name,'_brfs001'); % no extension b/c added on getLFP ==> ??
    end
end
cd(sessionDay)

% 3.c. load grating, NEV, 
    % 3.c.i. load grating
        grating = readBRFS(filenameGrating);
    % 3.c.ii. load NEV
        NEV = openNEV(filenameNEV,'noread','nomat','nosave'); 
        EventCodes = NEV.Data.SerialDigitalIO.UnparsedData - 128;
        EventTimes = double(NEV.Data.SerialDigitalIO.TimeStamp); 
        evtFs = double(NEV.MetaTags.SampleRes);
        [pEvC,pEvT] = parsEventCodesML(EventCodes,EventTimes);

% 3.d. get LFP for full session day (Channels x timepoints) and filter
[xLFP,EventCodes,EventTimes] = getLFP(filenameNs2,'ns2','eD',SessionParamsforMatched.sortdirection(a),sessionDay);
fc = 50;
fs = 1000;
[butter_b,butter_a] = butter(4,fc/(fs/2)); 
filtLFP = filtfilt(butter_b,butter_a,xLFP);
LFP = filtLFP'; %(Channels x tiempoints)

% 3.e. stimLFP (ConditionTrial x chan x timepoints)
%   if dicopNOsoa start_nosoa, if dicopWsoa_PSflash or dicopWsoa_NPSflash, start2.
%   MAJOR CONDITIONAL STATEMENT HERE. CALL FROM SESSIONPARAMS FOR PS/NPS

% 3.e.i. get gratingOnsets, then 3.e.ii. gratingOnsets goes into stimLFP

    % 3.e.i. get gratingOnsets. A structure with fields for each condition.
    cctmPS      = 0;
    cctmNPS     = 0;
    cctsoaPS    = 0;
    cctsoaNPS   = 0;
    clear gratingOnsets


       for  e = 1:length(pEvC)
           % 'monocPS'
           if strcmp('Binocular',grating.stim(e))           && ... 
                    grating.soa(e)          ==  800       && ...
                    grating.s1_contrast(e)  >= .5     && ...
                    grating.s2_contrast(e)  >= .5     && ...
                    grating.s1_tilt(e)      ==  SessionParamsforMatched.PS(a)  && ...         
                    grating.s2_tilt(e)      ==  SessionParamsforMatched.PS(a)        
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
                cctmPS = cctmPS +1;
                gratingOnsets.monocPS(cctmPS,:) = [start1 start2 stop];
                column.monocPS = 1;

           % 'monocNPS'     
           elseif strcmp('Binocular',grating.stim(e))           && ... 
                    grating.soa(e)          ==  800       && ...
                    grating.s1_contrast(e)  >= .5     && ...
                    grating.s2_contrast(e)  >= .5     && ...
                    grating.s1_tilt(e)      ==  SessionParamsforMatched.NPS(a)  && ...         
                    grating.s2_tilt(e)      ==  SessionParamsforMatched.NPS(a)        
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
                cctmNPS = cctmNPS +1;
                gratingOnsets.monocNPS(cctmNPS,:) = [start1 start2 stop];
                column.monocNPS = 1;



            % 'diWsoaPS'
           elseif strcmp('dCOS',grating.stim(e))           && ... 
                    grating.soa(e)          ==  800     && ...
                    grating.s1_contrast(e)  >= .5     && ...
                    grating.s2_contrast(e)  >= .5     && ...
                    grating.s2_tilt(e)      ==  SessionParamsforMatched.PS(a)
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
                cctsoaPS = cctsoaPS +1;
                gratingOnsets.diWsoaPS(cctsoaPS,:) = [start1 start2 stop];
                column.diWsoaPS = 2;

            % 'diWsoaNPS'
            elseif strcmp('dCOS',grating.stim(e))           && ... 
                    grating.soa(e)          ==  800     && ...
                    grating.s1_contrast(e)  >= .5     && ...
                    grating.s2_contrast(e)  >= .5     && ...
                    grating.s2_tilt(e)      ==  SessionParamsforMatched.NPS(a) 
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
                cctsoaNPS = cctsoaNPS +1;
                gratingOnsets.diWsoaNPS(cctsoaNPS,:) = [start1 start2 stop];
                column.diWsoaNPS = 2;        
            end

       end
   
    % 3.e.ii. gratingOnsets goes into stimLFP
    for y_mPS = 1:size(gratingOnsets.monocPS,1)
        stimtm = round(gratingOnsets.monocPS(y_mPS,column.monocPS)/30);% divide by 30 to convert to 1kHz. Note, LFP already in 1kHZ
        refwin = stimtm-pre:stimtm+post;
        stimLFP.monocPS(y_mPS,:,:) = LFP(:,refwin); % LFP is (chan x timepoints) stimLFP is (ConditionTrial x chan x timepoints)     
    end
    for y_mNPS = 1:size(gratingOnsets.monocNPS,1)
        stimtm = round(gratingOnsets.monocNPS(y_mNPS,column.monocNPS)/30);% divide by 30 to convert to 1kHz. Note, LFP already in 1kHZ
        refwin = stimtm-pre:stimtm+post;
        stimLFP.monocNPS(y_mNPS,:,:) = LFP(:,refwin); % LFP is (chan x timepoints) stimLFP is (ConditionTrial x chan x timepoints)     
    end
    for y_soaPS = 1:size(gratingOnsets.diWsoaPS,1)
        stimtm = round(gratingOnsets.diWsoaPS(y_soaPS,column.diWsoaPS)/30);% divide by 30 to convert to 1kHz. Note, LFP already in 1kHZ
        refwin = stimtm-pre:stimtm+post;
        stimLFP.diWsoaPS(y_soaPS,:,:) = LFP(:,refwin); % LFP is (chan x timepoints) stimLFP is (ConditionTrial x chan x timepoints)     
    end
    for y_soaNPS = 1:size(gratingOnsets.diWsoaNPS,1)
        stimtm = round(gratingOnsets.diWsoaNPS(y_soaNPS,column.diWsoaNPS)/30);% divide by 30 to convert to 1kHz. Note, LFP already in 1kHZ
        refwin = stimtm-pre:stimtm+post;
        stimLFP.diWsoaNPS(y_soaNPS,:,:) = LFP(:,refwin); % LFP is (chan x timepoints) stimLFP is (ConditionTrial x chan x timepoints)     
    end

% 3.f. calcCSD (ConditionTrials x Chan x timepoints)
transLFP = permute(stimLFP,[1,3,2]); % flip for the correct calcCSD dimension
 for f = 1:size(transLFP,1)
     oneSessionLFP(:,:) = transLFP(f,:,:);
    CSD(f,:,:) = calcCSD(oneSessionLFP); %CSD is (chan[reduced by 1 on each end],timepionts)
 end
 
 % 3.g. Average CSD across trials (Chan x timepoints)
avgCSD = squeeze(nanmean(CSD,1));

% 3.h. baseline subtract (Chan x timepoints)
bl = pre-50:pre-1;
for h = 1:size(avgCSD,1)
    blofChan = nanmean(avgCSD(h,bl),2);
    blCSD(h,:) = (avgCSD(h,:)-blofChan); 
end

% 3.i. Align CSD outputs
sinkindex = SessionParamsforMatched.EvalSink(a);
alignedMatrixIndex = 50-sinkindex+1+1; %needs +1 so it is an inclusive count, and another +1 to account for the fact that the index was set up for padded arrays.
for n = 1:size(blCSD,1)
    AllCSDaligned(alignedMatrixIndex,:,a) = blCSD(n,:); %SAVEALLCSD is (session,aligned-chan,timepoints)
    alignedMatrixIndex = alignedMatrixIndex+1;
end


end %end of 'a' loop, #BigLoop.
clearvars -except savetitle AllCSDaligned a allfolders condition pre post SessionParamsforMatched 

% 5. Average CSD Effect (ACE)
ACE = squeeze(nanmean(AllCSDaligned,3));

cd('E:\LaCie\VSS 2019 figs\190429 figs post MC meeting\filteredMatVar')
save(strcat(savetitle,'.mat'))

load gong.mat;
soundsc(y);