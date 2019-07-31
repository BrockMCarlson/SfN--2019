function [newTP,trigger,STIM] = photoReTrigger_Brock(TP,filename)

close all; 
photolabel = 'ainp1';
newTP      = nan(size(TP));


% correct for certain bhv file issues
[brpath,BRdatafile,~] = fileparts(filename);

bhvfile = [filename '.bhv']; 
if ~exist(bhvfile,'file')
bhvfile = dir([brpath '/' BRdatafile(1:8) '*.bhv']); 
bhvfile = bhvfile(1).name; 
end

% check for BHV and NS6
bhv = concatBHV([brpath filesep bhvfile]);

ns6file = [filename '.ns6'];

% extract visual stim info from BHV and ypos
Fs        = 30000; % always using NS6
refresh   = ceil((1/bhv.ActualVideoRefreshRate) * Fs); % time for a screen refresh in sampels


% look at photodiode signal
ns_header  = openNSx(ns6file,'noread');
ns_labels  = cellfun(@(x) x(1:5),{ns_header.ElectrodesInfo.Label},'UniformOutput',0);
photoidx   = find(strcmp(ns_labels,photolabel));
if isempty(photoidx)
    newTP   = [];
    trigger = sprintf('no channel "%s"',photolabel);
    return
end

channel    = sprintf('c:%u', photoidx);

%  get photodiode signal around event
clear NS
data       = [];
if ~ns_header.RawData.PausedFile
    NS                 = openNSx(ns6file,channel,'read','precision','double','sample');
    data               = NS.Data;
else
    segfiles           = dir([ns6file(1:end-4) '-s' '*']);
    for s = 1:length(segfiles)
        NS   = openNSx([brpath '/' segfiles(s).name],channel,'read','precision','double','sample');
        data = [data NS.Data];
    end
    warning('paused file! double check this code\n'); 
end

thresh   = 1.4; 


for tr = 1%:size(TP,2)
    clear trstart trend refwin dat detected tf fewer
    
    trstart      = floor(TP{tr}(1) - (refresh*10));
    trend        = floor(TP{tr}(end) + (refresh*10));
    refwin       = trstart:trend; 
    dat          = data(refwin);
    dat          = abs((dat - mean(dat)) / std(dat));
    
    figure, plot(dat)
    
    if ~isfield(STIM,'temporal_freq')
        detected  = refwin(find(dat > thresh,1,'first')); % onset of first cycle
        STIM.photo_on{tr} = detected;
        newTP(tr,:)  = [detected nan];
    elseif STIM.temporal_freq(tr) <= 1
        detected  = refwin(find(dat > thresh,1,'first')); % onset of first cycle
        STIM.photo_on{tr} = detected;
        newTP(tr,:)  = [detected nan];
    else
        
        high         = refwin(dat > thresh);
        fewer        = high(diff([refwin(1) high]) > refresh);
        if isempty(fewer)
           thresh = 2;  
           high      = refwin(dat > thresh);
           fewer     = high(diff([refwin(1) high]) > refresh);
        end
        
        if length(fewer) > STIM.temporal_freq(tr)
            newTP(tr,:)       = [fewer(1) nan]; 
            STIM.photo_on{tr} = fewer(1:STIM.temporal_freq(tr));
        elseif isempty(fewer)
            newTP(tr,:)       = [nan nan]; 
            STIM.photo_on{tr} = nan; 
            warning('cant trigger photodiode\n'); 
        else
            newTP(tr,:)       = [fewer(1) nan]; 
            STIM.photo_on{tr} = fewer;
        end
        
        
        
        if any(strcmp(BRdatafile(1:8),'181207_B'))
            if isempty(fewer)
                newTP(tr,:)  = [TP(tr,1) nan];
            else
                newTP(tr,:)  = [fewer(1) nan];
            end
        end
        
%         if  mod(tr,10) == 0 & tr < 500
%         figure, plot(refwin,dat); v = vline(fewer(1)); set(v,'color','k','linewidth',2); 
%         end

    end
end


