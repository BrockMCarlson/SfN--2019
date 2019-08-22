%BMC_photoTrigger
clear
close all

brdrname = 'G:\LaCie\all BRFS\'; 
brdrnamedir   = dir([brdrname '*_*']); 
expToAnaly = '_brfs001';

for i = 1:length(brdrnamedir)
    dirSet(i,:) = strcat(brdrname,brdrnamedir(i).name);
    fileSet(i,:) = brdrnamedir(i).name;

end
   


for a = 1:length(dirSet)
    disp(a)
clearvars -except a dirSet fileSet expToAnaly

brdrname = dirSet(a,:);
BRdatafile = fileSet(a,:);

% FILENAME FOR NEV FILE :
filename        = [brdrname filesep BRdatafile expToAnaly];
cd(brdrname)
NEV             = openNEV(strcat(filename,'.nev'),'noread','overwrite');
EventCodes      = NEV.Data.SerialDigitalIO.UnparsedData - 128;
EventSampels    = NEV.Data.SerialDigitalIO.TimeStamp;
[pEvC, pEvT]    = parsEventCodesML(EventCodes,EventSampels);

TP = getTP(pEvC,pEvT);


photolabel = 'ainp1';
newTP      = nan(size(TP));


% correct for certain bhv file issues
[brpath,BRdatafile,~] = fileparts(filename);

bhvfile = [BRdatafile '.bhv']; 
if ~exist(bhvfile,'file')
warning(strcat(BRdatafile,'_does not have bhv file'))
continue
end

% check for BHV and NS6
bhv = concatBHV([brpath filesep bhvfile]);

ns6file = [BRdatafile '.ns6'];
cd(brdrname)

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


% % for tr = 1:size(TP,2) %PROBLEM IN THIS LOOP, newTP IS SUPPOSED TO BE A Nx2 ARRAY????
% %     clear trstart trend refwin dat detected tf fewer
% %     
% %     trstart      = floor(TP{tr}(1) - (refresh*10));
% %     trend        = floor(TP{tr}(end) + (refresh*10));
% %     refwin       = trstart:trend; 
% %     dat          = data(refwin);
% %     dat          = abs((dat - mean(dat)) / std(dat));
% %     
% %     figure, plot(dat)
% %     
% %         high         = refwin(dat > thresh);
% %         fewer        = high(diff([refwin(1) high]) > refresh);
% %         if isempty(fewer)
% %            thresh = 2;  
% %            high      = refwin(dat > thresh);
% %            fewer     = high(diff([refwin(1) high]) > refresh);
% %         end
% %         
% %         if isempty(fewer)
% %             newTP(tr,:)       = [nan nan]; 
% %             STIM.photo_on{tr} = nan; 
% %             warning('cant trigger photodiode\n'); 
% %         end
% %   
% % end

hundreth = round([200+size(data,2)/100]);
startHu = hundreth*4;
endHu  = hundreth*5;
cutdata = data(:,startHu:endHu);
figure
plot(cutdata)
% plot title
titletext = {strcat('RawPhotoDiode...',BRdatafile)};
title(gca,titletext,'interpreter','none'); 
% Append a pdf
cd('G:\LaCie\SfN 2019\SfN 2019 figs\PhotoDiode')
if a ==1
    export_fig('RawPhotoDiode','-pdf','-nocrop') 
else
    export_fig('RawPhotoDiode','-pdf','-nocrop','-append')
end



normdat  = abs((cutdata - mean(cutdata)) / std(cutdata));
figure
plot(normdat)
% plot title
titletext = {strcat('NormPhotoDiode...',BRdatafile)};
title(gca,titletext,'interpreter','none'); 
% Append a pdf
cd('G:\LaCie\SfN 2019\SfN 2019 figs\PhotoDiode')
if a ==1
    export_fig('NormPhotoDiode','-pdf','-nocrop') 
else
    export_fig('NormPhotoDiode','-pdf','-nocrop','-append')
end



end