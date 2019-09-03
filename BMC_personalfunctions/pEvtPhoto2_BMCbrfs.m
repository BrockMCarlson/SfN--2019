function [pEvT_photo,phototrigger] = pEvtPhoto2_BMCbrfs(filename,pEvC,pEvT,ypos,bhv,photo_lORd,flag_plottriggeredBNC,method)

% filename     : BRdatafile with path, ONLY required input
% pEvC, pEvT   : cell arrays with event codes and times for every trial
% bhv          : ML's bhv varaible
% ypos         : grating location in DVA relative to the center of each 1/2 screen (grating.ypos)
% photo_lORd   : what the photodiode signal is called OR the actual BNC data

% pEvT_photo   : same as pEvT, but times are adjusted with photodiode signal (stimulus onset only)
% phototrigger : logical indicating if photodiode signal was used or not

% will retun pEvT_photo = [] and phototrigger = 'message'
% if photodiode signal is unavaible or not able to be analyze


% parse filename and check that it exists on disk
[brdrname,BRdatafile,extension] = fileparts(filename);
if isempty(extension) || ~strcmp(extension,'.ns6')
    extension = '.ns6';
    filename = fullfile(brdrname,[BRdatafile extension]);
end
bhvfile = fullfile(brdrname,[BRdatafile '.bhv']);
if ~exist(bhvfile,'file')
    pEvT_photo   = [];
    phototrigger = ('.bhv cannot be found');
elseif  ~exist(filename,'file')
    pEvT_photo   = [];
    phototrigger = sprintf('%s cannot be found',extension);
    return
end

% default input values
if ~exist('photo_lORd','var') || isempty(photo_lORd)
    photo_lORd = 'ainp1';
end

if ~exist('flag_plottriggeredBNC','var')  || isempty(flag_plottriggeredBNC)
    flag_plottriggeredBNC = 1;
end

% read from disk if not passed:
% pEvC and pEvT
if ~exist('pEvC','var') || ~exist('pEvT','var') || isempty(pEvC) || isempty(pEvT)
    nevfile = fullfile(brdrname,[BRdatafile '.nev']);
    NEV = openNEV(nevfile,'noread','nosave');
    EventCodes = NEV.Data.SerialDigitalIO.UnparsedData - 128;
    EventSampels = NEV.Data.SerialDigitalIO.TimeStamp;
    [pEvC, pEvT] = parsEventCodesML(EventCodes,EventSampels);
end
% ypos
if ~exist('ypos','var')  || isempty(ypos)
    if ~isempty(strfind(BRdatafile,'rom'))
        ext = '.gROMGrating_di';
    elseif ~isempty(strfind(BRdatafile,'rsf'))
        ext = '.gRSFGrating_di';
    elseif  ~isempty(strfind(BRdatafile,'rsiz'))
        ext = '.gRSZGrating_di';
    end
    grating = readgGrating([brdrname filesep BRdatafile ext]);
    ypos = mode(grating.ypos);
end
% bhv
if ~exist('bhv','var') || isempty(bhv)
    if ~exist(bhvfile,'file') % && str2num(BRdatafile(1:6)) < ...
            % 151222 
        warning('bhv does not exist, using sub')
        bhv = concatBHV( '/Volumes/harma/raw/151205_E_rom005.bhv' ); 
        bhv.PhotoDiodePosition{1} = 'None';
    else
        bhv = concatBHV(bhvfile);
    end
end
if ~exist('method','var') || isempty(method)
    n = datenum(BRdatafile(1:6),'yymmdd');
    if n < datenum('170216','yymmdd')
        method = 'default';
    else
        method = 'custom';
    end
end


% extract BHV info
RefreshRate = bhv.ActualVideoRefreshRate;
GratingConv = 0.5 + (ypos*bhv.PixelsPerDegree/bhv.ScreenYresolution); % accounts for CRT scan path
try
    FixONCode   = bhv.CodeNumbersUsed(strcmp(bhv.CodeNamesUsed,'Fixation spot ON'));
catch
    FixONCode = 35;
end
TrlENCode   = 18; %ML deault;
switch method
    case 'default'
        switch bhv.PhotoDiodePosition{1}
            case {'None'}
                phototrigger = false;
            case {'Lower left','Lower right'}
                phototrigger = true;
            otherwise
                pEvT_photo   = [];
                phototrigger = sprintf('bhv.PhotoDiodePosition returns unexpected value: %s',bhv.PhotoDiodePosition);
                return
        end
    case 'custom'
        phototrigger = true;
end

if ~phototrigger
    BNC = [];
else
    if ischar(photo_lORd)
        % Read in NS Header
        NS_header = openNSx(filename,'noread');
        
        % find relevant channel to get BNC data
        nidx = find(cell2mat(cellfun(@(x) ~isempty(strfind(x,photo_lORd)),{NS_header.ElectrodesInfo.Label},'UniformOutput',0)));
        channel    = sprintf('c:%u', nidx);
        NS = openNSx(filename,channel,...
            'read','precision','double');
        BNC = NS.Data;
    else
        BNC = photo_lORd;
    end
    
    % transform data to make it more usable
    BNC = -1*(BNC - median(BNC));
    
end

Fs         = 30000; % always using NS6
refresh    =  ceil((1/RefreshRate) * Fs); % time for a screen refresh in sampels
pEvT_photo = cell(size(pEvT)); X = []; Y = []; Z = []; failct = 0;
for t = 1: length(pEvC)
    
    codes = pEvC{t};
    times = nan(size(codes));
    onsets = find(pEvC{t} == 23 | pEvC{t} == 25 | pEvC{t} == 27 | pEvC{t} == 29 | pEvC{t} == 31);
    
    if isempty(onsets)
        pEvT_photo{t} = times;
    else
        
        if str2num(BRdatafile(1:6)) < 151231
            phototrigger = false;
        end
        
        if ~phototrigger
            % no photodiode to use for triggering,
            % assume constent offset relative to event code
            times(onsets) = pEvT{t}(onsets) + refresh*GratingConv;
            pEvT_photo{t} = times;
        else
            % use photodiode signal to trigger data
            
            % get threhsold values for trial
            zerotm = pEvT{t}(find(codes==FixONCode,1,'first'));
            endtm  = pEvT{t}(find(codes==TrlENCode,1,'first'));
            trlbnc = BNC(zerotm:endtm);
            trlsig = [0 diff(diff(trlbnc)) 0];
            thresh = 3*std(trlsig) + mean(trlsig);
            
            for p = 1:length(onsets)
                
                evt = pEvT{t}(onsets(p));
                
                switch method
                    case 'default'
                        
                        % determin if we're looking for a positive or negative deflection
                        win_pre    = evt - 3*refresh : evt;
                        win_post   = evt + refresh : evt + 4*refresh;
                        energy(1)  = std(diff(BNC(win_pre)));
                        energy(2)  = std(diff(BNC(win_post)));
                        
                        if diff(energy) < 0
                            % from high amplitude -> low amplitude
                        else
                            error('needs development')
                        end
                        
                        % find last inversion (this should be the last refresh before stimulus onset)
                        pre = []; k = 0;
                        while isempty(pre) && k < 2
                            k = k + 1;
                            win    = evt - k*refresh : evt - (k-1)*refresh;
                            signal = [0 fliplr(diff(diff(BNC(win)))) 0];
                            lim    = find(signal < thresh,1,'first');
                            if ~isempty(lim)
                                signal  = fliplr(BNC(win));
                                [~,idx] = max(signal);
                                pre = idx + (k-1)*refresh ;
                            end
                        end
                        
                        % adjust "evt" to TP with timing from photodiode signal
                        if isempty(pre)
                            failct = failct + 1;
                            disp('failed to re-trigger %u events',failct)
                            continue
                        else
                            TP = evt  - pre + refresh*GratingConv;
                        end
                        times(onsets(p)) = TP;
                        y = evt  - pre;
                        
                    case 'custom'
                        
                        % find first inversion after EVT
                        post = []; k = -1;
                        while isempty(post) && k < 2
                            k = k + 1;
                            win    = evt + k*refresh : evt + (k+1)*refresh;
                            signal = [0 fliplr(diff(diff(BNC(win)))) 0];
                            lim    = find(signal > thresh,1,'first') - 1;
                            if ~isempty(lim)
                                signal  = fliplr(BNC(win));
                                [~,idx] = max(signal);
                                post = idx + k*refresh ;
                            end
                        end
                        
                        % adjust "evt" to TP with timing from photodiode signal
                        if isempty(post)
                            failct = failct + 1;
                            fprintf('\nfailed to re-trigger %u events\n',failct)
                            continue
                        else
                            TP = evt  + post + refresh*GratingConv;
                        end
                        times(onsets(p)) = TP;
                        y = evt  + post;
                        
                end
                
                % save trigger points for plotting
                X = [X evt];
                Y = [Y y];
                Z = [Z TP];
                
            end
            
            pEvT_photo{t} = times;
            
            % save some data for plotting
            if flag_plottriggeredBNC
                plotdat{t,1} = pEvT{t}(onsets) - zerotm;
                plotdat{t,2} = trlbnc;
                plotdat{t,3} = trlsig;
            end
        end
    end
end


% if flag_plottriggeredBNC && phototrigger
%     
%     empty = cellfun(@isempty,plotdat(:,1));
%     plotdat(empty,:) = [];
%     
%     figure('Units','Inches','Position',[0 0 11 8.5]);
%     subplot(2,2,1)
%     plot(BNC)
%     axis tight; box off;
%     if ischar(photo_lORd)
%         title(sprintf('%s\n%s',BRdatafile,upper(photo_lORd)),'interpreter','none')
%     else
%         title(BRdatafile,'interpreter','none')
%     end
%     ylabel('|Photodiode Signal|')
%     
%     subplot(2,2,3)
%     plot([NaN diff(diff(BNC)) NaN])
%     axis tight; box off;
%     xlabel('Time (sampels)')
%     ylabel('d^2(|Photodiode Signal|)')
%     
%     subplot(2,2,2)
%     idx    = randi(size(plotdat,1),10,1);
%     offset = range(trlbnc);
%     for i = 1:length(idx)
%         plot(plotdat{idx(i),2} + offset*i); hold on;
%         plot(plotdat{idx(i),1},  offset*i, 'r*'); hold on;
%     end
%     axis tight; box off;
%     title('Example Trials')
%     ylabel('|Photodiode Signal|')
%     set(gca,'YTick',[])
%     
%     subplot(2,2,4)
%     idx    = randi(size(plotdat,1),10,1);
%     offset = range(trlsig);
%     for i = 1:length(idx)
%         plot(plotdat{idx(i),3} + offset*i); hold on;
%         plot(plotdat{idx(i),1},  offset*i, 'r*'); hold on;
%     end
%     axis tight; box off;
%     xlabel('Time (sampels)')
%     ylabel('d^2(|Photodiode Signal|)')
%     set(gca,'YTick',[])
%     
%     
%     figure('Units','Inches','Position',[0 0 11 8.5]);
%     p = 0;
%     for trig = 1:3
%         for dat = 1:2
%             p = p + 1;
%             clear DAT TP datstr trigstr
%             if dat == 1
%                 DAT =  BNC';
%                 datstr = '|Photodiode Signal|';
%             else
%                 DAT = [NaN; diff(diff(BNC')); NaN];
%                 datstr = 'd^2(|Photodiode Signal|)';
%             end
%             if trig == 1
%                 TP = X;
%                 trigstr = 'Triggered to Event Code';
%             elseif trig == 2
%                 TP = Y;
%                 trigstr = 'Triggered to Photodiode Peak';
%             else
%                 TP = Z;
%                 trigstr = 'Triggered to Photodiode Peak + Conv*Refresh';
%             end
%             
%             subplot(3,2,p)
%             [lfpDAT, TM] = trigData(DAT, TP, refresh*3, refresh*3);
%             TM = TM / Fs * 1000;
%             plot(TM,squeeze(lfpDAT));
%             axis tight; box off; hold on
%             plot([0 0],ylim,'-k');
%             for k = 1:2
%                 plot([-k -k]*refresh/Fs*1000,ylim,':k'); plot([k k]*refresh/Fs*1000,ylim,':k');
%             end
%             title(sprintf('%s\n%s',BRdatafile,trigstr),'interpreter','none')
%             xlabel('Time (ms)')
%             ylabel(datstr)
%         end
%     end
% end


%             if t == 4
%                 disp('bad tr')
%             end
%             clf; hold on
%             win = evt - 2*refresh : evt + 2*refresh;
%             x1  = BNC(win);
%             x2  = [0 diff(diff(win)) 0 ];
%             plot(x1);
%             plot(x2);
%             plot(2*refresh+1,0,'r*')
%            if isempty(post)
%                 plot(2*refresh+1 - pre,0,'g*');
%             else
%                 plot(2*refresh+1  + post,0,'g*');
%            end
%            title(sprintf('t = %u, p = %u',t,p));
% %             pause


% post = [];% k = 4;
%             while isempty(post) && k > 1
%                 k = k - 1;
%                 win = evt + (k-1)*refresh : evt + k*refresh;
%                 signal     = diff(diff(BNC(win)));
%                 if any(signal > thresh)
%                     [~,idx] = max(signal);
%                     post = idx + (k-1)*refresh + 1;
%                 end
%             end
%
%             adjust "evt" to TP with timing from photodiode signal
%             if isempty(post) && isempty(pre)
%                 failct = failct + 1;
%                 disp('failed to re-trigger %u events',failct)
%                 continue
%             elseif isempty(post)
%                 TP = evt  - pre;
%             else
%                 TP = evt  + post;
%             end
