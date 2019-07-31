function [STIM] = getEventTimeInfo(brdrname,BRdatafile)

% LOAD TEXT FILE WITH STIMULUS CONDITIONS :
if any(strfind(BRdatafile,'drfori'))
    ext         = '.gRFORIDRFTGrating_di';
    grating     = readgDRFTGrating([brdrname BRdatafile ext]);
elseif any(strfind(BRdatafile,'rforidrft'))
    ext         = '.gRFORIDRFTGrating_di';
    grating     = readgDRFTGrating([brdrname BRdatafile ext]);
elseif any(strfind(BRdatafile,'rfsfdrft'))
    ext         = '.gRFSFDRFTGrating_di';
    grating     = readgDRFTGrating([brdrname BRdatafile ext]);
elseif any(strfind(BRdatafile,'posdisparitydrft'))
    ext         = '.gPOSDISPARITYDRFTGrating_di';
    grating     = readgDRFTGrating([brdrname BRdatafile ext]);
elseif any(strfind(BRdatafile,'disparitydrft'))
    ext         = '.gDISPARITYDRFTGrating_di';
    grating     = readgDRFTGrating([brdrname BRdatafile ext]);
elseif any(strfind(BRdatafile,'cinterocdrft'))
    ext         = '.gCINTEROCDRFTGrating_di';
    grating     = readgDRFTGrating([brdrname BRdatafile ext]);
elseif any(strfind(BRdatafile,'coneinterocdrft'))
    ext         = '.gCONEINTEROCDRFTGrating_di';
    grating     = readgDRFTGrating([brdrname BRdatafile ext]);
elseif any(strfind(BRdatafile,'conedrft'))
    ext         = '.gCONEDRFTGrating_di';
    grating     = readgDRFTGrating([brdrname BRdatafile ext]);
elseif any(strfind(BRdatafile,'colorflicker'))
    ext         = '.gCOLORFLICKERDRFTGrating_di';
    grating     = readgDRFTGrating([brdrname BRdatafile ext]);
elseif any(strfind(BRdatafile,'bwflicker'))
    ext         = '.gBWFLICKERDRFTGrating_di';
    grating     = readgDRFTGrating([brdrname BRdatafile ext]);
elseif any(strfind(BRdatafile,'rfori'))
    ext         = '.gRFORIGrating_di';
    grating     = readgGrating([brdrname BRdatafile ext]);
elseif any(strfind(BRdatafile,'rfsize'))
    ext         = '.gRFSIZEGrating_di';
    grating     = readgGrating([brdrname BRdatafile ext]);
elseif any(strfind(BRdatafile,'cinteroc'))
    ext         = '.gCINTEROCGrating_di';
    grating     = readgGrating([brdrname BRdatafile ext]);
elseif any(strfind(BRdatafile,'color'))
    ext         = '.gCOLORGrating_di';
    grating     = readgGrating([brdrname BRdatafile ext]);
elseif any(strfind(BRdatafile,'rfsf'))
    ext         = '.gRFSFGrating_di'; 
    grating     = readgGrating([brdrname BRdatafile ext]);
elseif any(strfind(BRdatafile,'mcosinteroc'))
    ext         = '.gMCOSINTEROCGrating_di';
    grating     = readgGrating([brdname BRdatafile ext]);
elseif any(strfind(BRdatafile,'dotmapping'))
    ext         = '.gDotsXY_di';
    grating     = readgDotsXY([brdname BRdatafile ext]);
elseif any(strfind(BRdatafile,'brfs'))
    ext         = '.gBrfsGratings';
    grating     = readBRFS([brdrname BRdatafile ext]);
elseif any(strfind(BRdatafile,'evp'))
    grating     = [];  
end

% FILENAME FOR NEV FILE :
filename        = [brdrname '\' BRdatafile];
NEV             = openNEV(strcat(filename,'.nev'),'noread','overwrite');
EventCodes      = NEV.Data.SerialDigitalIO.UnparsedData - 128;
EventSampels    = NEV.Data.SerialDigitalIO.TimeStamp;
[pEvC, pEvT]    = parsEventCodesML(EventCodes,EventSampels);


if ~any(strfind(BRdatafile,'evp'))
[STIM]          = sortStimandTimeData(grating,pEvC,pEvT,'stim');
TP              = [STIM.onsets STIM.offsets];
else
    offset_idx  = find(EventCodes == 24 | EventCodes == 26);
    onset_idx   = find(EventCodes == 23 | EventCodes == 25); 
    for tr = 1:length(offset_idx)
       match_on_id(tr) = closestVal(onset_idx,offset_idx(tr)); 
    end
    TP(:,1)     = EventSampels(match_on_id); 
    TP(:,2)     = EventSampels(offset_idx); 
    TP          = double(TP); 
    STIM.ypos   = repmat(0, [size(TP,1) 1]); 
    STIM.xpos   = repmat(0, [size(TP,1) 1]); 
end


if any(strfind(BRdatafile,'evp'))
    [newTP,~,STIM]       = photoReTrigger(TP,[brdrname BRdatafile],nanunique(STIM.ypos),'default',STIM);
else
    [newTP,~,STIM]       = photoReTrigger(TP,[brdrname BRdatafile],nanunique(STIM.ypos),[],STIM);
end

if any(isnan(newTP(:,1)))
    [try_newTP]          = photoReTrigger(TP,[brdrname BRdatafile],nanunique(STIM.ypos),'constant',STIM);
    newTP(isnan(newTP))  = try_newTP(isnan(newTP));
end

if datenum(BRdatafile(1:6),'yymmdd') < datenum('160607','yymmdd')
    STIM.trg_photo = 0;
else
    STIM.trg_photo = 1;
end


if STIM.trg_photo && (any(strfind(BRdatafile,'drft')) || any(strfind(BRdatafile,'flicker'))) 
        tf           = mode(grating.temporal_freq); 
        period       = ceil(1000/tf/2).*double(NEV.MetaTags.SampleRes)./1000;
        newTP        = newTP - period;
end

switch BRdatafile
    case {'151205_E_rfori002','151205_E_rfori003','151205_E_cosinteroc001'}
        bhvfile =  '/Volumes/Toshiba External/151205_E_rfori005.bhv';
    case '161005_E_rfori001'
        bhvfile =  '/Volumes/Toshiba External/161005_E_rfori002.bhv';
    case '161006_E_evp001'
        bhvfile =  '/Volumes/Toshiba External/161006_E_rfori001.bhv';
    case '170724_I_evp003'
        bhvfile =  '/Volumes/Toshiba External/170724_I_rfori002.bhv';
    case '170724_I_mcosinteroc001'
        bhvfile =  '/Volumes/Toshiba External/170724_I_rfori002.bhv';
    case '180808_I_evp001'
        bhvfile = '/Volumes/TOSHIBA EXT/neurophys_data_for_ana/180808_I/180808_I_posdisparitydrft004.bhv';
    case {'160609_I_bwflicker001','160609_I_cinterocdrft013','160609_I_colorflicker003'}
        bhvfile = '/Volumes/TOSHIBA EXT/neurophys_data_for_ana/160609_I/160609_I_cinterocdrft023.bhv';
    case '180810_I_evp001'
        bhvfile = '/Volumes/TOSHIBA EXT/neurophys_data_for_ana/180810_I/180810_I_conedrft003.bhv';
    case {'180830_I_evp001','180830_I_evp002'}
        bhvfile = '/Volumes/TOSHIBA EXT/neurophys_data_for_ana/180830_I/180830_I_color003.bhv';
        case '180906_I_evp001'
            bhvfile = '/Volumes/TOSHIBA EXT/neurophys_data_for_ana/180906_I/180906_I_color001.bhv';
    case '181212_B_evp001'
        bhvfile = '/Volumes/TOSHIBA EXT/neurophys_data_for_ana/181212_B/181212_B_cinterocdrft003.bhv';
    case '181217_B_evp001'
        bhvfile = '/Volumes/TOSHIBA EXT/neurophys_data_for_ana/181217_B/181217_B_cinterocdrft002.bhv';
    case '190118_B_evp001'
        bhvfile = '/Volumes/TOSHIBA EXT/neurophys_data_for_ana/190118_B/190118_B_cinterocdrft001.bhv';
    case '190119_B_evp001'
        bhvfile = '/Volumes/TOSHIBA EXT/neurophys_data_for_ana/190119_B/190119_B_cinterocdrft002.bhv';
    case '190120_B_evp001'
        bhvfile = '/Volumes/TOSHIBA EXT/neurophys_data_for_ana/190120_B/190120_B_cinterocdrft001.bhv';
    case '190123_B_evp001'
        bhvfile = '/Volumes/TOSHIBA EXT/neurophys_data_for_ana/190123_B/190123_B_cinterocdrft001.bhv';
    case '190209_B_evp002'
        bhvfile = '/Volumes/TOSHIBA EXT/neurophys_data_for_ana/190209_B/190209_B_cinterocdrft001.bhv';
    case '190210_B_evp001'
        bhvfile = '/Volumes/TOSHIBA EXT/neurophys_data_for_ana/190210_B/190210_B_cinterocdrft001.bhv';
    case '190319_B_evp002'
        bhvfile = '/Volumes/TOSHIBA EXT/neurophys_data_for_ana/190319_B/190319_B_cinterocdrft001.bhv';
    otherwise
        bhvfile = [filename '.bhv'];
end


% check for BHV and NS6
bhv                   = concatBHV(bhvfile); 
STIM.refresh          = bhv.VideoRefreshRate; 
STIM.measured_refresh = bhv.ActualVideoRefreshRate;
STIM.onsets           = [];
STIM.offsets          = [];
STIM.onsets           = newTP(:,1);
STIM.offsets          = newTP(:,2);
[STIM]                = removeTrialsfromStruct(STIM,find(isnan(STIM.onsets))); 








