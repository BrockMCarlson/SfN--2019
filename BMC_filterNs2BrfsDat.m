%BMC_filterNs2BrfsDat.m
%GOAL: process and trigger brfs data
%   Version 1.2
%   Brock Carlson -- created 9/9/19
%   Taken from BMC_plotSessionDat_brfs_NS2ONLY.m 
%
%   DOES NOT WORK WITH NS6 FILES 
%   DOES NOT TRIGGER TO PHOTO DIODE
%       please see BMC_rforiPlot_STANDARD.m if you need either of these
%       functions.


clear

%% EDITABLE VARIABLES
filename = {'190319_B_rfori001'}';
% % filename = {'160102_E_rfori002','160102_E_brfs001','160427_E_rfori002','160427_E_brfs001','160510_E_rfori002','160510_E_brfs001'}';
% % % filename = {'151221_E_rfori001','151221_E_brfs001',...
% % %     '151222_E_rfori003','151222_E_brfs002',...
% % %     '151223_E_rfori002','151223_E_brfs001',...
% % %     '151231_E_rfori002','151231_E_brfs001',...
% % %     '160102_E_rfori002','160102_E_brfs001',...
% % %     '160104_E_rfori002','160104_E_brfs001',...
% % %     '160108_E_rfori004','160108_E_brfs002',...
% % %     '160111_E_rfori002','160111_E_brfs001',...
% % %     '160128_I_rfori003','160128_I_brfs001',...
% % %     '160130_I_rfori003','160130_I_brfs001',...
% % %     '160131_I_rfori003','160131_I_brfs001',...
% % %     '160204_I_rfori004','160204_I_brfs002',...
% % %     '160211_I_rfori004','160211_I_brfs002',...
% % %     '160212_I_rfori003','160212_I_brfs001',...
% % %     '160215_I_rfori002','160215_I_brfs003',...
% % %     '160425_E_rfori004','160425_E_brfs002',...
% % %     '160427_E_rfori002','160427_E_brfs001',...
% % %     '160510_E_rfori002','160510_E_brfs001',...
% % %     '160523_E_rfori002','160523_E_brfs001',...
% % %     '161005_E_rfori003','161005_E_brfs001',...
% % %     '161007_E_rfori002','161007_E_brfs002',...
% % %     '161011_E_rfori001','161011_E_brfs002'}';
sinkAllocate = 'BMC_DfS';
nameSaveType = 'LFPandCSDof';


% Computer-specific editable variables 
[ret, hostname] = system('hostname');
if strcmp(getenv('USER'),'maierav')
    % @Alex -- fill in necessary information for your system here.
     savefiledir = 'G:\LaCie\SfN 2019--figsAndMatVars\SfN 2019 figs\brfs conditions diagnostics';
    sessionParamDir = 'G:\LaCie';
elseif strcmp(getenv('USERNAME'),'Brock Carlson')
    % variables for end of script
    savefiledir = 'G:\LaCie\SfN 2019--figsAndMatVars\SfN 2019 figs\SfN_fig1\Fig1wAllBrfsSessions';
    sessionParamDir = 'G:/LaCie';
elseif ~ispc && contains(hostname,'Brocks-MacBook-Air')
    savefiledir = 'G:\LaCie\SfN 2019--figsAndMatVars\SfN 2019 figs\SfN_fig1 -- MacBook\Fig1wAllBrfsSessions -- MacBook';
    sessionParamDir = '/Volumes/SfN_2019/';
end


%% LOOP THROUGH ALL FILES

for a = 1:size(filename,1)
    clearvars -except a filename sinkAllocate pre post TM savefiledir nameSaveType sessionParamDir hostname
    
    disp(filename{a})
    

%% Computer-specific directories 
if strcmp(getenv('USER'),'maierav')
    % @Alex -- fill in necessary information for your system here.
        addpath(genpath('/Users/alex 1/Desktop/LAB/Brock'));
    dataDirectory = strcat('/Users/alex 1/Desktop/LAB/Brock',filesep,filename{a}(1:8));    
elseif strcmp(getenv('USERNAME'),'Brock Carlson')
    addpath(genpath('G:\LaCie\all BRFS'));
    dataDirectory = strcat('G:\LaCie\all BRFS\',filename{a}(1:8));    
elseif ~ispc && contains(hostname,'Brocks-MacBook-Air')
    dataDirectory = strcat('/Volumes/SfN_2019/all BRFS/',filename{a}(1:8));   
end




%% LOAD IN SESSION DATA
tic

cd(sessionParamDir)
load('SessionParams.mat')
    
switch sinkAllocate
    case 'BMC_DfS'
        SessionParams.EvalSink = SessionParams.BMC_DfS;
    case 'Old_DfS'
        SessionParams.EvalSink = SessionParams.Old_DfS;
end

% get info needed for day I'm analyzing
dayID = strfind(SessionParams.Date',str2double(filename{a}(1:6)));
PARAMS = SessionParams(dayID,:);


%% make names to load later
% Neural data filenames
cd(dataDirectory)
readNEVfile = strcat(dataDirectory,filesep,filename{a},'.nev');
readNS2file = strcat(dataDirectory,filesep,filename{a},'.ns2');



%% LOAD NEURAL DATA

%% Load pin-by-pin
% Load pin-by-pin and order
%'noread' for Header
NS2_Header      = openNSx(readNS2file,'noread');
% % % % NS6_Header      = openNSx(readNS6file,'noread');
% create a logical array indexing the position of the contacts for the
% electrode penetrated into V1 for the day
contactLogical = strcmp({NS2_Header.ElectrodesInfo.ConnectorBank},PARAMS.V1bank{1,1}(2));
contactLabels = {NS2_Header.ElectrodesInfo(contactLogical).Label};
contactNum = size(contactLabels,2);
if contactNum ~= PARAMS.el
    error('Houston we have a problem') %I'm not totally sure this set-up will always work... BMC 8/26/19
end

count = 0;
for i = 1:size(contactLogical,2)
    if contactLogical(i) == 1
        count = count+1;
        disp(strcat('i=',num2str(count)))

        %sort electrode contacts for Plexon DfS
        contactName = strcat(PARAMS.V1bank{1,1},sprintf('%02d',i));
        idx = contains(contactLabels,contactName);
        contactPosition = find(idx);
        pinNum = sprintf('c:%u',contactPosition);
        NS2         = openNSx(readNS2file,pinNum,'read','uV');
        %preallocate
            if count == 1
               sampleNumNS2 = length(NS2.Data); 
               ns2DAT_predivide = zeros(contactNum,sampleNumNS2);
               elLabelsOut = strings(contactNum,1);
            end
        ns2DAT_predivide(count,:) = NS2.Data;
        elLabelsOut(count,:) = contactName;
        clear pinNum NS2 NS6 contactName idx
    else
        continue
    end
end
ns2DAT = ns2DAT_predivide./4;%convert units to  uV
%flip if NN,  most superficial channel on top, regardless of number
if strcmp(string(PARAMS.SortDirection), 'descending')
    disp('flipped for NN array')
    ns2DAT = flipud(ns2DAT);
    elLabelsOut = flipud(elLabelsOut);
end

    
%% FILTER AND DOWNSAMPLE

%% filter  LFP and downsample the ns6LFP immediatly to align with STIM.onsets in 1kHz
    fc = 50;
    fs = 1000;
    [butter_b,butter_a] = butter(4,fc/(fs/2));
    ns2LFP = nan(size(ns2DAT));
    for j = 1:contactNum
        disp(strcat('j=',num2str(j)))
        clear ns2FiltVec ns6FiltVec
        ns2FiltVec = ns2DAT(j,:)';
        ns2LFP(j,:) = filtfilt(butter_b,butter_a,ns2FiltVec);
    end
    
%% CUT based on cortex alignment

ns2LFP_CUT = ns2LFP(PARAMS.rftop:PARAMS.rfbtm,:);

    
%% CALCULATE CSD
 % use mod_iCSD for correct units
ns2CSD_CUT   = padarray(mod_iCSD(ns2LFP_CUT),[1 0],NaN,'replicate');
ns2CSD = nan(size(ns2LFP));
ns2CSD(PARAMS.rftop:PARAMS.rfbtm,:) = ns2CSD_CUT; % place the CSD matrix within the correct cortical depth.

%% END OF CODE

cd(savefiledir)
savename = strcat(nameSaveType,filename{a},'.mat');
save(savename,'ns2LFP','ns2CSD','-v7.3','-nocompression');

%% SCRIPT END
toc

end

load gong.mat;
sound(y);





