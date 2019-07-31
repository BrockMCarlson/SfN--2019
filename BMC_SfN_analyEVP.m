
addpath(genpath('G:\LaCie\all BRFS'));
drname        = {'G:\LaCie\all BRFS\190319_B\','G:\LaCie\all BRFS\190327_B\',...
    'G:\LaCie\all BRFS\190410_B\','G:\LaCie\all BRFS\190715_B\'};
BRdatafile    = {'190319_B_evp002','190327_B_evp003','190410_B_evp006','190715_B_evp003'};

extension     = 'ns2'; % THIS CODE DOES NOT DOWNSAMPLE OR FILTER DATA
el            = {'eD','eD','eD','eA'};
sortdirection = 'ascending'; %  descending (NN) or ascending (Uprobe)
pre           = 50;
post          = 250;
chans         = [1:32];
trls          = [1:200];


flag_subtractbasline = false;
flag_halfwaverectify = false;

for a = 1:4
cd(drname{a})
clear LFP EventCodes EventTimes DAT TM CSD CSDf corticaldepth y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[LFP, EventCodes, EventTimes]= getLFP(BRdatafile{a},extension,el{a},sortdirection,drname{a});
triggerpoints = EventTimes(EventCodes == 23 | EventCodes == 25 | EventCodes == 27 | EventCodes == 29| EventCodes == 31);
if isempty(chans)
   chans = [1:size(LFP,2)];
end
LFP = LFP(:,chans);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[DAT, TM] = trigData(LFP, triggerpoints , pre, post);
if isempty(trls)
   EVP = mean(DAT,3);
else
EVP = mean(DAT(:,:,trls),3);
end
% deal w/ bad channes
switch BRdatafile{a}
   case {'151208_E_rfori001' '151208_E_rfori002','151218_I_evp002'}
       EVP(:,17) = mean([EVP(:,18), EVP(:,16)],2);
   case '151205_E_dotmapping003'
       EVP(:,18) = mean([EVP(:,17), EVP(:,19)],2);
   case {'160115_E_evp001', '160115_E_rfori001', '160115_E_rfori002','160115_E_mcosinteroc001'}
       EVP(:,end-3:end) = [];
   case '160831_E_evp001'
       EVP(:,21) = mean([EVP(:,20), EVP(:,22)],2);
end

%%
%figure;
switch sortdirection
   case 'ascending'
       corticaldepth = [chans] ;
   case 'descending'
       corticaldepth = fliplr([chans]);
end
%f_ShadedLinePlotbyDepth(EVP,corticaldepth,TM,[],1)
%title(BRdatafile{a},'interpreter','none')

%%
CSD = calcCSD(EVP);
if flag_subtractbasline
   CSD = bsxfun(@minus,CSD,mean(CSD(:,TM<0),2));
end
if flag_halfwaverectify
   CSD(CSD > 0) = 0;
end
CSD = padarray(CSD,[1 0],NaN,'replicate');
figure
subplot(1,2,1)
f_ShadedLinePlotbyDepth(CSD,corticaldepth,TM,[],1)
title(BRdatafile{a},'interpreter','none')
set(gcf,'Position',[1 40  700   785]);
%%
CSDf = filterCSD(CSD);

subplot(1,2,2)
switch sortdirection
   case 'ascending'
       y = [chans];
       ydir = 'reverse';
   case 'descending'
       y = fliplr([chans]);
       ydir = 'normal';
end
imagesc(TM,y,CSDf); colormap(flipud(jet));
climit = max(abs(get(gca,'CLim'))*.8);
set(gca,'CLim',[-climit climit],'Ydir',ydir,'Box','off','TickDir','out')
hold on;
plot([0 0], ylim,'k')
c = colorbar;
caxis([-200 200])

%% Save and Append a pdf
cd('G:\LaCie\SfN 2019\SfN 2019 figs\L4 sink 2019 brfs sessions')
if a ==1
    export_fig('NonSubBaseline_L4NewSessions','-pdf','-nocrop') 
else
    export_fig('NonSubBaseline_L4NewSessions','-pdf','-nocrop','-append')
end

end