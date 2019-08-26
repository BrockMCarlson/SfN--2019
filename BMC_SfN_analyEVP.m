
addpath(genpath('G:\LaCie\all BRFS'));
drname        = {'G:\LaCie\all BRFS\160102_E\'};%%,'G:\LaCie\all BRFS\160427_E\',...
   %% 'G:\LaCie\all BRFS\160510_E\'};
BRdatafile    = {'160102_E_rfori002'};%%'160427_E_brfs001','160510_E_brfs001'};
exportfigtext = {'rfori_160102'};%%'brfs_60427','brfs_160510'};


extension     = 'ns2'; % THIS CODE DOES NOT DOWNSAMPLE OR FILTER DATA
el            = {'eD'};%,'eD','eD'};
sortdirection = 'descending'; %  descending (NN) or ascending (Uprobe)
pre           = 50;
post          = 250;
chans         = [1:32];
trls          = [1:100];


flag_subtractbasline = true;
flag_halfwaverectify = false;

for a = 1:size(BRdatafile,2)
cd(drname{a})
clear LFP EventCodes EventTimes DAT TM CSD CSDf corticaldepth y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[LFP, EventCodes, EventTimes]= getLFP(BRdatafile{a},extension,el{a},sortdirection);
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
predivideCSD = calcCSD(EVP);
CSD = predivideCSD./4;
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


% % % % %% Save and Append a pdf
% % % % cd('G:\LaCie\SfN 2019--figsAndMatVars\SfN 2019 figs\L4 sink 2016 brfs sessions')
% % % % % % % if a ==1
% % % %     export_fig(exportfigtext{a},'-jpg','-nocrop') 
% % % % % % else
% % % % % %     export_fig(exportfigtext,'-pdf','-nocrop','-append')
% % % % % % % end

end