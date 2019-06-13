% BMC edit of analyze EVP from MLanalysisONLINE

BRdatafile = '151231_E_brfs001';
drname = 'E:\LaCie\all BRFS\151231_E';
extension = 'ns2'; % THIS CODE DOES NOT DOWNSAMPLE OR FILTER DATA
el = 'eD';
sortdirection = 'ascending'; %  descending (NN) or ascending (Uprobe)
pre = 100;
post = 500;

flag_subtractbasline = true;
flag_halfwaverectify = false;

clear LFP EventCodes EventTimes DAT TM CSD CSDf corticaldepth y

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[LFP, EventCodes, EventTimes]= getLFP(BRdatafile,extension,el,sortdirection,drname);
triggerpoints = EventTimes(EventCodes == 23 | EventCodes == 25 | EventCodes == 27 | EventCodes == 29| EventCodes == 31);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[DAT, TM] = trigData(LFP, triggerpoints , pre, post);
EVP = median(DAT,3);
% deal w/ bad channes
switch BRdatafile
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
figure;
switch sortdirection
    case 'ascending' 
        corticaldepth = [1:size(LFP,2)] ;
    case 'descending'
        corticaldepth = fliplr([1:size(LFP,2)]);
end
f_ShadedLinePlotbyDepth(EVP,corticaldepth,TM,[],1)
title(BRdatafile,'interpreter','none')

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
f_ShadedLinePlotbyDepth(CSD,corticaldepth,TM,[],1)
title(BRdatafile,'interpreter','none')

%%
CSDf = filterCSD(CSD);

figure
switch sortdirection
    case 'ascending'
        y = [1:size(CSDf,1)]./10;
        ydir = 'reverse';
    case 'descending'
        y = fliplr([1:size(CSDf,1)]./10);
        ydir = 'normal';
end
imagesc(TM,y,CSDf); colormap(flipud(jet));
climit = max(abs(get(gca,'CLim'))*.8);
set(gca,'CLim',[-climit climit],'Ydir',ydir,'Box','off','TickDir','out')
hold on;
plot([0 0], ylim,'k')
c = colorbar;
title(BRdatafile,'interpreter','none')

