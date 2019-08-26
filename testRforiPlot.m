% testrforiplots

figure
subplot(2,2,1)
f_ShadedLinePlotbyDepth(BLavg.ns6CSD,chans,TM,[],1)
plot([0 0], ylim,'k')
title({'pEvT codes ns6CSD ',filename}, 'Interpreter', 'none')
xlabel('time (ms)')
ylabel('Contacts indexed down from surface')


subplot(2,2,2)
imagesc(TM,chans,BLavg.ns6fCSD); colormap(flipud(jet));
climit = max(abs(get(gca,'CLim'))*.8);
if PARAMS.SortDirection == 'descending'
    ydir = 'reverse';
elseif PARAMS.SortDirection == 'ascending'
    ydir = 'normal';
end
set(gca,'CLim',[-climit climit],'Ydir',ydir,'Box','off','TickDir','out')
hold on;
plot([0 0], ylim,'k')
clrbar = colorbar;
title({'pEvT codes ns6fCSD',filename}, 'Interpreter', 'none')
xlabel('time (ms)')
ylabel('Contacts indexed down from surface')
clrbar.Label.String = 'nA/mm^3';

subplot(2,2,3)
f_ShadedLinePlotbyDepth(BLavg_photo.ns6CSD,chans,TM,[],1)
title({'photo diode ns6CSD ',filename}, 'Interpreter', 'none')
xlabel('time (ms)')
ylabel('Contacts indexed down from surface')

subplot(2,2,4)
imagesc(TM,chans,BLavg_photo.ns6fCSD); colormap(flipud(jet));
climit = max(abs(get(gca,'CLim'))*.8);
if PARAMS.SortDirection == 'descending'
    ydir = 'reverse';
elseif PARAMS.SortDirection == 'ascending'
    ydir = 'normal';
end
set(gca,'CLim',[-climit climit],'Ydir',ydir,'Box','off','TickDir','out')
hold on;
plot([0 0], ylim,'k')
clrbar = colorbar;
title({'photo diode ns6fCSD',filename}, 'Interpreter', 'none')
xlabel('time (ms)')
ylabel('Contacts indexed down from surface')
clrbar.Label.String = 'nA/mm^3';








% % % set(gcf,'Position',[1 40 331 662]);
