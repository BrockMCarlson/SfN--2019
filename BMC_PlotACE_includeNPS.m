%% BMC_PlotACE_includeNPS
% BMC adjusted from BMC_PlotACE_SSDP


clear

pre = 100;
post = 800;
TM = -pre:1:post;
AnalyzeSink = 'upper'; % 'lower' or 'upper'
exportfigtext = 'PerceptPlot2_upper';

color.biNOsoa = [5,113,176]/255; % light purple
color.diWsoaNPS = [202,0,32]/255;% light green
color.biWsoa = [123,50,148]/255 ; %dark purple
color.diWsoa = [0,136,55]/255; %dark green

%load session params for the title of each image
cd('G:\LaCie')
load('SessionParams.mat')
    count = 0;
 
for x = 1:size(SessionParams.MatchedExists,1)
    if SessionParams.PSexist(x) && SessionParams.NPSexist(x)
        count = count+1;
        SessionParamsForCondition.Date(count,1) = SessionParams.Date(x);
        SessionParamsForCondition.el(count,1)   = SessionParams.el(x);
        SessionParamsForCondition.sortdirection(count,1) = SessionParams.SortDirection(x);
        SessionParamsForCondition.EvalSink(count,1) = SessionParams.EvalSink(x);
        SessionParamsForCondition.PS(count,1) = SessionParams.PS(x);
        SessionParamsForCondition.NPS(count,1) = SessionParams.NPS(x);
        SessionParamsForCondition.V1bank(count,1) = SessionParams.V1bank(x);
    end
end   

% load and create structures based on conditions
cd('G:\LaCie\SfN 2019--figsAndMatVars')
    biWsoaPS	= load('PerceptPlot_biPSWsoafiltered.mat');    
    biWsoaNPS	= load('PerceptPlot_biNPSWsoafiltered.mat');
    diWsoaPS	= load('PerceptPlot_dicopWsoa_NPSflashfiltered.mat');
    diWsoaNPS	= load('PerceptPlot_dicopWsoa_PSflashfiltered.mat');    



cutAllCSDaligned.biWsoaNPS	= biWsoaNPS.AllCSDaligned(38:55,:,:);
cutAllCSDaligned.diWsoaNPS      = diWsoaNPS.AllCSDaligned(38:55,:,:);  
cutAllCSDaligned.biWsoaPS     = biWsoaPS.AllCSDaligned(38:55,:,:); 
cutAllCSDaligned.diWsoaPS     = diWsoaPS.AllCSDaligned(38:55,:,:);

%% Compare Sink Lineplots,
% cutACE row 13 is sink bottom, cortical depth of 0
% cutACE row 12 is also the sink(??), cortical depth of 0.1
switch AnalyzeSink
    case 'lower'
        sinkAllCSDaligned.biWsoaNPS	= cutAllCSDaligned.biWsoaNPS(12:13,:,:);
        sinkAllCSDaligned.diWsoaNPS     = cutAllCSDaligned.diWsoaNPS(12:13,:,:);  
        sinkAllCSDaligned.biWsoaPS	= cutAllCSDaligned.biWsoaPS(12:13,:,:); 
        sinkAllCSDaligned.diWsoaPS     = cutAllCSDaligned.diWsoaPS(12:13,:,:); 
        
        
        
        for i = 1:size(sinkAllCSDaligned.biWsoaNPS,3)
           sinkAvgAllCSDaligned.biWsoaNPS(i,:) = mean(sinkAllCSDaligned.biWsoaNPS(1:2,:,i),1);
           sinkAvgAllCSDaligned.diWsoaNPS(i,:) = mean(sinkAllCSDaligned.diWsoaNPS(1:2,:,i),1);
           sinkAvgAllCSDaligned.biWsoaPS(i,:) = mean(sinkAllCSDaligned.biWsoaPS(1:2,:,i),1);
           sinkAvgAllCSDaligned.diWsoaPS(i,:) = mean(sinkAllCSDaligned.diWsoaPS(1:2,:,i),1);
        end        
        titletext1 = 'CSD of lower sink. Cortical depth of 0 & 0.1 averaged.';
        
    case 'upper'
        sinkAvgAllCSDaligned.biWsoaNPS	= squeeze(cutAllCSDaligned.biWsoaNPS(6,:,:))';
        sinkAvgAllCSDaligned.diWsoaNPS     = squeeze(cutAllCSDaligned.diWsoaNPS(6,:,:))';  
        sinkAvgAllCSDaligned.biWsoaPS	= squeeze(cutAllCSDaligned.biWsoaPS(6,:,:))'; 
        sinkAvgAllCSDaligned.diWsoaPS     = squeeze(cutAllCSDaligned.diWsoaPS(6,:,:))';  
        titletext1 = 'CSD of upper sink. Cortical depth of 0.7.';
end



%% Process

% Pull out timecourse
postStimSink.diWsoaNPS  = sinkAvgAllCSDaligned.diWsoaNPS(:,:);
postStimSink.diWsoaPS   = sinkAvgAllCSDaligned.diWsoaPS(:,:);
postStimSink.biWsoaPS   = sinkAvgAllCSDaligned.biWsoaPS(:,:);
postStimSink.biWsoaNPS  = sinkAvgAllCSDaligned.biWsoaNPS(:,:);

%bl average each session
    % 3.h. baseline subtract (Chan x timepoints)
    bl = pre-50:pre-1;
    for g = 1:size(postStimSink.diWsoaNPS,1)
        blofChan = nanmean(postStimSink.diWsoaNPS(g,bl),2);
        blSink.diWsoaNPS(g,:) = (postStimSink.diWsoaNPS(g,:)-blofChan); 
    end
    for h = 1:size(postStimSink.diWsoaPS,1)
        blofChan = nanmean(postStimSink.diWsoaPS(h,bl),2);
        blSink.diWsoaPS(h,:) = (postStimSink.diWsoaPS(h,:)-blofChan); 
    end
        for h = 1:size(postStimSink.biWsoaPS,1)
        blofChan = nanmean(postStimSink.biWsoaPS(h,bl),2);
        blSink.biWsoaPS(h,:) = (postStimSink.biWsoaPS(h,:)-blofChan); 
        end
        for h = 1:size(postStimSink.biWsoaNPS,1)
        blofChan = nanmean(postStimSink.biWsoaNPS(h,bl),2);
        blSink.biWsoaNPS(h,:) = (postStimSink.biWsoaNPS(h,:)-blofChan); 
        end
    
        
% % % % Mean across trials
% % % AVG.diWsoaNPS = nanmean(blSink.diWsoaNPS,1);
% % % AVG.diWsoaPS = nanmean(blSink.diWsoaPS,1);
% % % AVG.biWsoaPS = nanmean(blSink.biWsoaPS,1);
% % % AVG.biWsoaNPS = nanmean(blSink.biWsoaNPS,1);
% % % 
% % % % bootci across trials
% % % ci.diWsoaNPS = bootci(4000,@mean,postStimSink.diWsoaNPS);
% % % ci.diWsoaPS = bootci(4000,@mean,postStimSink.diWsoaPS);
% % % ci.biWsoaPS = bootci(4000,@mean,postStimSink.biWsoaPS);
% % % ci.biWsoaNPS = bootci(4000,@mean,postStimSink.biWsoaNPS);
% % % 
% % % 

%% plot postSOA timecourse

postStimTM = -100:1:800;

for n = 1:size(blSink.diWsoaNPS,1)
titletext2 = SessionParamsForCondition.Date(n);
% plot sink for day with NPS day 1 only
h = figure;
plot(postStimTM,blSink.diWsoaNPS(n,:),'color','k','DisplayName','diWsoaNPS'); hold on;
plot(postStimTM,blSink.diWsoaPS(n,:),'color','r','DisplayName','diWsoaPS'); hold on;
plot(postStimTM,blSink.biWsoaPS(n,:),'color','g','DisplayName','monocPS'); hold on;
plot(postStimTM,blSink.biWsoaNPS(n,:),'color','b','DisplayName','monocNPS'); hold on;



    title({titletext2,titletext1})
    vline(0)
    hline(0)
    legend('Location','best')
    box off
    ax = gca;
    ax.YRuler.Exponent = 0;

cd('G:\LaCie\SfN 2019--figsAndMatVars\SfN 2019 figs\PerceptPlot')  
if n == 1
    export_fig(exportfigtext,'-pdf','-nocrop') 
else
    export_fig(exportfigtext,'-pdf','-nocrop','-append')
end
end
