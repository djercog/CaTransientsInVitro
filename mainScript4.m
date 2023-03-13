clear all
cd('C:\Users\GRO\Documents\Groc Lab\MATLAB');  %<-folder with codes

%1.Load the movie info...
% [roiV,tm]=procImageData('filename','pathToFilename',samplingRate);
% Notice: roiV{1}= data from background
fileName='ResultsILBuff';%'Basal_Env_APV_4n1'
filePath='C:\Users\GRO\Documents\Groc Lab\Experiments\Env\Ca_imaging\Results_Ca_ 2017\170110';
sampRate=20;
[roiV,tm]=procImageData(fileName,filePath,sampRate);

%2. nrBlocks=Number of blocks of the loaded movie to analyze
nrBlocks=3; %<-this could be defined automatically, if each block always contains 3000 frames
blockSize=length(tm)/nrBlocks;
colAux=rand(2*length(roiV),3);


%%plot de las raw trazas (not by block!)
clear strAux;
figure('color','w');hold on;
for i=2:length(roiV)
    plot(tm,roiV{i}.mean,'color',colAux(i,:),'Displayname',num2str(i));
    strAux{i-1}=['ROI ',num2str(i)];
end
legend(strAux);
ylabel('raw intensity');

%%%%% ANALYSIS
% first, segregate blocks from ROIs
for i=1:length(roiV)
    for j=1:nrBlocks
        roiVB{i,j}=roiV{i}.mean((j-1)*blockSize+1:blockSize*j);
    end
end

%transform deltaFF(t):
%(intensity(t)-mean(median(intensity([t-tmPrev,t])))/mean(median(intensity([t-tmPrev,t]))
clear deltaFF;
tmPrev=5;       %<-time window to consider as reference before actual point 
nrSampPrev=tmPrev*sampRate;
tmDeltaFF=tm(nrSampPrev+1:blockSize);
for i=1:length(roiV)
    for j=1:nrBlocks
        deltaFFB{i,j}=nan*zeros(1,length(tmDeltaFF));
        for s=nrSampPrev+1:length(roiVB{i,j})
            auxV=mean(median(roiVB{i,j}(s-nrSampPrev:s-1)));
            deltaFFB{i,j}(s-nrSampPrev)=(roiVB{i,j}(s)-auxV)/auxV;
        end
    end
end

%% FIGURES
%ALL deltaFF
clear strAux;
figure('color','w');hold on;
for i=2:size(deltaFFB,1)
    xAux2=[];auxV=[];
    for j=1:nrBlocks
        xAux=(j-1)*(blockSize/sampRate)+tmDeltaFF;
        xAux2=[xAux2,xAux];
        auxV=[auxV,deltaFFB{i,j}];
    end
    plot(xAux2,auxV,'color',colAux(i,:),'Displayname',num2str(i));
    strAux{i-1}=['ROI ',num2str(i)];
end
legend(strAux);

% 
% 
% %figure: raw, deltaF/F, derivative..
% i=5;        %<-choose a ROI
% nrStd=5;    %nr of std from current block
% clear ax;
% figure('color','w');hold on;
% ax(1)=subplot(3,1,1);
% plot(tm,roiV{i}.mean);
% ylabel('raw intensity');
% axis tight;
% ax(2)=subplot(3,1,2);hold on;
% for j=1:nrBlocks
%     xAux=(j-1)*(blockSize/sampRate)+tmDeltaFF;
%     devAux=quantile(deltaFFB{i,j},0.99);
%     plot([xAux(1),xAux(end)],[devAux,devAux],'-m');
%     devAux=nrStd*std(deltaFFB{i,j});
%     plot([xAux(1),xAux(end)],[mean(deltaFFB{i,j})+devAux,mean(deltaFFB{i,j})+devAux],'-r');
%     plot(xAux,deltaFFB{i,j},'b');
%     plot([xAux(1),xAux(end)],[mean(deltaFFB{i,j}),mean(deltaFFB{i,j})],'g');
% end
% legend('99%pc',[num2str(nrStd),' std']);
% ylabel('\DeltaF/F');
% axis tight;
% ax(3)=subplot(3,1,3);hold on;
% for j=1:nrBlocks
%     xAux=(j-1)*(blockSize/sampRate)+tmDeltaFF;
%     xAux=xAux(1:end-1)+(xAux(2)-xAux(1))/2;
%     deriv=diff(deltaFFB{i,j})*sampRate;
%     plot(xAux,deriv,'b');
%     devAux=quantile(deriv,0.99);
%     plot([xAux(1),xAux(end)],[devAux,devAux],'-m');
%     devAux=nrStd*std(deriv);
%     plot([xAux(1),xAux(end)],[mean(deriv)+devAux,mean(deriv)+devAux],'-r');
%     plot([xAux(1),xAux(end)],[mean(deriv),mean(deriv)],'g');
% end
% axis tight;
% ylabel('(\DeltaF/F)''');
% linkaxes(ax(1:3),'x');


% idea: if you block all the activity with APV, you should see the fluctuations
% reflecting the noise in your signal (theoretically, being zero).
% based on the noise from last block
% threshold on deltaF/F
minTmBetweenEvents=1;   %min amount of time between 2 bursts (in seconds)
nrStdFromLastBlock=5;   %choose a value for nr of std from last block (threshold)
i=8;                    %choose a value for Roi

clear ax;
figure('color','w');hold on;
ax(1)=subplot(3,1,1);
plot(tm,roiV{i}.mean);
ylabel('raw intensity');
axis tight;
ax(2)=subplot(3,1,2);hold on;
j=nrBlocks;
devAux=nrStdFromLastBlock*std(deltaFFB{i,j});
clear tmCrossAll;
kernelSize=10;  %smooth of signal
for j=1:nrBlocks
    xAux=(j-1)*(blockSize/sampRate)+tmDeltaFF;
    plot(xAux,deltaFFB{i,j},'b');
    deltaFFBSm=conv(deltaFFB{i,j}',ones(kernelSize,1)/kernelSize,'same');
    plot(xAux,deltaFFBSm,'m');
    %ev detection: crossing devAux
    halfTmStep=(xAux(2)-xAux(1))/2;
    indAux=int16(deltaFFBSm>devAux);
    indAux=indAux(2:end)'-indAux(1:end-1)';
    indAux=find(indAux==1);
    tmCross=xAux(indAux);  %tm of threshold crossing (rough)
    tmCross(find(diff(tmCross)<minTmBetweenEvents)+1)=[];
    if ~isempty(tmCross)
        plot(tmCross,devAux,'ok');
    end
    tmCrossAll(j)=length(tmCross);
end
title(['# Events per block =',num2str(tmCrossAll)]);
axis tight;
xlims=xlim;
plot([xlims(1),xlims(end)],[devAux,devAux],'-r');
legend({[num2str(nrStdFromLastBlock),' std (last block)'],['kernelSize=',num2str(kernelSize)]});
ylabel('\DeltaF/F');
axis tight;
ax(3)=subplot(3,1,3);hold on;
j=nrBlocks;
deriv=diff(deltaFFB{i,j})*sampRate;
devAux=nrStdFromLastBlock*std(deriv);   %compute deviation from last block
for j=1:nrBlocks
    xAux=(j-1)*(blockSize/sampRate)+tmDeltaFF;
    xAux=xAux(1:end-1)+(xAux(2)-xAux(1))/2;
    deriv=diff(deltaFFB{i,j})*sampRate;
    plot(xAux,deriv,'b');
end
axis tight;
xlims=xlim;
ylabel({'derivative','(\DeltaF/F)'''});
linkaxes(ax(1:3),'x');


%all individual deltaFFB derivative thresholded...
minTmBetweenEvents=1;   %min amount of time between 2 bursts (in seconds)
nrStdFromLastBlock=4; %<-choose a value for nr of std from last block (threshold)
kernelSize=10;  %smooth of signal
nrColumnsPlot=5;      %<-plot: number of columns and number of rows (figures)
nrRowsPlot=5;         % notice that if # ROIS > nrColumnsPlot*nrRowsPlot will give an error, then you have to change this number
if size(deltaFFB,1)>nrRowsPlot*nrColumnsPlot
    error('ERROR : nr ROIS > nrColumnsPlot*nrColumnsPlot');
end
clear tmCross*;
figure('color','w');k=1;
for i=1:size(deltaFFB,1)
    subplot(nrRowsPlot,nrColumnsPlot,i);hold on;
    j=nrBlocks; 
    devAux=nrStdFromLastBlock*std(deltaFFB{i,j});   %compute deviation from last block
    for j=1:nrBlocks
        xAux=(j-1)*(blockSize/sampRate)+tmDeltaFF;
        plot(xAux,deltaFFB{i,j},'b');
        deltaFFBSm=conv(deltaFFB{i,j}',ones(kernelSize,1)/kernelSize,'same');
        plot(xAux,deltaFFBSm,'m');
        
        %ev detection: crossing devAux
        halfTmStep=(xAux(2)-xAux(1))/2;
        indAux=int16(deltaFFBSm>devAux);
        indAux=indAux(2:end)'-indAux(1:end-1)';
        indAux=find(indAux==1);
        tmCross=xAux(indAux);  %tm of threshold crossing (rough)
        tmCross(find(diff(tmCross)<minTmBetweenEvents)+1)=[];
        if ~isempty(tmCross)
            plot(tmCross,devAux,'ok');
        end
        tmCrossAll{i,j}=tmCross;            %time of threshold crossing events
        tmCrossAllQty(i,j)=length(tmCross); %amount of crossing events
    end
    axis tight;
    xlims=xlim;
    plot([xlims(1),xlims(end)],[devAux,devAux],'-r');
    if rem(i,3)==1 
        if i==1
            ylabel({strrep(fileName,'_','\_'),'\DeltaF/F'});
        else
            ylabel('\DeltaF/F');
        end
    end
    if i==1
        title(['Roi',num2str(i),' (BG)']);
    else
        title(['Roi',num2str(i),'(#ev=',num2str(tmCrossAllQty(i,:)),')']);
    end
end
    



%autocorrelation x block
clear acgDeltaF dataByBlock;
acgT=0.05;
acgMaxLag=20;
for b=1:size(deltaFFB,2)
    for i=1:size(deltaFFB,1)
        [c,lag]=xcov(deltaFFB{i,b},floor(acgMaxLag/acgT), 'unbiased');
        v=var(deltaFFB{i,b});
        acgDeltaF{b}(i,:)=c/v;
        acgTm=acgT*lag;
    end
end

figure('color','w','position',[130 246 1137 195]);
for b=1:size(deltaFFB,2)
    subplot(1,size(deltaFFB,2),b);hold on;
    for i=1:size(acgDeltaF{b},1)
        plot(acgTm,acgDeltaF{b}(i,:),'color',colAux(i,:));
    end
    plot(acgTm,mean(acgDeltaF{b},1),'color','k','linewidth',2);
    title(['ACG',]);
    ylabel('AutoCorrelation (\DeltaF/F)');
    xlabel('lag (s)');
    ylim([-0.5,1]);
end

