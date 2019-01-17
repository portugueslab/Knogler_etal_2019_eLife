%% paths and directories

clear all
close all
clc

addpath('J:\Laura Knogler\matlab\');% make sure 'abfload.m' in folder
addpath('C:\Users\lknogler\Documents\code\knogler_et_al_2');
%cd('J:\_Shared\experiments\E0022_multistim\v01_ephys\'); %choose your directory
cd('C:\Users\lknogler\Documents\code\Knogler-et-al.-example-datasets\electrophysiology\ephys_example_PC1');

%% ephys traces and parameter loading

% select h5 file
info = h5info('PC_type1_DS_motion_onset.h5');
pc=h5read('PC_type1_DS_motion_onset.h5',strcat('/',info.Datasets.Name))';
info = h5info('VR_type1_DS_motion_onset.h5');
vr=h5read('VR_type1_DS_motion_onset.h5',strcat('/',info.Datasets.Name))';

sampling=8.3333; %in points/ms (KHz)
[s1,s3]=size(pc);
% s1 = length of recording in pts; data sampled at 8.333 KHz
% s3 = number of trials

trig=13672; % trigger between ephys trace start and stimulus start occurs with small fixed delay 

% set up time vector in ms with respect to stimulus onset
t=1:1:s1; 
t=t-trig;
t=t*si/1000; 

% define the ranges of interest for Knogler et al. 2019 visual stimuli
omrgo=(0:10000:50000)+2500;
omrstop=omrgo+5000;
okrgo=[60000 70000 80000];
okrstop=[70000 80000 90000];
blkon=[0 2000 4000 6000 8000]+90000+1000; 
blkoff=blkon+1000;
stimon=[omrgo okrgo blkon];
stimoff=[omrstop okrstop blkoff];

%% extract behavior from VR and save in 'beh.mat' structure

%this uses my VRextract.m function
[VRV,VRVb,allbouts,swims]=VRextract(vr,sampling);
s3=length(swims);

% if any trials look bad (poor S:N or bad noise artifacts), remove these trials from d [all channels] and rerun behavior

%% make sure everything is zeroed

% zero the baseline of all traces
p=floor(s3/8); if p==0; p=1; end;
for z = 1:s3
    tempmean=mean(pc(13000:200000,z)); 
    pc(:,z)=pc(:,z)-tempmean;
end

% view some overlaid unscaled trials to get sense of amplitudes, activity (Channel 1 only)
figure
for z = 1:p:s3
    p2=plot(pc(:,z)); 
    p2.Color(4) = 0.3; 
    axis([0 inf min(pc(50000:end,1))*2 max(pc(50000:end,1))*2]);
    hold on
    waitforbuttonpress
end
hold off
close(gcf)

% scale the traces to have the same max amplitude of CS, 
% to account for changes in seal noise/resistance during the recording
% note the baseline must be zeroed for this to work!
max_all=max(max(pc(trig+700:end-5,:)));

pc_mult=pc;
for z = 1:s3 %scaling factor
    max_temp=max(max(pc(100000:end-5,z)));
    scaling=max_all/max_temp;
    pc_mult(:,z)=pc_mult(:,z)*scaling;
end

for z = 1:s3 %rezeroing
    tempmean=mean(pc_mult(13000:200000,z)); 
    pc_mult(:,z)=pc_mult(:,z)-tempmean;
end

%view scaled overlaid trials - Channel 1
figure
for z = 1:s3
    p2=plot(pc_mult(:,z)); 
    p2.Color(4) = 0.2; 
    axis([0 inf min(pc_mult(50000:end,1))*2 max(pc_mult(50000:end,1))*2]);
    title(z)
    hold on
    waitforbuttonpress
end
hold off
close(gcf)

%% set the threshold to calculate spikes

f1=figure
for z = 1:p:s3
    plot(t,pc_mult(:,z));
    hold on
end
%axis([4000 12000 -1 30]); %get cropped view for good threshold setting
axis([35000 80000 min(pc_mult(50000:end,1)) max(pc_mult(50000:end,1))])
title('select SS threshold by dragging line and double clicking')

%get SS threshold
l=imline(gca,[5000 0.5; 80000 0.5]);
position=wait(l); %[x1 y1; x2 y2]
threshold_SS=(position(1,2)+position(2,2))/2;
close(gcf)

f1=figure
for z = 1:s3 %every 10th trial starting from first pairing
    h1=plot(t,pc_mult(:,z));
    h1.Color(4) = 0.3; 
    hold on
end
axis([2000 80000 min(pc_mult(50000:end,1)) max(pc_mult(50000:end,1))*1.5]); %get cropped view for good threshold setting
title('select CS threshold by dragging line and double clicking')

%get CS threshold
l=imline(gca,[2000 0.5; 80000 0.5]);
position=wait(l); %[x1 y1; x2 y2]
threshold_CS=(position(1,2)+position(2,2))/2;
close(gcf)

figure; plot(pc_mult(:,1,end)); hold on; 
l1=line([0 s1],[threshold_SS threshold_SS]); l1.Color='r'; l1.LineWidth=1.5;
title('zoom in to measure post-CS period to blank');
blankperiod=input('number of points to blank following a complex spike to avoid artifact?'); % adjust this as needed

%% find spikes, rates and save in 'PCspikes.mat' structure along with raw traces

rate_filter=1000*ones([1,1667])/(200); % a given sampling rate defines the rate filter; 1 spike in 200ms (1667 points) is 5Hz
t_trig=round(trig/sampling);

PCspikes=struct;
for z = 1:s3
    
    %first for CS
    traceCS=pc_mult(:,z);
    traceCS(1:trig+2500)=0; %blank points to avoid trigger onset artifact
    traceCS(traceCS<threshold_CS)=0;
    traceCS(traceCS>=threshold_CS)=1;
    dtraceCS=[0 diff(traceCS)'];
    spike_starts=find(dtraceCS==1);
    spike_ends=find(dtraceCS==-1);
    if numel(spike_ends)==0
        spike_starts=[];
        spike_ids=[];
    else
        if spike_ends(1)<spike_starts(1)
            spike_ends(1)=[];
        end
        if spike_starts(end)>spike_ends(end)
            spike_starts(end)=[];
        end
        spike_ids=floor((spike_starts+spike_ends)/2);
    end
    spikes=zeros(size(traceCS));
    spikes(spike_ids)=1;
    PCspikes(z).CS=spike_ids; 

    rate=conv(spikes,rate_filter,'same');
    PCspikes(z).CSrate=rate; %collect CS spike rates in structure
    
    %second for SS
    traceSS=pc_mult(:,z);
    traceSS(1:trig+3350)=0; %blank points to avoid trigger onset artifact
    PCspikes(z).trace=pc_mult(:,z);
    traceSS(traceSS<threshold_SS)=0; 
    traceSS(traceSS>=threshold_SS)=1;
    dtraceSS=[0 diff(traceSS)'];
    spike_starts=find(dtraceSS==1);
    spike_ends=find(dtraceSS==-1);
    if numel(spike_ends)==0
        spike_starts=[];
        spike_ids=[];
    else
        if spike_ends(1)<spike_starts(1)
            spike_ends(1)=[];
        end
        if spike_starts(end)>spike_ends(end)
            spike_starts(end)=[];
        end
        spike_ids=floor((spike_starts+spike_ends)/2);
    end

    for k=1:length(PCspikes(z).CS) %% this will remove CS artifacts from SS counts
        CS=PCspikes(z).CS;
        blankon=CS(k)-4; blankoff=CS(k)+blankperiod; %250 pts = 30 ms
        artifact=find(spike_ids<blankoff & spike_ids>blankon);
        spike_ids(artifact)=[];
    end
        
    spikes=zeros(size(traceSS));
    spikes(spike_ids)=1;
    PCspikes(z).SS=spike_ids; % collect SS spike times in struct
    rate=conv(spikes,rate_filter,'full');
    PCspikes(z).SSrate=rate(1:s1); % collect SS spike rates in struct
    
    %convert SS and CS spike times to real times (trig=0)
    PCspikes(z).SSt=round(PCspikes(z).SS/sampling)-t_trig;
    PCspikes(z).CSt=round(PCspikes(z).CS/sampling)-t_trig; 
    
end

save('PCspikes.mat','PCspikes','blankperiod','threshold_SS','threshold_CS');


%% run some checks to make sure that the analysis looks good;

% if not, reset the CS and or SS threshold, and/or the blankperiod

%check
z=randi(s3);
figure; 
hold on; plot(PCspikes(z).SSrate/5,'k');
plot(zscore(PCspikes(z).trace) +20,'b');
plot(zscore(VRV(:,z)) -20,'c');
plot(PCspikes(z).CS,ones(size(PCspikes(z).CS))*40,'r.','MarkerSize',20);
plot(PCspikes(z).SS,ones(size(PCspikes(z).SS))*10,'g.','MarkerSize',20);
%line([0 875000],[threshold_SS threshold_SS]);
axis([20000 s1 -25 max(zscore(PCspikes(z).SSrate(50000:end)))*10]);
ylabel('SS rate / 10');
legend('SSrate','raw trace','VR','CS raster','SS raster');
hold off;

% SS rate plots
figure;
k=0;
ymax=s3+1;
SSrates=[]; CSrates=[];
for z=1:length(PCspikes) 
    SS=smooth(PCspikes(z).SSrate,200); 
    SSrates(z,:)=PCspikes(z).SSrate;
    plot(t,SS+k,'k');
    hold on
    k=k+40;
end
ylabel('trials');
title('simple spike rates across trials');
axis([0 100000 -5 k+5]);
for j=1:14
    if j<4
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-5 -5 k+5 k+5],'b');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);
    elseif j>3 & j<7
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-5 -5 k+5 k+5],'g');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);    
    elseif j==7
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-5 -5 k+5 k+5],'r');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);    
    elseif j==8
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-5 -5 k+5 k+5],'m');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);  
    elseif j==9
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-5 -5 k+5 k+5],'y');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);  
    else
        h=patch([stimoff(j-1) stimon(j) stimon(j) stimoff(j-1)],[-5 -5 k+5 k+5],'k');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);  
    end
end
hold off

%% summary plot

load 'C:\Users\lknogler\Documents\code\knogler_et_al_2\Exp022_params.mat'

% choose your directory if you are not currently in the right one
load 'PCspikes.mat'; load 'beh.mat';
s3=length(PCspikes);

figure;
subplot(4,6,1:2) % swimming bouts across the experiment
imagesc(allbouts(:,trig:855338)); colorbar; 
title('bouts across the trial');

subplot(4,6,3:4) % mean SSrates
hold off; SSrates=[];
for z=1:length(PCspikes) 
    SS=smooth(PCspikes(z).SSrate,200); 
    SSrates(z,:)=PCspikes(z).SSrate;
    CSrates(z,:)=PCspikes(z).CSrate;
end
if size(SSrates,1)>1; plot(t,nanmean(SSrates),'k'); else plot(t,(SSrates),'k'); end; hold on;
thismax=max(SSrates(1,:))+1;
for j=1:14
    if j<4
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-1 -1 thismax thismax],'b');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);
    elseif j>3 & j<7
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-1 -1 thismax thismax],'g');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);    
    elseif j==7
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-1 -1 thismax thismax],'r');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);    
    elseif j==8
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-1 -1 thismax thismax],'m');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);  
    elseif j==9
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-1 -1 thismax thismax],'y');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);  
    else
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-1 -1 thismax thismax],'k');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);  
    end
end   
if size(SSrates,1)>1 plot(t,nanmean(allbouts)*-10,'r'); %else plot(t,(allbouts)*-10,'r');
end; hold off;
axis([0 100000 -10 thismax]);
title(['avg SSrate ', num2str(mean(nanmean(SSrates)),2) ' Hz']);

subplot(4,6,5:6) % CS raster plot
k=0; ymax=length(PCspikes) +1;
for z=1:length(PCspikes) 
    CS=PCspikes(z).CSt; 
    plot(CS,(ones(length(CS))+k),'k.');
    hold on
    k=k+1;
end
title('CS raster plot'); 
axis([0 t(end) 0 k+1]);
for j=1:14
    if j<4
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-1 -1 k+1 k+1],'b');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);
    elseif j>3 & j<7
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-1 -1 k+1 k+1],'g');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);    
    elseif j==7
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-1 -1 k+1 k+1],'r');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);    
    elseif j==8
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-1 -1 k+1 k+1],'m');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);  
    elseif j==9
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-1 -1 k+1 k+1],'y');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);  
    else
        h=patch([stimoff(j-1) stimon(j) stimon(j) stimoff(j-1)],[-1 -1 k+1 k+1],'k');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);  
    end
end  
hold off

subplot(4,6,11:12) %%  CS spike hist
histval=[]; 
for z=1:length(PCspikes) 
    CS=PCspikes(z).CSt;
    histval = [histval CS];
end
edges=[500:500:105000]; 
histogram(histval,edges);
yy=ylim(gca); yy=yy(2);
thismax=z*2.2; 
axis([0 t(end) 0 k*2]); 
title('CS histogram');
hold on
for j=1:14
    if j<4
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-1 -1 k*2 k*2],'b');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);
    elseif j>3 & j<7
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-1 -1 k*2 k*2],'g');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);    
    elseif j==7
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-1 -1 k*2 k*2],'r');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);    
    elseif j==8
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-1 -1 k*2 k*2],'m');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);  
    elseif j==9
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-1 -1 k*2 k*2],'y');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);  
    else
        h=patch([stimoff(j-1) stimon(j) stimon(j) stimoff(j-1)],[-1 -1 k*2 k*2],'k');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);  
    end
end
hold off

subplot(4,6,9:10) %  PSTH raster
hold off; allplots=[];k=1; allSSt=[];
ptime=linspace(-1000,1000,16667);
for j=1:length(PCspikes)
    CS=PCspikes(j).CS;
    SSrate=PCspikes(j).SSrate;
    SS=PCspikes(j).SS;
    for p=1:length(CS) 
        center=CS(p);
        if CS(p)<8334 
            CS(p)=nan;
            continue
        end
        if CS(p)>s1-8333
            CS(p)=[];
            break
        end
        on=CS(p)-8333;
        off=CS(p)+8333;
        thisrate=SSrate(on:off);
        thisSS=SS(SS>on & SS<off);
        thisSS(thisSS>CS(p)-15 & thisSS<CS(p))=0;
        thisSS=thisSS-center;
        allplots(k,:)=thisrate;
        plot(thisSS,ones(length(thisSS))+k,'k.')        
        allSSt=[allSSt round(thisSS/sampling)];
        hold on
        k=k+1;
        hold on;
    end
end
allplots(:,8334:8334+blankperiod)=nan;
h2=patch([-4 blankperiod blankperiod -4],[0 0 k+1 k+1],'r');
title('PSTH raster'); axis([round(-500*sampling) round(1000*sampling) 0 inf]); 

subplot(4,6,16) % PSTH hist
hold off; edges=[-500:10:1000]; 
histogram(allSSt,edges);
title('PSTH histogram');

subplot(4,6,22) % all bout durations by CS
hold off; CSs=struct;
timewindow=500; k1=0; k2=0;
latency1=[]; pct1=[]; num1=[]; n1=[]; dur1=[]; 
durneg1=[]; allCS1=[]; 
for z=1:length(swims) 
    ons=swims(z).ons; %time of ons in pts
    offs=swims(z).offs; 
    for o=1:length(ons)
        window1=ons(o); 
        window2=offs(o); 
        if window1>0 & window2<s1
            list=find(PCspikes(z).CS>window1 & PCspikes(z).CS<window2); 
            if length(list)>0
                CSt=PCspikes(z).CSt(list);
                offset=round((ons(o)-trig)/sampling); 
                allCS1=[allCS1 (CSt-offset)]; 
                dur=round((offs(o)-ons(o))/sampling); num1=[num1 length(list)]; n1=[n1 length(list)/dur];
                dur1=[dur1 dur]; k1=k1+1;
                latency1=[latency1 CSt(1)-offset]; pct1=[pct1 (CSt(1)-offset)/dur];
            elseif length(list)==0 
                durneg1=[durneg1 round((offs(o)-ons(o))/sampling)]; 
            end
        end 
    end
end
CSs.pos.dur=dur1; CSs.pos.latency=latency1; CSs.neg.dur=durneg1; save('CSs.mat','CSs');
percentCSpos=(length(dur1)/((length(dur1)+length(durneg1))))*100;
ecdf(dur1); hold on; ecdf(durneg1); legend('CS+','CS-'); axis([0 1000 0 1]); legend('Location','southeast');
h=get(gca,'children'); set(h(1),'LineWidth',1.5,'Color','g'); set(h(2),'LineWidth',1.5,'Color','b'); 
title(['durations (' num2str(percentCSpos,2) '% CS+)']);

subplot(4,6,7) % fOMR
hold off; timewindow=2500; 
start=omrgo(1:3); stop=omrstop(1:3);
start=start-timewindow; stop=stop+timewindow;
timebase=linspace(-timewindow,timewindow+5000,length(start(1):stop(1)));
omrrates=[]; 
k=1;
for z=1:length(PCspikes)
    for stim=1:3
        window=start(stim):stop(stim);
        window=round(window*sampling)+trig;
        omrrates(k,:)=PCspikes(z).SSrate(window);
        k=k+1;
    end
end
plot(timebase,nanmean(omrrates),'k');
thismax=max(nanmean(omrrates))*1.2;
h=patch([0 5000 5000 0],[0 0 thismax thismax],'b');
set(gca,'children',flipud(get(gca,'children')));
set(h,'EdgeColor','none','facealpha',0.3);
axis([-timewindow 5000+timewindow 0 thismax]);
title('forward OMR-triggered spike rate');

subplot(4,6,8) % flashes
timewindow=500; %window in ms
start=blkon(1:5); stop=blkoff(1:5); 
start=start-timewindow+1000;
stop=stop+timewindow+1000;
timebase=linspace(-timewindow,timewindow+1000,round((1000+(2*timewindow))*sampling));

whiterates=[]; k=1
for z=1:length(PCspikes)
    for stim=1:5
        window=start(stim):stop(stim);
        window=round(window*sampling)+trig;
        whiterates(k,:)=PCspikes(z).SSrate(window);
        k=k+1;
    end
end
timebase=linspace(-timewindow,timewindow+1000,length(whiterates));
plot(timebase,nanmean(whiterates),'k');
thismax=max(nanmean(whiterates))*1.2;
h=patch([-timewindow 0 0 -timewindow],[0 0 thismax thismax],'k');
set(gca,'children',flipud(get(gca,'children')));
set(h,'EdgeColor','none','facealpha',0.3);
h2=patch([1000 1000+timewindow 1000+timewindow 1000],[0 0 thismax thismax],'k');
set(gca,'children',flipud(get(gca,'children')));
set(h2,'EdgeColor','none','facealpha',0.3);
axis([-timewindow 1000+timewindow 0 thismax]);
title('flash-triggered spike rate');

subplot(4,6,13:14) % side and reverse OMR
hold off; timewindow=2500; %window in ms
start=omrgo(4:6); stop=omrstop(4:6);
start=start-timewindow;
stop=stop+timewindow;
timebase=linspace(-timewindow,timewindow+5000,length(start(1):stop(1)));
omrratesside1=[];
omrratesrev=[]; 
omrratesside2=[];
for z=1:length(PCspikes)
    for stim=1:3
        window=start(stim):stop(stim);
        window=round(window*sampling)+trig;
        if stim==1
            omrratesside1(z,:)=PCspikes(z).SSrate(window);
        elseif stim==2
            omrratesrev(z,:)=PCspikes(z).SSrate(window);
        elseif stim==3
            omrratesside2(z,:)=PCspikes(z).SSrate(window);
        end
    end
end
h1=plot(timebase,nanmean(omrratesside1),'r'); hold on;
h2=plot(timebase,nanmean(omrratesrev),'k');
h3=plot(timebase,nanmean(omrratesside2),'g');
thismax=max(nanmean(omrratesside2))*2;
h=patch([0 5000 5000 0],[0 0 thismax thismax],'b');
set(gca,'children',flipud(get(gca,'children')));
set(h,'EdgeColor','none','facealpha',0.3);
axis([-timewindow 5000+timewindow 0 thismax]);
legend([h1 h2 h3],{'leftward','reverse','rightward'});
title('other OMR-triggered spike rates');

subplot(4,6,19:20) % OKR sensory
hold off; timewindow=0; 
start=okrgo; stop=okrstop;
start=start-timewindow;
stop=stop+timewindow;
timebase=linspace(-timewindow,timewindow+10000,length(start(1):stop(1)));
OKRfull=[]; OKRside1=[]; OKRside2=[];
for z=1:length(PCspikes)
    for stim=1:3
        window=start(stim):stop(stim);
        window=round(window*sampling)+trig;
        if stim==1
            OKRfull(z,:)=PCspikes(z).SSrate(window);
        elseif stim==2
            OKRside1(z,:)=PCspikes(z).SSrate(window);
        elseif stim==3
            OKRside2(z,:)=PCspikes(z).SSrate(window);
        end
    end
end
h1=plot(timebase,nanmean(OKRfull),'r'); hold on;
h2=plot(timebase,nanmean(OKRside1),'m');
h3=plot(timebase,nanmean(OKRside2),'y'); 
thismax=max(nanmean(OKRfull))*2;
h=patch([0 10000 10000 0],[0 0 thismax thismax],'k');
set(gca,'children',flipud(get(gca,'children')));
set(h,'EdgeColor','none','facealpha',0.3);
axis([-timewindow 10000+timewindow 0 thismax]);
legend([h1 h2 h3],{'full','side1','side2'});
title('OKR-triggered spike rates');

subplot(4,6,15) % latency as cumulative histogram
hold off; ecdf(latency1); hold on;  ecdf(pct1*max(latency1)); title('%time to CS');
legend('abs latency','%latency scaled'); legend('Location','southeast');
title('latency to CS'); axis([0 max(latency1) 0 1]);

subplot(4,6,21) % IBIs following CS+ bouts   
hold off; CSposIBIs=[]; CSnegIBIs=[];
for z=1:length(swims) 
    ons=swims(z).ons; 
    offs=swims(z).offs; 
    for o=1:length(ons)-1
        window1=ons(o); 
        window2=offs(o); 
        if window1>0 & window2<s1
            list=find(PCspikes(z).CS>window1 & PCspikes(z).CS<window2); 
            if length(list)>0
                CSposIBIs=[CSposIBIs (swims(z).ons(o+1)-swims(z).offs(o))];
            elseif length(list)==0 
                CSnegIBIs=[CSnegIBIs (swims(z).ons(o+1)-swims(z).offs(o))];
            end
        end 
    end
end
CSposIBIs=CSposIBIs/sampling; CSnegIBIs=CSnegIBIs/sampling;
ecdf(CSposIBIs); hold on; ecdf(CSnegIBIs); title('IBIs by CSbout');
h=get(gca,'children'); set(h(1),'LineWidth',1.5,'Color','g'); set(h(2),'LineWidth',1.5,'Color','b'); 
legend('CS+','CS-'); legend('Location','southeast');

subplot(4,6,17) % SS bout on-triggered rates
hold off; SSbouts=[]; allswims=[]; 
everydur1=[];everydurhalf=[];
k=1; count=0; p=[];
timewindow=2000; %in ms, the amount of time to look at rate before and after bout start
for z=1:length(swims) 
    ons=swims(z).ons; 
    for o=1:length(ons)
        window1=ons(o)-round(timewindow*sampling); 
        window2=ons(o)+round(timewindow*sampling); 
        if window1>0 & window2<length(PCspikes(z).SSrate)
            thisrate=PCspikes(z).SSrate(window1:window2);
            thisswim=allbouts(z,window1:window2);
            SSbouts(k,:)=thisrate;
            allswims(k,:)=thisswim;
            k=k+1; count=count+1;
            tempdur=round(swims(z).duration(o)/sampling);
            everydur1=[everydur1 tempdur];
        end 
    end
end
t_bout=linspace(-timewindow,timewindow,length(thisrate));
h1=shadedErrorBar(downsample(t_bout,100),downsample(nanmean(SSbouts),100),downsample(std(SSbouts)/sqrt(count),100),'k',1);hold on; 
line([0 0],[0 30]); title('Bout on-triggered rate'); ylabel('Firing rate (Hz)'); xlabel('Time from onset (ms)'); SSbouts2=SSbouts;
axis([-500 1000 0 (1.3*max(nanmean(SSbouts2)))]); hold off

subplot(4,6,18) % SS bout off-triggered rates
SSbouts=[]; allswims=[];
k=1; count=0; p=[];
for z=1:length(swims) 
    offs=swims(z).offs; 
    for o=1:length(offs)
        window1=offs(o)-round(timewindow*sampling); 
        window2=offs(o)+round(timewindow*sampling); 
        if window1>0 & window2<length(PCspikes(z).SSrate)
            thisrate=PCspikes(z).SSrate(window1:window2);
            thisswim=allbouts(z,window1:window2);
            SSbouts(k,:)=thisrate;
            allswims(k,:)=thisswim;
            k=k+1;
            count=count+1;
        end 
    end
end
h1=shadedErrorBar(downsample(t_bout,100),downsample(nanmean(SSbouts),100),downsample(std(SSbouts)/sqrt(count),100),'k',1); hold on;
line([0 0],[0 30]); axis([-500 500 0 (1.3*max(nanmean(SSbouts2)))]);
title('Bout off-triggered rate'); ylabel('Firing rate (Hz)'); xlabel('Time from offset (ms)');
hold off

subplot(4,6,23) % CS bout on-triggered rates
CSbouts=[]; allswims=[]; 
k=1; count=0; p=[];
timewindow=1000; %in ms, the amount of time to look at rate before and after bout start
for z=1:length(swims) 
    ons=swims(z).ons; 
    for o=1:length(ons)
        window1=ons(o)-round(timewindow*sampling); 
        window2=ons(o)+round(timewindow*sampling); 
        if window1>0 && window2<length(PCspikes(z).CSrate)
            thisrate=PCspikes(z).CSrate(window1:window2);
            thisswim=allbouts(z,window1:window2);
            CSbouts(k,:)=thisrate;
            allswims(k,:)=thisswim;
            k=k+1;
            count=count+1;
        end 
    end    
end
t_bout=linspace(-timewindow,timewindow,length(CSbouts));
h1=shadedErrorBar(downsample(t_bout,100),downsample(nanmean(CSbouts),100),downsample(std(CSbouts)/sqrt(count),100),'k',1);
thismax=1.3*max(nanmean(CSbouts)); if thismax<1; thismax=1; end;
line([0 0],[0 30]); title('bout ON-triggered CSrate')
axis([-500 1000 0 thismax]); hold off

subplot(4,6,24) % CS bout off-triggered rates
CSbouts=[]; allswims=[]; 
k=1; count=0; p=[];
for z=1:length(swims) 
    offs=swims(z).offs; 
    for o=1:length(offs)
        window1=offs(o)-round(timewindow*sampling); 
        window2=offs(o)+round(timewindow*sampling); 
        if window1>0 & window2<length(PCspikes(z).CSrate)
            thisrate=PCspikes(z).CSrate(window1:window2);
            thisswim=allbouts(z,window1:window2);
            CSbouts(k,:)=thisrate;
            allswims(k,:)=thisswim;
            k=k+1;
            count=count+1;
        end 
    end
end
h1=shadedErrorBar(downsample(t_bout,100),downsample(nanmean(CSbouts),100),downsample(std(CSbouts)/sqrt(count),100),'k',1);
line([0 0],[0 30]); axis([-500 500 0 thismax]);
title('bout OFF-triggered CSrate');
hold off

set(gcf, 'Position', get(0,'Screensize'));
saveas(gcf,'summary.jpg');