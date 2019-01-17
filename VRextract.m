function [VRV,VRVb,allbouts,swims]=VRextract(VR,sampling)

[s1,s3]=size(VR);
sampling=8.3333;
si=(1/sampling)*1000;

trig=13672; 
s1=length(VR);
t=1:1:s1; 
t=t-trig;
t=t*si/1000; 


for z=1:size(VR,2)
    VRV(:,z)=movingstd(VR(:,z),10,'central');
    VRV(:,z)=VRV(:,z)-nanmean(VRV(:,z));
end

f1=figure
for z = 1:size(VRV,2) 
    VRV(:,z)=zscore(VRV(:,z));
    plot(VRV(:,z));
    hold on
end
title('select spike threshold')
thismax=max(VRV(:,1))*1.5;
axis([0 length(VRV) -10 thismax]);

%get spike threshold
l=imline(gca,[0 thismax; length(VRV) thismax]);
position=wait(l); %[x1 y1; x2 y2]
threshold=(position(1,2)+position(2,2))/2;
close(gcf)

VRVb=VRV;
VRVb(VRVb<threshold)=0;
VRVb(VRVb>threshold)=1;

%figure; plot(VRV(:,1)/10); hold on; plot(VRVb(:,1)+1.5);

omrgo=(0:10000:50000)+2500;
omrstop=omrgo+5000;
okrgo=[60000 70000 80000];
okrstop=[70000 80000 90000];
blkon=[0 2000 4000 6000 8000]+90000; 
blkoff=blkon+1000;

stimon=[omrgo okrgo blkon];
stimoff=[omrstop okrstop blkoff];

figure;
subplot(1,3,1)
k=0;
for z=1:size(VRV,2)
    plot(t,VRV(:,z)+k,'k');
    hold on
    k=k+20;
end
ymax=k;
axis([0 100000 -10 ymax]);
title('Trial behavior');
for j=1:14
    if j<4
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-60 -60 ymax ymax],'b');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);
    elseif j>3 & j<7
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-60 -60 ymax ymax],'g');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);    
    elseif j==7
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-60 -60 ymax ymax],'r');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);    
    elseif j==8
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-60 -60 ymax ymax],'m');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);  
    elseif j==9
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-60 -60 ymax ymax],'y');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);  
    else
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-60 -60 ymax ymax],'k');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);  
    end
end

subplot(1,3,2)
k=0;
for z=1:size(VRV,2)
    plot(t,VRVb(:,z)+k,'k');
    hold on
    k=k+1.5;
end
ymax=k;
axis([0 100000 -0.2 ymax]);
title('Trial binary');
for j=1:14
    if j<4
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-60 -60 ymax ymax],'b');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);
    elseif j>3 & j<7
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-60 -60 ymax ymax],'g');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);    
    elseif j==7
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-60 -60 ymax ymax],'r');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);    
    elseif j==8
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-60 -60 ymax ymax],'m');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);  
    elseif j==9
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-60 -60 ymax ymax],'y');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);  
    else
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-60 -60 ymax ymax],'k');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);  
    end
end

%saveas(gcf,'beh_summary.jpg');

 %VR extract bouts
 
 allbouts=[];
 alldurations=[];
 swims=struct;
 IBI=200; %define a "hard" inter-bout interval minimum in ms
 % increasing IBI means that bouts are more likely to be split up
 for z = 1:size(VRV,2)
    clear diffa theseons allons gaps alloffs bouts duration theseoffs
    temp=VRVb(:,z);
    activity=find(temp==1);
    diffa=diff(temp);
    allons=find(diffa==1);
    alloffs=find(diffa==-1);
    if alloffs(1)<allons(1)
        alloffs(1)=[]; %remove things that start "on"
    end
    if isempty(allons)==0;
        theseons(1)=allons(1); %find first on
        for p=1:length(allons)-1;
            gaps(p)=allons(p+1)-alloffs(p);
        end
        list=find(gaps>round(IBI*sampling))+1;
        theseons=[allons(1) allons(list)'];
        lastoff=alloffs(end);
        theseoffs=alloffs(list-1);
        theseoffs=[theseoffs' lastoff];

        bouts=zeros(size(temp));
        duration=[]; list=[];
        
        if length(theseons)~=length(theseoffs)
           % do something 
        end
        
        for j=1:length(theseons)
            duration(j)=theseoffs(j)-theseons(j); %in pts
            meanb=mean(VRVb(theseons(j):theseoffs(j),z)); %for a real bout, many of the points should be above threshold!
            if duration(j)<round(50*sampling); %50ms for min "bout"?
                %duration(j)=nan;
                %theseoffs(j)=nan; theseons(j)=nan;
                list=[list j];
                %continue
            %bouts(theseons(j):theseoffs(j))=1;       
            
            elseif meanb<0.05
                list=[list j];
            end
        end
        
        theseoffs(list)=[]; theseons(list)=[]; duration(list)=[];
        for j=1:length(theseons)
            bouts(theseons(j):theseoffs(j))=1;
        end
        swims(z).ons=theseons; swims(z).offs=theseoffs;
        alldurations=[alldurations duration];

        swims(z).duration=duration;
    
    else
        bouts=zeros(size(temp));
    end
    
%     figure; plot(t,bouts); axis([0 100000 -1.2 1.2]);
%     hold on; plot(t,VRVb(:,z)*-1); 
%     title(z);
%     hold off;
%     waitforbuttonpress
%     close(gcf);
    
    allbouts(z,:)=bouts;
    
 end


subplot(1,3,3)
k=0;
for z=1:size(VRV,2)
    plot(t,allbouts(z,:)+k,'k');
    hold on
    k=k+1.5;
end
ymax=k;
axis([0 100000 -0.2 ymax]);
title('Trial bouts');
for j=1:14
    if j<4
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-60 -60 ymax ymax],'b');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);
    elseif j>3 & j<7
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-60 -60 ymax ymax],'g');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);    
    elseif j==7
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-60 -60 ymax ymax],'r');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);    
    elseif j==8
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-60 -60 ymax ymax],'m');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);  
    elseif j==9
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-60 -60 ymax ymax],'y');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);  
    else
        h=patch([stimon(j) stimoff(j) stimoff(j) stimon(j)],[-60 -60 ymax ymax],'k');
        set(gca,'children',flipud(get(gca,'children')));
        set(h,'EdgeColor','none','facealpha',0.3);  
    end
end

set(gcf, 'Position', get(0,'Screensize')); % Maximize figure
saveas(gcf,'beh_summary.jpg');
 

alldurations=round(alldurations/sampling); %to convert pts to ms
save('beh.mat','VR','VRV','VRVb','allbouts','alldurations','swims');

figure; histogram(alldurations,100); title('bout duration (ms) histogram');
saveas(gcf,'bout_hist.jpg');


end