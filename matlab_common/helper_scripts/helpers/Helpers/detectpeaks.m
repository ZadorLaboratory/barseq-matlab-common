function [refinedpeaknumber,refinedpeaks]=detectpeaks(input,cortex,bulb,plotting,normalize)

% detect peaks

peakdist=3; %distance between peaks
for i=1:size(input,1)
    minpeakheight=max(input(i,[cortex,bulb]))/2; %set minimum peak height criterium
    
    if max(input(i,cortex))<minpeakheight %find peaks in the bulb, if there is nothing good in cortex
        peaknumber(i)=1;
        peak(i).position=bulb;
        peak(i).height=input(i,bulb);
    else %find peaks in cortex, including the endpoints
        [peak(i).height,peak(i).position]=findpeaks([0,input(i,cortex),0],'MINPEAKHEIGHT',minpeakheight,'MINPEAKDISTANCE',peakdist);
%         [peak(i).height,peak(i).position,peak(i).width,peak(i).prominence]=findpeaks2015([0,input(i,cortex),0],'MINPEAKHEIGHT',minpeakheight,'MINPEAKDISTANCE',peakdist,'WidthReference','halfheight');

        peak(i).position=peak(i).position-1;
        if isempty(bulb)==0
        if input(i,bulb)>=minpeakheight; %check for peaks in the bulb
            peak(i).height=[peak(i).height,input(i,bulb)];
            peak(i).position=[peak(i).position,bulb];
        end
        end
    end
    peaknumber(i)=size(peak(i).position,2);
end

%refine peaks in cortex by checking for deep valleys between them
refinedpeaknumber=[];
j=1;
for i=1:size(input,1);
    refinedpeaks(i).position=[];
    if isempty(bulb)==0
    cortexpeaks=peak(i).position~=bulb; %number of peaks in cortex
    else
        cortexpeaks=logical(ones(size(peak(i).position)));
    end
        m=1;
        cortexpeakheight=peak(i).height(cortexpeaks);
        cortexpeakposition=peak(i).position(cortexpeaks);
        for k=1:sum(cortexpeaks)
            cutoff=cortexpeakheight(k)/2;
            %left interval
            if cortexpeakposition(k)==1
                leftminval=0;
            else
                
            leftcrossing_tmp=find(input(i,1:cortexpeakposition(k))>=cortexpeakheight(k),2,'last');
            if size(leftcrossing_tmp,2)==2;
            leftcrossing=leftcrossing_tmp(1);
            else
                leftcrossing=[];
            end
            if isempty(leftcrossing)==1
                leftminval=0;
            else
            [leftminval,leftminpos]=min(input(i,leftcrossing:cortexpeakposition(k)));
            end
            end
            
            %right interval
            if cortexpeakposition(k)==max(cortex)
                rightminval=0;
            else
            rightcrossing_tmp=cortexpeakposition(k)-1+find(input(i,cortexpeakposition(k):max(cortex))>=cortexpeakheight(k),2,'first');
            if size(rightcrossing_tmp,2)==2;
            rightcrossing=rightcrossing_tmp(2);
            else
                rightcrossing=[];
            end
            
            if isempty(rightcrossing)==1;
                rightminval=0;
            else
            [rightminval,rightminpos]=min(input(i,cortexpeakposition(k):rightcrossing));
            end
            end
            prominence=cortexpeakheight(k)-max(leftminval,rightminval);
            
            if prominence>cutoff
                refinedpeaks(i).position(m)=cortexpeakposition(k);
                m=m+1;
                
                
            end
        end
        if isempty(bulb)==0
        bulbpeaks=peak(i).position==bulb; %number of peaks in cortex
        refinedpeaks(i).position=[refinedpeaks(i).position,peak(i).position(bulbpeaks)];
        end

        refinedpeaknumber(i)=size(refinedpeaks(i).position,2);

        
end




    %plot histogram of peak numbers
    figure;
    h=hist(refinedpeaknumber,0:5);
    
    
    if isempty(normalize)==1
    bar(0:5,h,0.5)
    else
        bar(0:5,h/sum(h),0.5)
    end
    
    xlim([-0.5 5.5]);
    xlabel('number of peaks')
    ylabel('barcodes')
  if plotting==1
  
    %plot traces with peaks indicated as circles
    figure;
    for i=1:size(input,1)
        plot(input(i,[cortex,bulb]));
        hold on;
        plot(refinedpeaks(i).position,input(i,refinedpeaks(i).position),'o');
        title(['barcode ',int2str(i)]);
        hold off
        pause();
    end
end








