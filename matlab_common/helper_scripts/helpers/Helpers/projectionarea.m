function [meanprojectionwidth,stdprojectionwidth]=projectionarea(bigmatrix,bigmatrixcounts,threshold)


cortexpeak=max(bigmatrixcounts(:,1:22),[],2)>30;

%find max
m_brain=max(bigmatrix(cortexpeak,:),[],2);


j=1;
for minimum=threshold
c=sum(bigmatrix(cortexpeak,1:22)<minimum*repmat(m_brain,1,22),2);
h_width(j,:)=cumsum(hist(22-c,0:22));
j=j+1;

end


plot(0:300:22*300,(h_width))
% xlim([-0.5,22.5]);
xlabel('projection width (um)')
ylabel('cumulative number of neurones')
xlim([0 22*300])
ylim([0 size(m_brain,1)])
meanprojectionwidth=mean((22-c))/22*100;
stdprojectionwidth=std((22-c))/22*100;


