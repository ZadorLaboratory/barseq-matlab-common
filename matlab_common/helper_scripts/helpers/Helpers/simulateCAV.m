%simultate bulk injection of CAV-cre into OB or AC and counting of labeled
%axons in OB, CC, SC, AC. replicates experiment done in Schwarz et al. 2015
function [output]=simulateCAV(inputmatrix,inputnormalizedmatrix,infectionthreshold,plotting)

%map brain areas on MAPseq target areas
OB=23; %olfactory bulb
CC=7; % CC is at 1.4 to 0.7, first slice is at 3, so this is sections 6 to 7
SC=11; % SC is at 1.1 to -2, so sections 7 to 16
AC=19; % AC is at -2.5, so section 19
    
    



density=[];
j=1;
for i=[OB,AC]
    infected=inputmatrix(:,i)>=infectionthreshold;
    density(j,:)=sum(inputnormalizedmatrix(infected,[OB,AC,CC,SC]),1);
    j=j+1;
end
output=density./repmat(sum(density,2),1,4); %rows are injections, columns are target areas

if plotting==1
figure;
bar(output');    
xlabel('areas')
ylabel('fraction of total projection')  
end
end
    
