%function to plot the output of the previous CAVsimulation
function [m2,s2]=plotCAVsimulation(output1,output2,output3,output4)


s=std([output1(:),output2(:),output3(:),output4(:)],0,2);
m=mean([output1(:),output2(:),output3(:),output4(:)],2);

m2=[m(1:2),m(3:4),m(5:6),m(7:8)];
s2=[s(1:2),s(3:4),s(5:6),s(7:8)];

figure;errorbar_groups(m2,s2,'errorbar_width',0.5,'bar_width',0.5,'bar_colors',[1 0 0;0 0 1]);
xlabel('areas')
ylabel('fraction of axon')
legend('injection in OB', 'injection in AC')

end

