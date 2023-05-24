function checkrolalignment(rol,rolnr)
%plot the displacement between 10000 randomly chosen rol and rolnr rolonies

for i=2:length(rol)
    a=rol{i};b=rolnr{i};
    I=randperm(size(a,1),min(10000,size(a,1)));
    figure;hold on;
    for n=1:min(10000,size(a,1))
        plot([a(I(n),1) b(I(n),1)],[a(I(n),2) b(I(n),2)],'r','LineWidth',2);
    end
    title(['Cycle ',num2str(i)]);
end
