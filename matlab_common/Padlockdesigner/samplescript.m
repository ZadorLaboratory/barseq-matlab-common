%designing padlock probes for layer specific IT markers.
accnum={}; %copy accession numbers from excel into accnum.
save accnum;
RPI=cell(1,1); %copy RPI from excel (6nt)
save('RPI6nt.mat','RPI');

[allprimers,allpadlocks,allTm,allloc,allgenes]=batchmakepadlocks(accnum, RPI,1);
save temp1

validtargetsnm32=validatepadlockdesigns(allpadlocks,allgenes(:,3));
save temp2
finalist=filteredvalidtargets(validtargetsnm32,allpadlocks,allprimers,allgenes);
[~,I]=sort_nat(finalist(:,2));
finalist=finalist(I,:);
save final
save('finalprobes.mat','finalist');
