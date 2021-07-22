

k = 3;
Pfinter = sort(Pfinter);
Psinter = sort(Psinter);
finterX = Pf(:,Pfinter);
sinterX = Ps(:,Psinter);
fTosIdx = knnsearch(sinterX',finterX','K',k);
sTofIdx = knnsearch(finterX',sinterX','K',k);

figure(1);clf; hold on;
plot3(Pf(1,Pfinter),Pf(2,Pfinter),Pf(3,Pfinter),'.');
plot3(Ps(1,Psinter),Ps(2,Psinter),Ps(3,Psinter),'.');

seei = 1;
plot3(sinterX(1,fTosIdx(seei,:)),sinterX(2,fTosIdx(seei,:)),sinterX(3,fTosIdx(seei,:)),'d');
hold off;

fout = fopen('FSIInter_TRI_K12.txt', 'w');
fprintf(fout,'F_TO_S_NUM %d K %d\n', numel(Pfinter),k);
fprintf(fout,['%d\t',repmat('%d\t',1,k),'\n'], [Pfinter',Psinter(fTosIdx)]' - 1);
fprintf(fout,'S_TO_F_NUM %d K %d\n', numel(Psinter),k);
fprintf(fout,['%d\t',repmat('%d\t',1,k),'\n'], [Psinter',Pfinter(sTofIdx)]' - 1);



fclose(fout);