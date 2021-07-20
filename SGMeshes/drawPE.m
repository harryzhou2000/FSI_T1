PE = readPE('sol_UNI0S1.txt');
plot3(PE.x,PE.y,PE.ux,'.');
hold on;
plot3(PE.x,PE.y,PE.uy,'.');
hold off;

%%
quiver(PE.x,PE.y,PE.ux,PE.uy,0.05);
%%
clf
exg = 0.1;
rep = getReport('report_UNI0S9.txt');
tx = [rep.x1,rep.x2,rep.x3];
ty = [rep.y1,rep.y2,rep.y3];
dtx = [rep.u1,rep.u2,rep.u3];
dty = [rep.v1,rep.v2,rep.v3]; 
c = rep.vms;
vmsmax = max(rep.vms,[],'all');
% patch((tx+exg*dtx)' ,(ty+exg*dty)',min(c',vmsmax/2))
patch((tx+exg*dtx)' ,(ty+exg*dty)',min(c',vmsmax*0.8),'EdgeColor','k')
colorbar;
axis equal
xlim([-3.1,-1.9]);