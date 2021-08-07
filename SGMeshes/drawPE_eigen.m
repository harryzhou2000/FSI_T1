PE = readPE('resultPE3.txt');
plot3(PE.x,PE.y,PE.ux,'.');
hold on;
plot3(PE.x,PE.y,PE.uy,'.');
hold off;

%%
quiver(PE.x,PE.y,PE.ux,PE.uy);
%%
clf
exg = 1e-1;
rep = getReportEigen('resultPE2_Eigen0.txt');
tx = [rep.x1,rep.x2,rep.x3];
ty = [rep.y1,rep.y2,rep.y3];
dtx = [rep.u1,rep.u2,rep.u3];
dty = [rep.v1,rep.v2,rep.v3]; 


% patch((tx+exg*dtx)' ,(ty+exg*dty)',min(c',vmsmax/2))
patch((tx+exg*dtx)' ,(ty+exg*dty)',rep.x1','EdgeColor','none')
colorbar;
axis equal