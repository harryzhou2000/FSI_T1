model = createpde;
a = 0.2;
b = 1;
c1 = 0.3;
c2 = 0.7;
d = 1;
R1 = [3,4,0,a,a,0,0,0,b,b]';
R2 = [3,4,-d,0,0,-d,c1,c1,c2,c2]';

gm = [R1,R2]    ;
sf = 'R1+R2';

ns = [char('R1');char('R2')];
ns = ns';
[g,bt] = decsg(gm,sf,ns);
[g,bt] = csgdel(g,bt);

geometryFromEdges(model,g);
pdegplot(model,'EdgeLabels','on')
axis equal
% xlim([-0.1,1.1]);
% ylim([-0.1,1.1]); 


%%
hm = 0.05;
generateMesh(model,'Hmax',hm,'GeometricOrder','linear');
pdeplot(model,'NodeLabels','off');

lim = 1e-6;
nodesa =  findNodes(model.Mesh,'box',[0-lim,0+lim],[-2,2]);
nodesb =  findNodes(model.Mesh,'box',[a-lim,a+lim],[-2,2]);
nodesc =  findNodes(model.Mesh,'box',[-2,2],[0-lim,0+lim]);
nodesd =  findNodes(model.Mesh,'box',[-2,2],[b-lim,b+lim]);
%%  
e = [nodesa,nodesc,nodesb,nodesd];
t = model.Mesh.Elements;
p = model.Mesh.Nodes;


%%
xs = p(1,:);
ys = p(2,:);
ts = t(1:3,:);
file = fopen('PEmesh3.txt','w');
% type 0: internal 1: disp 2: force 3: dispx 4:dispy
fprintf(file, '%d\t%d\t%d\t4\n',size(p,2),size(ts,2),size(e,2));
fprintf(file,'pTypes\n');
fprintf(file,'0 0 0 0\n');
fprintf(file,'1 2 0 0\n');
fprintf(file,'2 1 0 0\n');
fprintf(file,'3 1 -3e-5 0 \n');


fprintf(file,'P\n');



for it =  1:size(xs,2)
    etype = 0;
    if( Within(ys(it),b,lim) || Within(ys(it),0,lim))
        etype = 1;
    end
    if( Within(xs(it),0,lim) && (ys(it) >= b* 0.8 || ys(it) <= b*0.2))
        etype = 2;
    end
    if( Within(xs(it),-d,lim) && (ys(it) >= b*0.3 && ys(it) <= b*0.7))
        etype = 3;
    end
    fprintf(file,'%.15e\t%.15e\t%d\n',xs(it),ys(it),etype);
end
fprintf(file,'T\n');
for it =  1:size(ts,2)
    fprintf(file,'%d\t%d\t%d\n',ts(1,it),ts(2,it),ts(3,it));
end
fprintf(file,'ED\n');
% for it = 1:size(e,2)
%     
%     fprintf(file,'%d\t%d\n',e(1,it),e(2,it));
% end
fclose('all');