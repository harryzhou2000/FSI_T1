%block Fluid
model = createpde('structural','static-solid');
importGeometry(model,'BLK.STL');

%%
pdegplot(model,'FaceLabels','on')

%%
Finf = [1,2,4,5,6];
    
%%
m1 = model.generateMesh('GeometricOrder','linear','Hmax',0.05,'Hmin',0.01);


%%
pdeplot3D(model,'Mesh','on')
%%

pfinf = model.Mesh.findNodes('region','Face',(Finf));
pfwall = model.Mesh.findNodes('region','Face',setdiff(1:model.Geometry.NumFaces,Finf));
pfall = model.Mesh.findNodes('region','Face',(1:model.Geometry.NumFaces));
pfinter = setdiff(pfall,pfinf);

[mp,me,mt] = model.Mesh.meshToPet;
mp(:,mt(1:4,213))
np = size(mp,2); nt = size(mt,2);
% plot3(mp(1,pfinter),mp(2,pfinter),mp(3,pfinter),'.')
model.Mesh.findElements('attached',mt(1:4,213))
%%
file = fopen('BLK_MeshCAD.txt','W');

fprintf(file,'NUM_POINT %d\n',np);
fprintf(file,'%.15e\t%.15e\t%.15e\n',mp);% pointcoord
fprintf(file,'NUM_TET %d\n',nt);
fprintf(file,'%d\t%d\t%d\t%d\n',mt(1:4,:)-1);% pointindex*4
fprintf(file,'NUM_BOUND %d\n',numel(pfall));
% pointindex-btype-bsetkey
fprintf(file,'%d\t%d\t%d\n',[pfinf-1,pfinter-1; ones(size(pfinf)),1*ones(size(pfinter));...
    ones(size(pfinf)),1*ones(size(pfinter)) ]);% bound: pfinter is 66(surf + tosolid force), pf1 is 1(farfield)
fclose(file);
