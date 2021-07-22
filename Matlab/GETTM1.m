%nozzle
model = createpde('structural','static-solid');
importGeometry(model,'NuzzleL.STL');

%%
pdegplot(model,'FaceLabels','on')

%%
Finf = [1,2,3];
Finf2 = [5];
    
%%
m1 = model.generateMesh('GeometricOrder','linear','Hmax',0.1,'Hmin',0.1);


%%
pdeplot3D(model,'Mesh','on')
%%

pfinf = model.Mesh.findNodes('region','Face',(Finf));
pfinf2 = model.Mesh.findNodes('region','Face',(Finf2));
pfwall = model.Mesh.findNodes('region','Face',setdiff(1:model.Geometry.NumFaces,[Finf,Finf2]));
pfall = model.Mesh.findNodes('region','Face',(1:model.Geometry.NumFaces));
pfinter = setdiff(pfall,pfinf);
pfinter = setdiff(pfinter,pfinf2);

[mp,me,mt] = model.Mesh.meshToPet;
np = size(mp,2); nt = size(mt,2);
plot3(mp(1,pfinter),mp(2,pfinter),mp(3,pfinter),'.')

model.Mesh.findElements('attached',mt(1:4,213))
%%
file = fopen('NZL_MeshCADLRD.txt','W');

fprintf(file,'NUM_POINT %d\n',np);
fprintf(file,'%.15e\t%.15e\t%.15e\n',mp);% pointcoord
fprintf(file,'NUM_TET %d\n',nt);
fprintf(file,'%d\t%d\t%d\t%d\n',mt(1:4,:)-1);% pointindex*4
fprintf(file,'NUM_BOUND %d\n',numel(pfall));
% pointindex-btype-bsetkey
fprintf(file,'%d\t%d\t%d\n',[pfinf-1,pfinf2-1,pfinter-1; ones(size(pfinf)),ones(size(pfinf2)),66*ones(size(pfinter));...
    ones(size(pfinf)),3*ones(size(pfinf2)),2*ones(size(pfinter)) ]);% bound: pfinter is 66(surf + tosolid force), pf1 is 1(farfield)
fclose(file);
