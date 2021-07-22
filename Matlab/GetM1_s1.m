GetM1;
%%
[mp,me,mt] = model.Mesh.meshToPet;
np = size(mp,2); nt = size(mt,2);
pf1 = model.Mesh.findNodes('region','Face',Froot);
pfall = model.Mesh.findNodes('region','Face',1:model.Geometry.NumFaces);
pfinter = setdiff(pfall,pf1);



plot3(mp(1,pf1),mp(2,pf1),mp(3,pf1),'.')
axis equal
%%
file = fopen('SC_2_0714_MeshFEM.txt','W');
fprintf(file,'NUM_POINT %d\n',np);
fprintf(file,'%.15e\t%.15e\t%.15e\n',mp);
fprintf(file,'NUM_TET %d\n',nt);
fprintf(file,'%d\t%d\t%d\t%d\n',mt(1:4,:)-1);
fprintf(file,'NUM_BOUND %d\n',numel(pfall));
%bound: pfinter is 66(fluid force), pf1 is 2(fixed)
fprintf(file,'%d\t%d\n',[pf1-1,pfinter-1; 2*ones(size(pf1)),66*ones(size(pfinter));...
    ones(size(pfinf)),2*ones(size(pfinter))]);
fclose(file);

%%
pitS = pfinter;
xitS = mp(:,pfinter);