model = createpde('structural','static-solid');
importGeometry(model,'SC_2_0714_A1_F.STL');

%%
pdegplot(model,'FaceLabels','on')

%%
Froot = [1];
Ftip = [2];
    
%%
m1 = model.generateMesh('GeometricOrder','linear','Hmax',0.1);
pdeplot3D(model,'Mesh','on')

%% static load
structuralBC(model,'Face',Froot,'Constraint','fixed');
structuralBoundaryLoad(model,'Face',Ftip,'SurfaceTraction',[0;100e3;0]);
structuralProperties(model,'YoungsModulus',1e9,'PoissonsRatio',0.3);
model.SolverOptions.ReportStatistics = 'on';
res = model.solve();

%%
pdeplot3D(model,'ColorMapData',(res.Displacement.Magnitude),'Deformation',res.Displacement)

%%