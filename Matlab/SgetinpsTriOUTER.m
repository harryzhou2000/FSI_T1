%IMP TO FLUID
% dfile = fopen('SC_2_0714_OUTER-1.INP','r');
dfile = fopen('Triangle_I1_OUTER-K2.INP','r');
%% p,e get
lineA = '';
ps = nan(3,1000000);
es = nan(4,1000000);
while(~feof(dfile))
    line = fgetl(dfile);
    if(numel(strfind(line,'**')) > 0)
        continue;
    end
    %NODE
    if(numel(strfind(line,'*')) > 0)
        if(numel(strfind(line,'NODE')) > 0)
            while(1)
                lineA = fgetl(dfile);
                if(numel(strfind(lineA,'*'))>0)
                    break;
                end
                nums = split(lineA,',');
                ps(:,str2double(nums{1})) = [str2double(nums{2});str2double(nums{3});str2double(nums{4})];
                fprintf([nums{1},'\n']);
            end
            piscornor = false(1,size(ps,2));
        end
        if(numel(strfind(lineA,'ELEMENT')) > 0)
            while(1)
                lineA = fgetl(dfile);
                if(numel(strfind(lineA,'*'))>0)
                    break;
                end
                nums = split(lineA,',');
                es(:,str2double(nums{1})) = [str2double(nums{2});str2double(nums{3});str2double(nums{4});str2double(nums{5})];
                piscornor(es(:,str2double(nums{1}))) = 1;
                fprintf([nums{1},'\n']);
            end
        end
    end
end
%%
nps = sum(~isnan(ps),'all')/3;
nes = sum(~isnan(es),'all')/4;
ncorner = sum(piscornor);
pc = nan(3,ncorner);
pn2c = nan(1,nps);
mi = 0;
for i = 1:nps
    if(piscornor(i))
        mi = mi+1;
        pn2c(i) = mi;
        pc(:,mi) = ps(:,i);
    else
        pn2c(i) = -1;
    end
end
ec = nan(size(es));
for i = 1:nes
    ec(:,i) = pn2c(es(:,i))';
end
%%
fclose(dfile);

%% region get
% dfile = fopen('SC_2_0714_OUTER-bound-.INP','r');
 dfile = fopen('Triangle_I1_OUTER-bound-K2.INP','r');
nsets = 0;
sets = {};
hasend = 0;
while(~feof(dfile))
    if(hasend)
        line = lineA;
        hasend = 0;
    else
        line = fgetl(dfile);
        
    end
    
    if(numel(strfind(line,'NSET,NSET=BC')) > 0)
        nsets = nsets+1;
        sets{end+1} = nan(1,1000);
        iset = 1;
        while(1)
            lineA = fgetl(dfile);
            if(numel(strfind(lineA,'*'))>0)
                break;
            end
            nums = split(lineA,',');
            for it = 1:numel(nums)
                if((~isnan(str2double(nums{it}))) && pn2c(str2double(nums{it})) > 0)
                    sets{end}(iset) = pn2c(str2double(nums{it}));
                    iset = iset +1;
                end
            end
        end
        sets{end} = sets{end}(1:(sum(~isnan(sets{end}))));
        hasend = 1;
    end
    
    
end
fclose(dfile);
%%
figure(1);
clf; hold on;
for i = 1:numel(sets)
    if sum(i == [4 7 8 10])
        plot3(pc(1,sets{i}),pc(2,sets{i}),pc(3,sets{i}),'o');
    else
        plot3(pc(1,sets{i}),pc(2,sets{i}),pc(3,sets{i}),'.');
    end
    
end
axis equal;
legend


%% NOINFFIRST
Finf = [1 2 3 5 6 9];
Finf2 = [];
pfall = []; pfinf = []; pfinf2 = [];
for i = 1:numel(sets)
    pfall = union(sets{i}, pfall);
end
for i = Finf
    pfinf = union(sets{i}, pfinf);
end
for i = Finf2
    pfinf2 = union(sets{i}, pfinf2);
end
pfinter = setdiff(pfall,pfinf);
pfinter = setdiff(pfinter,pfinf2);
figure(1); clf; hold on;
plot3(pc(1,pfinter),pc(2,pfinter),pc(3,pfinter),'.');
plot3(pc(1,pfinf2),pc(2,pfinf2),pc(3,pfinf2));
plot3(pc(1,pfinf),pc(2,pfinf),pc(3,pfinf),'.');
set(gca,'Clipping','off');

mt = ec(:,~isnan(sum(ec,1))); mp = pc;
np = size(mp,2); nt = size(mt,2);
pfinf = pfinf'; pfinf2 = pfinf2'; pfinter = pfinter';
%% from GETTM1.m
file = fopen('MeshCAD_TRI_K2.txt','W');

fprintf(file,'NUM_POINT %d\n',np);
fprintf(file,'%.15e\t%.15e\t%.15e\n',mp);% pointcoord
fprintf(file,'NUM_TET %d\n',nt);
fprintf(file,'%d\t%d\t%d\t%d\n',mt(1:4,:)-1);% pointindex*4
fprintf(file,'NUM_BOUND %d\n',numel(pfall));
% pointindex-btype-bsetkey
fprintf(file,'%d\t%d\t%d\n',[pfinf-1,pfinf2-1,pfinter-1; ones(size(pfinf)),3*ones(size(pfinf2)),66*ones(size(pfinter));...
    ones(size(pfinf)),3*ones(size(pfinf2)),2*ones(size(pfinter)) ]);% bound: pfinter is 66(surf + tosolid force), pf1 is 1(farfield)
fclose(file);