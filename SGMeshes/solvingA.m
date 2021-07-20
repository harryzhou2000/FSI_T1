A = getA ('PE2_A.txt');
A = sparse(A(:,1)+1,A(:,2)+1,A(:,3));
b = getb('PE2_b.txt');
b = b.VarName1;
M = getA('PE2_M.txt');
M = sparse(M(:,1)+1,M(:,2)+1,M(:,3));
% ATA = full(A'*A);
% ATb = full(A'*b);
Ap = A(sum(A,2)~=0,sum(A,2)~=0);
bp = b(sum(A,2)~=0,:);
ps = symamd(Ap);
As = Ap(ps,ps);
bs = bp(ps);
dm = diag(M);
nzero = dm~=0;

Anz = A(nzero,nzero);
Mnz = M(nzero,nzero);
eigsAA =eigs(Anz,Mnz,10,'smallestabs');
MnzRatio = max(abs(Mnz),[],'all');
Mnz = Mnz/MnzRatio;
% dg = max(Ap,[],2);
% dg = diag(Ap);
% S = sqrt(diag(dg));
% Si = inv(S);


[V,D] = eigs(Anz+Mnz*1e-5,Mnz,10,'smallestabs');
D1 = full((D-eye(size(D)))/MnzRatio);
Anz = 0.5*(Anz+Anz');
%%
evenI = 2:2:size(Anz,1);
Anz = Anz(evenI,evenI);
%%
Ga = graph(Anz- diag(diag(Anz)));
Ga.plot('Layout','Force','Iterations',1000,'WeightEffect','none','UseGravity','off');
%%
% L =  tril(Ap);
% D = diag(diag(Ap));
% omega = 0.7;
% S = 1/sqrt(omega*(2-omega))*(D-omega*L)*diag(1./sqrt(diag(Ap)));
% Si = diag(sqrt(diag(Ap)))*inv((D-omega*L))*sqrt(omega*(2-omega));
% cond(Si*Ap*Si')

%%
[P,R,C] = equilibrate(Ap);
cond(R*P*Ap*C)