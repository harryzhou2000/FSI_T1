%% 0012 T1
clear;
exp0012;
DO2 = ImportExtracted('SG_0012_T1_O2.dat');
DO1 = ImportExtracted('SG_0012_T1_O1.dat');
DRK = ImportExtracted('SG_0012_T1_RK4.dat');

DO2c = DO2(DO2.x>=-2 & DO2.x<=-1,:);
DO1c = DO1(DO1.x>=-2 & DO1.x<=-1,:);
DRKc = DRK(DRK.x>=-2 & DRK.x<=-1,:);
[~,~,DO2r] = sep0012(DO2c);
[~,~,DO1r] = sep0012(DO1c);
[~,~,DRKr] = sep0012(DRKc);
clf; 
qinf = 0.775^2 * 1.4 /2;
subplot(1,2,1);
hold on;   
plot(DRKr.x+2,(DRKr.p-1)/qinf,'DisplayName','Order 2');
% plot(DO2r.x+2,(DO2r.p-1)/qinf,'DisplayName','Order 2 Euler');
plot(DO1r.x+2,(DO1r.p-1)/qinf,'DisplayName','Order 1');
plot(EXPup(:,1),EXPup(:,2),'o','DisplayName','Re 9.9e6 Experiment Upper');
plot(EXPlo(:,1),EXPlo(:,2),'o', 'DisplayName','Re 9.9e6 Experiment Lower');

hold off;
legend; xlabel('x');ylabel('Cp');
grid on; grid minor;

subplot(1,2,2);
hold on;
plot(DRKr.x+2,DRKr.ma,'DisplayName','Order 2');
% plot(DO2r.x+2,(DO2r.p-1)/qinf,'DisplayName','Order 2 Euler');
plot(DO1r.x+2,DO1r.ma,'DisplayName','Order 1');
plot(EXPup(:,1),EXPup(:,4),'o','DisplayName','Re 9.9e6 Experiment Upper');
plot(EXPlo(:,1),EXPlo(:,4),'o', 'DisplayName','Re 9.9e6 Experiment Lower');

hold off;
legend; xlabel('x');ylabel('Ma');
grid on; grid minor;


%% 20714
D714 = ImportExtracted('SG_20714_FaceTrace.dat');
D714c = D714(D714.x>=-1 & D714.x<=-0,:);
[~,~,D714r] = sep20714(D714c);

subplot(1,2,1)
plot(D714r.x + 1,(D714r.p - 1)/qinf);
xlabel('x');ylabel('Cp');
grid on; grid minor;
subplot(1,2,2)
plot(D714r.x + 1,D714r.ma);
xlabel('x');ylabel('Ma');
grid on; grid minor;
%% 0012 T2

DM04 = ImportExtracted('SG_0012_T2_A125_M04.dat');
DM08 = ImportExtracted('SG_0012_T2_A125_M08.dat');
DM04c = DM04(DM04.x>=-2 & DM04.x<=-1,:);
DM08c = DM08(DM08.x>=-2 & DM08.x<=-1,:);
[~,~,DM04r] = sep0012(DM04c);
[~,~,DM08r] = sep0012(DM08c);
clf; 
subplot(1,2,1);
hold on;   
plot(DM04r.x+2,DM04r.p,'DisplayName','Ma 0.4');
plot(DM08r.x+2,DM08r.p,'DisplayName','Ma 0.8');
hold off;
legend; xlabel('x');ylabel('p');grid on; grid minor;
subplot(1,2,2);
hold on;   
plot(DM04r.x+2,DM04r.ma,'DisplayName','Ma 0.4');
plot(DM08r.x+2,DM08r.ma,'DisplayName','Ma 0.8');
hold off;
legend; xlabel('x');ylabel('Ma');grid on; grid minor;


