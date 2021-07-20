x1 = [0,0];
x2 = [0,1];
x3 = [1,0];
x4 = [1,1];
x5 = [0.5,0.3];
x6 = [0.5,0.7];
a1 = 0.25;
a2 = 0.4*0.25;
a3 = 0.5*0.3;
a4 = a2;
a5 = a1;
a6 = a3;

p1p1 = -dot(x1-x2,x1-x2)/a1/4 ...
    -dot(x6-x2,x6-x2)/a2/4 ...
    -dot(x6-x4,x6-x4)/a4/4 ...
    -dot(x3-x4,x3-x4)/a5/4 ...
    -dot(x1-x3,x1-x3)/a6/4;
p2p2 = -dot(x2-x5,x2-x5)/a2/4 ...
     -dot(x5-x4,x5-x4)/a4/4 ...
     -dot(x4-x2,x4-x2)/a3/4;

     

   


A = [p1p1,2.3;2.3,p2p2];
b = [0.666666;1.2666666];
u0  = [0;0];
r0 = b-A*u0;