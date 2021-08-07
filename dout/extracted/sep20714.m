function [Tu,Tl,Tr] =sep20714(T)
Tu = T(T.y>0.012,:);
Tl =  T(T.y<=0.012,:);
Tu = sortrows(Tu,'x');
Tl = sortrows(Tl,'x','descend');
Tr = [Tu;Tl];