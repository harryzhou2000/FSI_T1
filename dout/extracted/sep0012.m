function [Tu,Tl,Tr] =sep0012(T)
Tu = T(T.y>0,:);
Tl =  T(T.y<=0,:);
Tu = sortrows(Tu,'x');
Tl = sortrows(Tl,'x','descend');
Tr = [Tu;Tl];