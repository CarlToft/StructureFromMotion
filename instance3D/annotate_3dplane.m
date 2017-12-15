function [tri,Utri,Uvirtual]=annotate_3dplane(tri,Utri,index)

v12 = Utri(:,index(2))-Utri(:,index(1));
Uvirtual = Utri(:,index(3))-v12;

nbr=size(Utri,2);

Utri = [Utri,Uvirtual];
tri = [tri,index,[index([1,3]);nbr+1]];
