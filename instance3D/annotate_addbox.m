function [tribox,Ubox]=annotate_addbox(Utmp, boxindex)


U0 = [0,0,0;1,0,0;1,1,0;0,1,0;0,0,1;1,0,1;1,1,1;0,1,1]';
tribox = [1,2,3;1,3,4;2,3,6;3,7,6;6,7,8;5,6,8;3,7,8;3,8,4;1,4,5;4,8,5]';

T=findaff(structure(U0(:,boxindex)),structure(Utmp));
Ubox = T(1:3,:)*pextend(U0);

