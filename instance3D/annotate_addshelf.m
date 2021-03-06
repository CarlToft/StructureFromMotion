function [trishelf,Ushelf]=annotate_addshelf(Utmp, shelfindex)


U0 = [0,0,0;80,0,0;80,0,-2;0,0,-2;...
          0,38.7,0;1.5,38.7,0;78.5,38.7,0;80,38.7,0;80,38.7,-2;0,38.7,-2;...
          0,2,0;1.5,2,0;0,2,77.2;1.5,2,77.2;...% 11-14
          0,38.7,77.2;1.5,38.7,77.2;... % 15-16
          78.5,38.7,77.2;80,38.7,77.2;... % 17-18
          78.5,2,77.2;80,2,77.2;...% 19-20
          0,0,77.2;80,0,77.2;...% 21-22
          78.5,2,0;80,2,0;...% 23-24
          1.5,3.2,36.2;78.5,3.2,36.2;78.5,38.7,36.2;1.5,38.7,36.2;...% 25-28
          1.5,3.2,34.2;78.5,3.2,34.2;78.5,38.7,34.2;1.5,38.7,34.2   % 29,32
          ]';

trishelf = [1,2,3;1,3,4;1,5,8;1,8,2;2,8,9;2,9,3;1,5,10;1,10,4;...
        11,12,14;11,14,13;12,6,16;12,16,14;...
        21,22,18;21,18,15;23,24,20;23,20,19;24,8,18;24,18,20;23,7,17;23,17,19;...
        29,30,26;29,26,25;25,26,27;25,27,28;29,30,31;29,31,32]';


T=findaff(structure(U0(:,shelfindex)),structure(Utmp));
Ushelf = T(1:3,:)*pextend(U0);

