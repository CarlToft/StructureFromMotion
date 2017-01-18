function p=getpluckercameras(P)
%MOTION/GETPLUCKERCAMERAS getpluckercameras(P) returns 
%cameras (with index) in plucker coordinates. 

% up row
ur = [2 3; 3 1; 1 2];
% down col
dc = [1 2; 1 3; 1 4; 3 4; 4 2; 2 3];

for i=1:size(P,2)
   p{i} = zeros(3,6);
   tt=P{i};
   for row = 1:3
      for col = 1:6
	 p{i}(row,col) = det([tt(ur(row,1),dc(col,1)) ... 
		    tt(ur(row,1),dc(col,2));
		    tt(ur(row,2),dc(col,1)) ... 
		    tt(ur(row,2),dc(col,2))]);
      end
   end
end

