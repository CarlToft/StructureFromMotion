function [p,pc]=getpoints(i,index)
%IMAGEDATA/GETPOINTS [p,pc]=getpoints(i,index) returns i's points p with index
%  in homog. coordinates. and covariances pc

if nargin==1,
  p=i.points;
  pc=i.pointcov;
else
  p=i.points(:,index);
  if nargout==2,
    if isempty(i.pointcov),
      %covariances are not set
      if size(index,2)==1,
        pc = [1 0 0;0 1 0];
      else
        pc = cell(1,size(index,2));
        for t=1:size(index,2);
          pc{i} = [1 0 0;0 1 0];
        end
      end
    else %covariances are set
      pc=i.pointcov(index);
      if size(pc,2)==1;
        pc=pc{1};
      end
    end
  end
end



