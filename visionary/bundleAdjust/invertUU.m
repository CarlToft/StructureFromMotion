function invUU = invertUU(UU,fixedPoint)
if nargin < 2
    fixedPoint = false;
end
% if fixedPoint is true, the first coordinate of the first point is fixed
% by omitting it in the variables. Replace so that all blocks are 3x3
if fixedPoint
    [i,j,data] = find(UU);
    UU = sparse(i+1,j+1,data);
    UU(1,1) = 1;
end

%pick diagonal blocks
[i,j,data] = find(UU);
%much faster if I use full matrices here
A = full(sparse(i,mod(j-1,3)+1,data));
invUU = speye(size(UU,1));
[i,j,data] = find(invUU);
invA = full(sparse(i,mod(j-1,3)+1,data));

%Gauss elimination downwards
for i = 1:3
    %div by diagonal element
    invA(i:3:end,:) = invA(i:3:end,:)./repmat(A(i:3:end,i),[1 3]);
    A(i:3:end,:) = A(i:3:end,:)./repmat(A(i:3:end,i),[1 3]);
    for j = i+1:3
        %subtract element
        invA(j:3:end,:) = invA(j:3:end,:) - repmat(A(j:3:end,i),[1 3]).*invA(i:3:end,:);
        A(j:3:end,:) = A(j:3:end,:) - repmat(A(j:3:end,i),[1 3]).*A(i:3:end,:);
    end
end
%backsubstitution
for i = 3:-1:2
    for j = i-1:-1:1
        invA(j:3:end,:) = invA(j:3:end,:) - repmat(A(j:3:end,i),[1 3]).*invA(i:3:end,:);
        A(j:3:end,:) = A(j:3:end,:) - repmat(A(j:3:end,i),[1 3]).*A(i:3:end,:);
    end
end
[i,j,data] = find(invA);
if fixedPoint
    invUU = sparse(i(2:end)-1,j(2:end)+(ceil(i(2:end)/3)-1)*3-1,data(2:end));
else
    invUU = sparse(i,j+(ceil(i/3)-1)*3,data);
end
