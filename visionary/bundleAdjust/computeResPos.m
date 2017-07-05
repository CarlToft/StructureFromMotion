function res = computeResPos(P,C_p)
nRes = length(P)*3;
res = zeros(nRes,1);
iRes = 0;
for i = 1:length(P)
    if all(isfinite(C_p(:,i)))
        res(iRes+(1:3)) = P{i}*[C_p(:,i); 1];
        iRes = iRes + 3;
    end
end
