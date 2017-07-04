function res = computeResCam(P,U,u)
nRes = 0;
for i = 1:length(P)
    nRes = nRes + length(u.index);
end
res = zeros(nRes,1);
iRes = 0;
for i = 1:length(P)
    uu = NaN*ones(3,u.pointnr);
    uu(:,u.index{i}) = pflat(u.points{i});
    vis = isfinite(uu(1,:));
    r = [ ...
        ((P{i}(1,:)*U(:,vis))./(P{i}(3,:)*U(:,vis)) - uu(1,vis)); ...
        ((P{i}(2,:)*U(:,vis))./(P{i}(3,:)*U(:,vis)) - uu(2,vis))];
    res(iRes+(1:numel(r))) = r(:);
    iRes = iRes + numel(r);
end
