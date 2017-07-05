function res = computeResRot(P,R_p)
res = zeros(9*length(P),1);
for i = 1:length(P)
    if ~isempty(R_p{i})
        dR = P{i}(1:3,1:3) - R_p{i};
        res((i-1)*9+(1:9)) = dR(:);
    end
end