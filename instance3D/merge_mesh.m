function [tri, Utri] = merge_mesh(tricell, Utricell)
    tri = zeros(3,0);
    Utri = zeros(3,0);
    for j = 1:length(tricell)
        tri = [tri,tricell{j}+size(Utri,2)];
        Utri = [Utri,Utricell{j}];
    end
end
