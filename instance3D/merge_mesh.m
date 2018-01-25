function [tri, Utri] = merge_mesh(mesh_data)
    tri = zeros(3,0);
    Utri = zeros(3,0);
    for j = 1:length(mesh_data)
        tri = [tri,mesh_data(j).tri+size(Utri,2)];
        Utri = [Utri,mesh_data(j).U];
    end
end
