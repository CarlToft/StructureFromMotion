function [tri, Utri, instance_labels] = merge_mesh(mesh_data)
    tri = zeros(3,0);
    Utri = zeros(3,0);
    instance_labels = zeros(1,0);
    for j = 1:length(mesh_data)
        tri = [tri,mesh_data(j).tri+size(Utri,2)];
        Utri = [Utri,mesh_data(j).U];
        instance_labels = [instance_labels, j*ones(1, length(mesh_data(j).tri))];
    end
end
