function [instance_labels, class_labels] = merge_labels(instance_labels_cell, class_labels_cell)
    instance_labels = zeros(1,0);
    class_labels = zeros(1,0);
    for j = 1:length(instance_labels_cell)
        instance_labels = [instance_labels,instance_labels_cell{j}];
        class_labels = [class_labels,class_labels_cell{j}];
    end
end
