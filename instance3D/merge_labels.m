function [labelobj, labeltype] = merge_labels(labelobjcell, labeltypecell)
    labelobj = zeros(1,0);
    labeltype = zeros(1,0);
    for j = 1:length(labelobjcell)
        labelobj = [labelobj,labelobjcell{j}];
        labeltype = [labeltype,labeltypecell{j}];
    end
end
