% Returns parameters for a particular plane in bookshelf
function [plane] = get_plane(Ushelf, plane_idx)
    if plane_idx == 1
        pt1 = 25;
        pt2 = 26;
        pt3 = 27;
    else
        error('Only top plane (plane_idx==1) implemented');
    end

    v1 = Ushelf(:,pt2)-Ushelf(:,pt1);
    v2 = Ushelf(:,pt3)-Ushelf(:,pt2);
    n = cross(v1,v2);n= n/norm(n);

    plane = [n;-n'*Ushelf(:,pt1)];
end
