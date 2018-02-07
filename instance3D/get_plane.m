% Returns parameters for a particular plane in bookshelf
function [plane, front_direction] = get_plane(Ushelf, plane_idx)
    if plane_idx == 1
        pt1 = 25;
        pt2 = 26;
        pt3 = 27;
        pt4 = 12;
        pt5 = 14;
    else
        error('Only top plane (plane_idx==1) implemented');
    end

    v1 = Ushelf(:,pt2)-Ushelf(:,pt1);
    v2 = Ushelf(:,pt3)-Ushelf(:,pt2);
    n = cross(v1,v2);n= n/norm(n);

    plane = [n;-n'*Ushelf(:,pt1)];

    v3 = Ushelf(:,pt5)-Ushelf(:,pt4);
    v4 = Ushelf(:,pt3)-Ushelf(:,pt2);
    front_direction = cross(v3,v4);front_direction=front_direction/norm(front_direction);
end
