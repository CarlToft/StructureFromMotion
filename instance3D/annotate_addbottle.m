function [tri, U] = annotate_addbottle(plane, lid_position)
    [tribottle, Ubottle] = annotate_get_bottle_model();

    n = plane(1:3);

    ll = -n'*lid_position-plane(4);
    Uplane = lid_position+ll*n;

    sc = norm(Uplane-lid_position);

    R = [null(n'), n]; if det(R)<0,R=[R(:,2),R(:,1),R(:,3)];end
    T = [sc*R,Uplane;[0,0,0,1]];

    Utmp=T(1:3,:)*pextend(Ubottle);

    tri = tribottle;
    U = Utmp;
end
