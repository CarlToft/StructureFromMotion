function [tri, U] = annotate_addsaftsoppa(plane, lid_position, front_direction)
    if size(front_direction, 2) == 3
        % Make sure we are dealing with column vectors
        front_direction = front_direction.';
    end
    if size(front_direction, 2) == 2
        % Two reference points whose difference defines the alignment vector
        front_direction = front_direction(:,2) - front_direction(:,1);
    end
    front_direction = front_direction / norm(front_direction);


    [trisaft, Usaft, Ucc, th] = annotate_get_saftsoppa_model();

    n = plane(1:3);

    % Project lid position onto the plane of the shelf
    ll = -n'*lid_position-plane(4);
    Uplane = lid_position+ll*n;

    sc = norm(Uplane-lid_position)/norm(-18.8-Ucc(3));

    R = [null(n'), n]; if det(R)<0,R=[R(:,2),R(:,1),R(:,3)];end
    n2ii=cross(front_direction,n);n2ii=n2ii/norm(n2ii);
    tmp = R'*n2ii;
    tmp(3)=0;tmp=tmp/norm(tmp);
    th2 = atan2(tmp(1),tmp(2));

    R = R*[cos(th2),sin(th2),0;-sin(th2),cos(th2),0;0,0,1];

    tt = Uplane - sc*R*[3.75;cos(th)*2.0;-18.8];
    T = [sc*R,tt;[0,0,0,1]];

    Utmp=T(1:3,:)*pextend(Usaft);

    tri = trisaft;
    U = Utmp;
end
