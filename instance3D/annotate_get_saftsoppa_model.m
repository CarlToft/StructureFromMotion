function [trisaft, Usaft, Ucc, th] = annotate_get_saftsoppa_model()
    %saft soppa
    alpha = 0:0.4:2*pi;rr = 1.5;height=1;nn=length(alpha);
    % Center of lid in local frame
    cc = [3.75,2.0];

    % Slope of upper plane
    th = atan(1.5/7.0);
    Usaft=[rr*cos(alpha),rr*cos(alpha),0;rr*sin(alpha),rr*sin(alpha),0;zeros(size(alpha)),height*ones(size(alpha)),1];
    ccindex = size(Usaft,2);

    Usaft(1,:)=Usaft(1,:)+cc(1);
    Usaft(2,:)=Usaft(2,:)+cc(2);

    R = [1,0,0;0,cos(th),sin(th);0,-sin(th),cos(th)];

    Usaft = R'*Usaft;
    trisaft = [1:nn,      nn+1:2*nn, 2*nn+1*ones(1,nn) ; ...
                 [2:nn,1],  [nn+2:2*nn,nn+1]  , nn+1:2*nn, ;...
                 nn+1:2*nn, [2:nn,1]       , [nn+2:2*nn,nn+1]];

    Uframe = [0,0,0;7.5,0,0;7.5,7.0,1.5;0,7.0,1.5;0,0,-18.8;7.5,0,-18.8;7.5,7.0,-18.8;0,7.0,-18.8]';
    trisaft(:,end+[1:12])=size(Usaft,2)+[1,2,3;1,3,4;5,6,7;5,7,8;1,2,5;2,5,6;2,3,7;2,7,6;3,4,7;4,7,8;1,5,8;1,8,4]';
    Usaft = [Usaft,Uframe];
    Ucc = Usaft(:,ccindex);
end
