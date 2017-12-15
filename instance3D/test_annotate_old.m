


bib='C:\Users\fredrik\Documents\MATLAB\sfm_library\data\InstanceRecognition\instance_sequences\';

file='experiment1\';

if 0,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
%load settings and result file
load([bib,file,'result.mat']);
nbrpoints_sfm = size(U,2);

mot = motion(P_uncalib);
figure(1);clf;plot(mot);axis equal;


%%%%%%%%%%%%%%%%%%%

%now we add points defining planes

[u_uncalib,U]=add_1point(settings,P_uncalib,U,u_uncalib,[1,4,44,48]);
[u_uncalib,U]=add_1point(settings,P_uncalib,U,u_uncalib,[1,4,44,48]);
[u_uncalib,U]=add_1point(settings,P_uncalib,U,u_uncalib,[1,4,44,48]);

[U,P_uncalib,lambda] = modbundle_sparse(U,P_uncalib,u_uncalib,10,10);

%save C:\Users\fredrik\Documents\MATLAB\sfm_library\data/instanceRecognition/instance_sequences/experiment1/result_anno.mat settings P_uncalib u_uncalib U nbrpoints_sfm

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN THE SCRIPT BELOW TO CREATE MESH and ANNOTATIONS
load C:\Users\fredrik\Documents\MATLAB\sfm_library\data/instanceRecognition/instance_sequences/experiment1/result_anno.mat



%construct mesh
nbrpoints_total = size(U,2);
nbrpoints_Utri = nbrpoints_total-nbrpoints_sfm;
Utri = zeros(3,0);
tri = zeros(3,0);

boxindex = [1,2,3,6];
[tribox, Ubox] = annotate_addbox(U(1:3,nbrpoints_sfm+[1:4]), boxindex);

tri = [tri,tribox+size(Utri,2)];
Utri = [Utri,Ubox];

boxindex = [1,2,3,6];
[tribox, Ubox] = annotate_addbox(U(1:3,nbrpoints_sfm+[6:9]), boxindex);

tri = [tri,tribox+size(Utri,2)];
Utri = [Utri,Ubox];

boxindex = [1,2,3,6];
[~, Ubox] = annotate_addbox(U(1:3,nbrpoints_sfm+[10:13]), boxindex);
tribox = [1,2,3;1,3,4;2,3,6;3,7,6;3,7,8;3,8,4;1,4,5;4,8,5]';

tri = [tri,tribox+size(Utri,2)];
Utri = [Utri,Ubox];

boxindex = [4,1,2,6];
[tribox, Ubox] = annotate_addbox(U(1:3,nbrpoints_sfm+[14:17]), boxindex);

tri = [tri,tribox+size(Utri,2)];
Utri = [Utri,Ubox];

boxindex = [1,2,3,6];
[tribox, Ubox] = annotate_addbox(U(1:3,nbrpoints_sfm+[18:21]), boxindex);

tri = [tri,tribox+size(Utri,2)];
Utri = [Utri,Ubox];

tribox = [-7,1,4;1,4,5;1,2,5;2,5,6;2,3,6;3,6,7]';

tri = [tri,tribox+size(Utri,2)];
Utri = [Utri,U(1:3,nbrpoints_sfm+[22:28])];



%figure(1);clf;trisurf(tri',Utri(1,:),Utri(2,:),Utri(3,:));axis equal;rotate3d on;

%for imindex = 1:48,imindex
%    annotate_plot(settings,P_uncalib,imindex, tri, Utri);
%    pause
%    close(imindex);
%end

nbrpoints_anno = nbrpoints_sfm+28;

%now we add points defining objects

%[u_uncalib,U]=add_1point(settings,P_uncalib,U,u_uncalib,[1,4,44,48]);

%[U,P_uncalib,lambda] = modbundle_sparse(U,P_uncalib,u_uncalib,10,10);

%save C:\Users\fredrik\Documents\MATLAB\sfm_library\data/instanceRecognition/instance_sequences/experiment1/result_anno.mat settings P_uncalib u_uncalib U nbrpoints_sfm nbrpoints_anno

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ADD OBJECTS
%bottles


alpha = 0:0.4:2*pi;rr = 0.13;height=0.8;nn=length(alpha);
Ubottle=[rr*cos(alpha),rr*cos(alpha),0;rr*sin(alpha),rr*sin(alpha),0;zeros(size(alpha)),height*ones(size(alpha)),1];
tribottle = [1:nn,      nn+1:2*nn, 2*nn+1*ones(1,nn) ; ...
             [2:nn,1],  [nn+2:2*nn,nn+1]  , nn+1:2*nn, ;...
             nn+1:2*nn, [2:nn,1]       , [nn+2:2*nn,nn+1]];
%figure(1);clf;trisurf(tribottle',Ubottle(1,:),Ubottle(2,:),Ubottle(3,:));axis equal;rotate3d on;

%labels
labelcnt=0;

%gravity vector

Uobj = zeros(3,0);
triobj = zeros(3,0);
labelobj = zeros(1,0);

%plane 1 of bottles
v1 = Utri(:,2)-Utri(:,1);
v2 = Utri(:,3)-Utri(:,2);

n = cross(v1,v2);n= n/norm(n);
pl = [n;-n'*Utri(:,1)];

for ii=nbrpoints_anno+[1:5],
    labelcnt=labelcnt+1;
    ll = -n'*U(1:3,ii)-pl(4);
    Uplane = U(1:3,ii)+ll*n;

    sc = norm(Uplane-U(1:3,ii));

    R = [null(n'), n]; if det(R)<0,R=[R(:,2),R(:,1),R(:,3)];end
    T = [sc*R,Uplane;[0,0,0,1]];

    Utmp=T(1:3,:)*pextend(Ubottle);

    triobj = [triobj,tribottle+size(Uobj,2)];
    Uobj=[Uobj,Utmp];
    labelobj=[labelobj,labelcnt*ones(1,size(tribottle,2))];
end

%plane 2 of bottles
v1 = Utri(:,26)-Utri(:,25);
v2 = Utri(:,27)-Utri(:,25);

n = cross(v1,v2);n= n/norm(n);
pl = [n;-n'*Utri(:,25)];

for ii=nbrpoints_anno+[6:14],
    if ii==nbrpoints_anno+15,
        Upp1 = 2*U(1:3,nbrpoints_anno+14)-U(1:3,nbrpoints_anno+13);
        Upp2 = U(1:3,nbrpoints_anno+14)+U(1:3,nbrpoints_anno+9)-U(1:3,nbrpoints_anno+8);
        Upp = (Upp1+Upp2)/2;
    else
        Upp = U(1:3,ii);
    end
    labelcnt=labelcnt+1;
    ll = -n'*Upp-pl(4);
    Uplane = Upp+ll*n;

    sc = norm(Uplane-Upp);

    R = [null(n'), n]; if det(R)<0,R=[R(:,2),R(:,1),R(:,3)];end
    T = [sc*R,Uplane;[0,0,0,1]];

    Utmp=T(1:3,:)*pextend(Ubottle);

    triobj = [triobj,tribottle+size(Uobj,2)];
    Uobj=[Uobj,Utmp];
    labelobj=[labelobj,labelcnt*ones(1,size(tribottle,2))];
end

alpha = 0:0.4:2*pi;rr = 0.24;height=1;nn=length(alpha);
Ubottle=[rr*cos(alpha),rr*cos(alpha),0;rr*sin(alpha),rr*sin(alpha),0;zeros(size(alpha)),height*ones(size(alpha)),1];
tribottle = [1:nn,      nn+1:2*nn, 2*nn+1*ones(1,nn) ; ...
             [2:nn,1],  [nn+2:2*nn,nn+1]  , nn+1:2*nn, ;...
             nn+1:2*nn, [2:nn,1]       , [nn+2:2*nn,nn+1]];
%figure(1);clf;trisurf(tribottle',Ubottle(1,:),Ubottle(2,:),Ubottle(3,:));axis equal;rotate3d on;


%plane of coke bottles
v1 = Utri(:,34)-Utri(:,33);
v2 = Utri(:,35)-Utri(:,33);

n = cross(v1,v2);n= n/norm(n);
pl = [n;-n'*Utri(:,33)];

for ii=nbrpoints_anno+[15:38],
    Upp = U(1:3,ii);
    labelcnt=labelcnt+1;
    ll = -n'*Upp-pl(4);
    Uplane = Upp+ll*n;

    sc = norm(Uplane-Upp);

    R = [null(n'), n]; if det(R)<0,R=[R(:,2),R(:,1),R(:,3)];end
    T = [sc*R,Uplane;[0,0,0,1]];

    Utmp=T(1:3,:)*pextend(Ubottle);

    triobj = [triobj,tribottle+size(Uobj,2)];
    Uobj=[Uobj,Utmp];
    labelobj=[labelobj,labelcnt*ones(1,size(tribottle,2))];
end

%saft soppa
alpha = 0:0.4:2*pi;rr = 1.5;height=1;nn=length(alpha);
cc = [3.75,2.0];

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
         
Uframe = [0,0,0;7.5,0,0;7.5,7.0,1.5;0,7.0,1.5;0,0,-18;7.5,0,-18;7.5,7.0,-18;0,7.0,-18]';
trisaft(:,end+[1:12])=size(Usaft,2)+[1,2,3;1,3,4;5,6,7;5,7,8;1,2,5;2,5,6;2,3,7;2,7,6;3,4,7;4,7,8;1,5,8;1,8,4]';
Usaft = [Usaft,Uframe];
Ucc = Usaft(:,ccindex);

%figure(1);clf;trisurf(trisaft',Usaft(1,:),Usaft(2,:),Usaft(3,:));axis equal;hold on;rotate3d on;plot(structure(Ucc),'r*');

%plane of saft soppa
%v1 = Utri(:,18)-Utri(:,17);
%v2 = Utri(:,19)-Utri(:,17);
v1 = Utri(:,10)-Utri(:,9);
v2 = Utri(:,11)-Utri(:,9);
v3 = Utri(:,10)-Utri(:,9);
v4 = Utri(:,13)-Utri(:,9);


n = cross(v1,v2);n= n/norm(n);
n2 = cross(v3,v4);n2= n2/norm(n2);
pl = [n;-n'*Utri(:,9)];

for ii=nbrpoints_anno+[39:58],
    Upp = U(1:3,ii);
    labelcnt=labelcnt+1;
    ll = -n'*Upp-pl(4);
    Uplane = Upp+ll*n;

    sc = norm(Uplane-Upp)/norm(-18-Ucc(3));

    R = [null(n'), n]; if det(R)<0,R=[R(:,2),R(:,1),R(:,3)];end
    
    tmp = R'*n2;
    tmp(3)=0;tmp=tmp/norm(tmp);
    th2 = atan2(tmp(1),tmp(2))+pi;
    R = R*[cos(th2),sin(th2),0;-sin(th2),cos(th2),0;0,0,1];
    
    tt = Uplane - sc*R*[3.75;cos(th)*2.0;-18];
    T = [sc*R,tt;[0,0,0,1]];

    Utmp=T(1:3,:)*pextend(Usaft);

    triobj = [triobj,trisaft+size(Uobj,2)];
    Uobj=[Uobj,Utmp];
    labelobj=[labelobj,labelcnt*ones(1,size(trisaft,2))];
end







figure(1);clf;trisurf(tri',Utri(1,:),Utri(2,:),Utri(3,:));axis equal;rotate3d on;hold on;
trisurf(triobj',Uobj(1,:),Uobj(2,:),Uobj(3,:));


pause;
         
         

for imindex = 1:48,imindex
    annotate_plot(settings,P_uncalib,imindex, tri, Utri, triobj, Uobj, labelobj);
    pause
    close(imindex);
end


























