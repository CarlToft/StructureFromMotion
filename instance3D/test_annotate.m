
%OBJECT TYPES
%1 - cocacola pet 1.5 litre
%2 - cocacola light pet 1.5 litre
%3 - cocacola can 33 cl
%4 - blueberry soup
%5 - blackberry soup


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
    box_points_index = 7974;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN THE SCRIPT BELOW TO CREATE MESH and ANNOTATIONS
load C:\Users\fredrik\Documents\MATLAB\sfm_library\data/instanceRecognition/instance_sequences/experiment1/result_anno.mat


%figure(1);clf;trisurf(trishelf',Ushelf(1,:),Ushelf(2,:),Ushelf(3,:));axis equal;rotate3d on;

%construct mesh
nbrpoints_total = size(U,2);
nbrpoints_Utri = nbrpoints_total-nbrpoints_sfm;
Utri = zeros(3,0);
tri = zeros(3,0);

%add bookshelf 1
Utmp = U(1:3,nbrpoints_sfm+[1,3,18,19,2,21]);
shelfindex = [25,27,12,23,26,30];
[trishelf, Ushelf] = annotate_addshelf(Utmp, shelfindex);

tri = [tri,trishelf+size(Utri,2)];
Utri = [Utri,Ushelf];

%add bookshelf 2
Utmp = U(1:3,nbrpoints_sfm+[87,88,89,90,91,92]);
shelfindex = [26,27,23,24,12,1];
[trishelf, Ushelf] = annotate_addshelf(Utmp, shelfindex);

tri = [tri,trishelf+size(Utri,2)];
Utri = [Utri,Ushelf];

%add bookshelf 3
Utmp = U(1:3,nbrpoints_sfm+[92:97]);
shelfindex = [2,25,29,28,26,1];
[trishelf, Ushelf] = annotate_addshelf(Utmp, shelfindex);

tri = [tri,trishelf+size(Utri,2)];
Utri = [Utri,Ushelf];


%blueberry box
boxindex = [1,2,3,6,7,8];
[~, Ubox] = annotate_addbox(U(1:3,nbrpoints_sfm+[10:13,98,99]), boxindex);
tribox = [1,2,3;1,3,4;2,3,6;3,7,6;3,7,8;3,8,4;1,4,5;4,8,5]';

tri = [tri,tribox+size(Utri,2)];
Utri = [Utri,Ubox];


%add cocacola bottle cover
tribox = [8,1,4;1,4,5;1,2,5;2,5,6;2,3,6;3,6,7]';

tri = [tri,tribox+size(Utri,2)];
Utri = [Utri,U(1:3,nbrpoints_sfm+[22:28,18])];



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

%save C:\Users\fredrik\Documents\MATLAB\sfm_library\data/instanceRecognition/instance_sequences/experiment1/result_anno.mat settings P_uncalib u_uncalib U nbrpoints_sfm nbrpoints_anno Utri tri Uobj triobj labelobj labeltype

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
labeltype = zeros(1,0);

%plane 1 of bottles
v1 = Utri(:,26)-Utri(:,25);
v2 = Utri(:,27)-Utri(:,26);

n = cross(v1,v2);n= n/norm(n);
pl = [n;-n'*Utri(:,25)];


for ii=nbrpoints_anno+[1:5],
    labelcnt=labelcnt+1;
    labeltype(labelcnt) = 2; %cocacola light pet 1.5 litre
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
v1 = Utri(:,64+26)-Utri(:,64+25);
v2 = Utri(:,64+27)-Utri(:,64+26);

n = cross(v1,v2);n= n/norm(n);
pl = [n;-n'*Utri(:,64+25)];

for ii=nbrpoints_anno+[6:14,72],
    if ii==nbrpoints_anno+15,
        Upp1 = 2*U(1:3,nbrpoints_anno+14)-U(1:3,nbrpoints_anno+13);
        Upp2 = U(1:3,nbrpoints_anno+14)+U(1:3,nbrpoints_anno+9)-U(1:3,nbrpoints_anno+8);
        Upp = (Upp1+Upp2)/2;
    else
        Upp = U(1:3,ii);
    end
    labelcnt=labelcnt+1;
    labeltype(labelcnt) = 1; %cocacola pet 1.5 litre
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
v1 = Utri(:,2)-Utri(:,1);
v2 = Utri(:,7)-Utri(:,2);

n = cross(v1,v2);n= n/norm(n);
pl = [n;-n'*Utri(:,1)];

for ii=nbrpoints_anno+[15:38],
    Upp = U(1:3,ii);
    labelcnt=labelcnt+1;
    labeltype(labelcnt) = 3; %cocacola can 33 cl
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
         
Uframe = [0,0,0;7.5,0,0;7.5,7.0,1.5;0,7.0,1.5;0,0,-18.8;7.5,0,-18.8;7.5,7.0,-18.8;0,7.0,-18.8]';
trisaft(:,end+[1:12])=size(Usaft,2)+[1,2,3;1,3,4;5,6,7;5,7,8;1,2,5;2,5,6;2,3,7;2,7,6;3,4,7;4,7,8;1,5,8;1,8,4]';
Usaft = [Usaft,Uframe];
Ucc = Usaft(:,ccindex);

%figure(1);clf;trisurf(trisaft',Usaft(1,:),Usaft(2,:),Usaft(3,:));axis equal;hold on;rotate3d on;plot(structure(Ucc),'r*');

%plane of saft soppa
v1 = Utri(:,96+2)-Utri(:,96+1);
v2 = Utri(:,96+3)-Utri(:,96+2);
v3 = Utri(:,96+3)-Utri(:,96+2);
v4 = Utri(:,96+6)-Utri(:,96+2);
%v1 = Utri(:,32+26)-Utri(:,32+25);
%v2 = Utri(:,32+27)-Utri(:,32+26);
%v3 = Utri(:,32+27)-Utri(:,32+26);
%v4 = Utri(:,32+14)-Utri(:,32+12);


n = cross(v1,v2);n= n/norm(n);
n2 = cross(v3,v4);n2= n2/norm(n2);

%n = n + 0.05*n2;n=n/norm(n);
pl = [n;-n'*Utri(:,96+1)-0.0003];

cc = 0;
for ii=nbrpoints_anno+[39:58],
    Upp = U(1:3,ii);
%    if ii==nbrpoints_anno+50,
%        Upp=Upp-0.00008*n2;
%    elseif ii==nbrpoints_anno+54,
%        Upp=Upp-0.0001*n2;
%    elseif ii==nbrpoints_anno+57,
%        Upp=Upp-0.0001*n2;
    if ii==nbrpoints_anno+58,
        Upp=Upp-0.0001*n2;
    end
    labelcnt=labelcnt+1;
    cc = cc +1;
    if cc==5,
        cc=1;
    end
    if cc<3,
        labeltype(labelcnt) = 5; %blackberry
    else
        labeltype(labelcnt) = 4; %blueberry
    end
        
    ll = -n'*Upp-pl(4);
    Uplane = Upp+ll*n;

    sc = norm(Uplane-Upp)/norm(-18.8-Ucc(3));

    R = [null(n'), n]; if det(R)<0,R=[R(:,2),R(:,1),R(:,3)];end
    
    tmp = R'*n2;
    tmp(3)=0;tmp=tmp/norm(tmp);
    th2 = atan2(tmp(1),tmp(2))+3*pi/2;
    R = R*[cos(th2),sin(th2),0;-sin(th2),cos(th2),0;0,0,1];
    
    tt = Uplane - sc*R*[3.75;cos(th)*2.0;-18.8];
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


























