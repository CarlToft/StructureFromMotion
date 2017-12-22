
%OBJECT TYPES
%1 - cocacola pet 1.5 litre
%2 - cocacola light pet 1.5 litre
%3 - cocacola can 33 cl
%4 - blueberry soup
%5 - blackberry soup


%bib='C:\Users\fredrik\Documents\MATLAB\sfm_library\data\InstanceRecognition\instance_sequences\';
bib='C:\Users\fredrik\Dropbox\tmp\Lucas\InstanceRecognition\instance_sequences\';

file='GOPR1830\';

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
    
    %save C:\Users\fredrik\Documents\MATLAB\sfm_library\data/instanceRecognition/instance_sequences/experiment2/result_anno.mat settings P_uncalib u_uncalib U nbrpoints_sfm

end
if true
%now we add points defining planes

%     4 Regular Coca-Cola bottles
    [u_uncalib,U]=add_1point(settings,P_uncalib,U,u_uncalib,[60,68,44,48]);
    [u_uncalib,U]=add_1point(settings,P_uncalib,U,u_uncalib,[60,68,44,48]);
    [u_uncalib,U]=add_1point(settings,P_uncalib,U,u_uncalib,[60,68,44,48]);
    [u_uncalib,U]=add_1point(settings,P_uncalib,U,u_uncalib,[60,68,44,48]);

%     3 Coca-Cola Zero bottles
    [u_uncalib,U]=add_1point(settings,P_uncalib,U,u_uncalib,[60,68,44,48]);
    [u_uncalib,U]=add_1point(settings,P_uncalib,U,u_uncalib,[60,68,44,48]);
    [u_uncalib,U]=add_1point(settings,P_uncalib,U,u_uncalib,[60,68,44,48]);
    
    [U,P_uncalib,lambda] = modbundle_sparse(U,P_uncalib,u_uncalib,10,10);
    
    %save C:\Users\fredrik\Documents\MATLAB\sfm_library\data/instanceRecognition/instance_sequences/experiment2/result_anno.mat settings P_uncalib u_uncalib U nbrpoints_sfm

end
STOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN THE SCRIPT BELOW TO CREATE MESH and ANNOTATIONS
eval(['load ',bib,file,'result_anno.mat']);


%figure(1);clf;trisurf(trishelf',Ushelf(1,:),Ushelf(2,:),Ushelf(3,:));axis equal;rotate3d on;

%construct mesh
nbrpoints_total = size(U,2);
nbrpoints_Utri = nbrpoints_total-nbrpoints_sfm;
Utri = zeros(3,0);
tri = zeros(3,0);

%add bookshelf 1
Utmp = U(1:3,nbrpoints_sfm+[1:6]);
shelfindex = [1,22,26,27,25,21];
[trishelf, Ushelf] = annotate_addshelf(Utmp, shelfindex);

tri = [tri,trishelf+size(Utri,2)];
Utri = [Utri,Ushelf];


%add bookshelf 2
Utmp = U(1:3,nbrpoints_sfm+[7:12]);
shelfindex = [22,2,1,26,25,7];
[trishelf, Ushelf] = annotate_addshelf(Utmp, shelfindex);

tri = [tri,trishelf+size(Utri,2)];
Utri = [Utri,Ushelf];

%add bookshelf 3
Utmp = U(1:3,nbrpoints_sfm+[13:17]);
shelfindex = [2,6,1,25,22];
[trishelf, Ushelf] = annotate_addshelf(Utmp, shelfindex);

tri = [tri,trishelf+size(Utri,2)];
Utri = [Utri,Ushelf];



% %blue berry covers
% tribox = [1,2,3;1,3,4;3,4,7;4,7,8;5,6,7;5,7,8;5,8,9;10,11,12;11,13,12;7,8,12;8,12,13]';
% 
% Utmp = U(1:3,nbrpoints_sfm+[20:31]);
% %small fix of one point
% v1=U(1:3,nbrpoints_sfm + 29)-U(1:3,nbrpoints_sfm + 30);v2=U(1:3,nbrpoints_sfm + 31)-U(1:3,nbrpoints_sfm + 30);n=cross(v1,v2);n=n/norm(n);
% Utmp(:,11)=Utmp(:,11)+0.00008*n;
% 
% 
% 
% 
% v1 = Utmp(:,11)-Utmp(:,10);
% v2 = Utmp(:,12)-Utmp(:,10);
% Utmp = [Utmp,(Utmp(:,11)+v2+Utmp(:,12)+v1)/2];
% tri = [tri,tribox+size(Utri,2)];
% Utri = [Utri,Utmp];
% 
% %coca cola cans cover
% tribox = [1,2,3;1,3,4;3,4,6;4,5,6;5,6,7;5,7,8]';
% 
% Utmp = U(1:3,nbrpoints_sfm+[32:39]);
% tri = [tri,tribox+size(Utri,2)];
% Utri = [Utri,Utmp];






% ADD OBJECTS
%bottles

nbrpoints_anno = nbrpoints_sfm+18;

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


labeltype = [labeltype,2,1,2,2,1,1];

% 4 Regular Coca-Cola bottles
for ii=nbrpoints_anno+[1:4],
    labelcnt = labelcnt+1;
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

% 3 Coca-Cola Zero bottles
for ii=nbrpoints_anno+[5:7],
    labelcnt = labelcnt+1;
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




% THE FOLLOWING IS OLD
STOP_HERE





%plane 2 of bottles
v1 = Utri(:,85+26)-Utri(:,85+25);
v2 = Utri(:,85+27)-Utri(:,85+26);

n = cross(v1,v2);n= n/norm(n);
pl = [n;-n'*Utri(:,85+25)];

labeltype = [labeltype,1,1,2,2,1,1,1,2,2];
for ii=nbrpoints_anno+[7:15],
    Upp = U(1:3,ii);
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


%plane of saft soppa
v1 = Utri(:,32+26)-Utri(:,32+25);
v2 = Utri(:,32+27)-Utri(:,32+26);
v3 = Utri(:,32+27)-Utri(:,32+26);
v4 = Utri(:,32+14)-Utri(:,32+12);

n = cross(v1,v2);n= n/norm(n);
n2 = cross(v3,v4);n2= n2/norm(n2);

%n = n + 0.05*n2;n=n/norm(n);
pl = [n;-n'*Utri(:,32+25)-0.0001];

cc = 0;
for ii=nbrpoints_anno+[16:30],
    Upp = U(1:3,ii);
    labelcnt=labelcnt+1;
    cc = cc +1;
    if cc==5,
        cc=1;
    end
    if cc<=8,
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



%three packages in random position
pl = [n;-n'*Utri(:,32+25)];
labeltype = [labeltype,4,5,5];
for ii=nbrpoints_anno+[31:2:36],
    Upp1 = U(1:3,ii);
    Upp2 = U(1:3,ii+1);
    jj=(ii-nbrpoints_anno-31)/2;
    Upp = U(1:3,nbrpoints_anno+37+jj);
    
    if ii==nbrpoints_anno+31,
        nn = Upp2-Upp1;
        Upp = Upp+0.015*nn;
    end
    
    
    labelcnt=labelcnt+1;
        
    ll = -n'*Upp-pl(4);
    Uplane = Upp+ll*n;
    
    ll1 = -n'*Upp1-pl(4);
    Uplane1 = Upp1+ll1*n;
    
    ll2 = -n'*Upp2-pl(4);
    Uplane2 = Upp2+ll2*n;

    sc = norm(Uplane-Upp)/norm(-18.8-Ucc(3));

    R = [null(n'), n]; if det(R)<0,R=[R(:,2),R(:,1),R(:,3)];end
    n2ii=cross(Upp1-Upp2,n);n2ii=n2ii/norm(n2ii);
    tmp = R'*n2ii;
    tmp(3)=0;tmp=tmp/norm(tmp);
    th2 = atan2(tmp(1),tmp(2));
    
    if ii==nbrpoints_anno+31,
        th2 = th2 + 3*pi/180;
    end
    R = R*[cos(th2),sin(th2),0;-sin(th2),cos(th2),0;0,0,1];
    
    tt = Uplane - sc*R*[3.75;cos(th)*2.0;-18.8];
%    tt = Uplane1 - sc*R*[0,7.0,-18.8]';
    T = [sc*R,tt;[0,0,0,1]];

    Utmp=T(1:3,:)*pextend(Usaft);

    triobj = [triobj,trisaft+size(Uobj,2)];
    Uobj=[Uobj,Utmp];
    labelobj=[labelobj,labelcnt*ones(1,size(trisaft,2))];
end

%three packages in random position
pl = [n;-n'*Utri(:,32+1)-0.0001];
labeltype = [labeltype,4,5,5];
for ii=nbrpoints_anno+[40:3:45],
    Upp1 = U(1:3,ii);
    Upp2 = U(1:3,ii+1);
    Upp = U(1:3,ii+2);
    
    labelcnt=labelcnt+1;
        
    ll = -n'*Upp-pl(4);
    Uplane = Upp+ll*n;
    
    ll1 = -n'*Upp1-pl(4);
    Uplane1 = Upp1+ll1*n;
    
    ll2 = -n'*Upp2-pl(4);
    Uplane2 = Upp2+ll2*n;

    sc = norm(Uplane-Upp)/norm(-18.8-Ucc(3));

    R = [null(n'), n]; if det(R)<0,R=[R(:,2),R(:,1),R(:,3)];end
    n2ii=cross(Upp1-Upp2,n);n2ii=n2ii/norm(n2ii);
    tmp = R'*n2ii;
    tmp(3)=0;tmp=tmp/norm(tmp);
    th2 = atan2(tmp(1),tmp(2));
    
    R = R*[cos(th2),sin(th2),0;-sin(th2),cos(th2),0;0,0,1];
    
    tt = Uplane - sc*R*[3.75;cos(th)*2.0;-18.8];
    T = [sc*R,tt;[0,0,0,1]];

    Utmp=T(1:3,:)*pextend(Usaft);

    triobj = [triobj,trisaft+size(Uobj,2)];
    Uobj=[Uobj,Utmp];
    labelobj=[labelobj,labelcnt*ones(1,size(trisaft,2))];
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


Utop = U(1:3,nbrpoints_anno+[46:69]);
[uu,ss,vv]=svd(pextend(Utop)');
toppl = vv(:,4);
toppl = toppl/norm(toppl(1:3));

for ii=nbrpoints_anno+[46:69],
    Upp = U(1:3,ii);
    if ii==nbrpoints_anno+46 || ii==nbrpoints_anno+51,
        Upp = U(1:3,ii)+0.06*(U(1:3,ii+6)-U(1:3,ii));
    end
    if ii==nbrpoints_anno+47 || ii==nbrpoints_anno+49 || ii==nbrpoints_anno+50,
        Upp = U(1:3,ii)+0.03*(U(1:3,ii+6)-U(1:3,ii));
    end
    ll = -toppl(1:3)'*Upp-toppl(4);
    Upp = Upp+ll*toppl(1:3);
    
    
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

    



figure(1);clf;trisurf(tri',Utri(1,:),Utri(2,:),Utri(3,:));axis equal;rotate3d on;hold on;
trisurf(triobj',Uobj(1,:),Uobj(2,:),Uobj(3,:));


pause;
         
         

for imindex = 1:48,imindex
    annotate_plot(settings,P_uncalib,imindex, tri, Utri, triobj, Uobj, labelobj);
    pause
    close(imindex);
end


























