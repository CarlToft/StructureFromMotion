
%OBJECT TYPES
%1 - cocacola pet 1.5 litre
%2 - cocacola light pet 1.5 litre
%3 - cocacola can 33 cl
%4 - blueberry soup
%5 - blackberry soup


%bib='C:\Users\fredrik\Documents\MATLAB\sfm_library\data\InstanceRecognition\instance_sequences\';
%bib='C:\Users\fredrik\Dropbox\tmp\Lucas\InstanceRecognition\instance_sequences\';
bib='/home/lucas/instance_deps/InstanceRecognition/instance_sequences/';

file='GOPR1830/';





%load settings and result file
% load([bib,file,'result.mat']);
% eval(['load ',bib,file,'result_anno.mat']);
load result_anno_updated_path


% Create data structure to determine the source of all points in U & u_uncalib
lookups = struct();
% Indices for all SFM points
lookups.sfm = create_lookups(U, u_uncalib);




% Add points for shelves and update reconstruction
% [u_uncalib,U]=add_1point(settings,P_uncalib,U,u_uncalib,[1,4,44,48]);
% [U,P_uncalib,lambda] = modbundle_sparse(U,P_uncalib,u_uncalib,10,10);
load result1_shelves_and_updated_cameras
lookups.shelves = struct('U_idx', {});
idx = size(U,2)+1 - (6+6+5);
lookups.shelves(1) = create_lookups(U, u_uncalib, 6, idx); idx = idx + 6;
lookups.shelves(2) = create_lookups(U, u_uncalib, 6, idx); idx = idx + 6;
lookups.shelves(3) = create_lookups(U, u_uncalib, 5, idx);




% Add points for upper right shelf of bottles.
% Bottles from right to left. 4 regular Coca-Cola, followed by 3 Coca-Cola Zero.
% [u_uncalib,U]=add_1point(settings,P_uncalib,U,u_uncalib,[60,68,44,48]);
load result2a_upper_right_bottles
lookups.right_shelf_bottles = struct('U_idx', {});
nbr_instances = 7;
for j = 1:nbr_instances
    idx = size(U,2)-nbr_instances + j;
    lookups.right_shelf_bottles(j) = create_lookups(U, u_uncalib, 1, idx);
end


% Add points for upper left shelf of bottles.
% Bottles from right to left. 6 regular Coca-Cola, followed by 2 Coca-Cola Zero.
% [u_uncalib,U]=add_1point(settings,P_uncalib,U,u_uncalib,[210,220,230,290]);
load result2b_upper_left_bottles
lookups.left_shelf_bottles = struct('U_idx', {});
nbr_instances = 8;
for j = 1:nbr_instances
    idx = size(U,2)-nbr_instances + j;
    lookups.left_shelf_bottles(j) = create_lookups(U, u_uncalib, 1, idx);
end


% Add points for saftsoppa on mid-upper shelf.
% [u_uncalib,U]=add_1point(settings,P_uncalib,U,u_uncalib,[210,220,230,290]);
load result3a_mid_upper_saftsoppa
lookups.midupper_shelf_saftsoppa = struct('U_idx', {});
nbr_instances = 9;
for j = 1:nbr_instances
    idx = size(U,2)-nbr_instances + j;
    lookups.midupper_shelf_saftsoppa(j) = create_lookups(U, u_uncalib, 1, idx);
end



mesh = struct();

mesh.shelves = struct('tri', {}, 'U', {});
%add bookshelf 1
[mesh.shelves(1).tri, mesh.shelves(1).U] = annotate_addshelf(U(1:3, lookups.shelves(1).U_idx), [1,22,26,27,25,21]);

%add bookshelf 2
% mesh.shelves(2) = struct();
[mesh.shelves(2).tri, mesh.shelves(2).U] = annotate_addshelf(U(1:3, lookups.shelves(2).U_idx), [22,2,1,26,25,7]);

%add bookshelf 3
% mesh.shelves(3) = struct();
[mesh.shelves(3).tri, mesh.shelves(3).U] = annotate_addshelf(U(1:3, lookups.shelves(3).U_idx), [2,6,1,25,22]);




% ADD OBJECTS
%labels
instance_cnt=0;


% Upper right shelf with bottles
mesh.right_shelf_bottles = struct('tri', {}, 'U', {});

% 4 Regular Coca-Cola bottles, 3 Coca-Cola Zero bottles
class_labels_tmp = [1,1,1,1,2,2,2];

plane = get_plane(mesh.shelves(1).U, 1);

for j=1:length(lookups.right_shelf_bottles),
    instance_cnt=instance_cnt+1;
    lid_position = U(1:3, lookups.right_shelf_bottles(j).U_idx);
    [mesh.right_shelf_bottles(j).tri, mesh.right_shelf_bottles(j).U] = annotate_addbottle(plane, lid_position);
    mesh.right_shelf_bottles(j).class_label = class_labels_tmp(j);
end




% Upper left shelf with bottles
mesh.left_shelf_bottles = struct('tri', {}, 'U', {});

% 6 Regular Coca-Cola bottles, 2 Coca-Cola Zero bottles
class_labels_tmp = [1,1,1,1,1,1,2,2];

plane = get_plane(mesh.shelves(3).U, 1);

for j=1:length(lookups.left_shelf_bottles),
    instance_cnt=instance_cnt+1;
    lid_position = U(1:3, lookups.left_shelf_bottles(j).U_idx);
    [mesh.left_shelf_bottles(j).tri, mesh.left_shelf_bottles(j).U] = annotate_addbottle(plane, lid_position);
    mesh.left_shelf_bottles(j).class_label = class_labels_tmp(j);
end



% 9 boxes of saftsoppa
[plane, front_direction] = get_plane(mesh.shelves(2).U, 1);
for j=1:length(lookups.midupper_shelf_saftsoppa),
    instance_cnt=instance_cnt+1;
    lid_position = U(1:3, lookups.midupper_shelf_saftsoppa(j).U_idx);
    [mesh.midupper_shelf_saftsoppa(j).tri, mesh.midupper_shelf_saftsoppa(j).U] = annotate_addsaftsoppa(plane, lid_position, front_direction);
    mesh.midupper_shelf_saftsoppa(j).class_label = 6;
end

% Merge mesh
[tri, Utri] = merge_mesh(mesh.shelves);

all_instances = [ ...
    mesh.right_shelf_bottles, ...
    mesh.left_shelf_bottles, ...
    mesh.midupper_shelf_saftsoppa];

[triobj, Uobj, instance_labels] = merge_mesh(all_instances);
class_labels = [all_instances.class_label];



% THE FOLLOWING IS OLD
STOP_HERE







alpha = 0:0.4:2*pi;rr = 0.24;height=1;nn=length(alpha);
Ubottle=[rr*cos(alpha),rr*cos(alpha),0;rr*sin(alpha),rr*sin(alpha),0;zeros(size(alpha)),height*ones(size(alpha)),1];
tribottle = [1:nn,      nn+1:2*nn, 2*nn+1*ones(1,nn) ; ...
             [2:nn,1],  [nn+2:2*nn,nn+1]  , nn+1:2*nn, ;...
             nn+1:2*nn, [2:nn,1]       , [nn+2:2*nn,nn+1]];
%figure(1);clf;trisurf(tribottle',Ubottle(1,:),Ubottle(2,:),Ubottle(3,:));axis equal;rotate3d on;



%plane of coke bottles
v1 = mesh.shelves(1).U(:,2)-mesh.shelves(1).U(:,1);
v2 = mesh.shelves(1).U(:,7)-mesh.shelves(1).U(:,2);

n = cross(v1,v2);n= n/norm(n);
pl = [n;-n'*mesh.shelves(1).U(:,1)];


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
    
    
    instance_cnt=instance_cnt+1;
    class_labels(instance_cnt) = 3; %cocacola can 33 cl
    ll = -n'*Upp-pl(4);
    Uplane = Upp+ll*n;

    sc = norm(Uplane-Upp);

    R = [null(n'), n]; if det(R)<0,R=[R(:,2),R(:,1),R(:,3)];end
    T = [sc*R,Uplane;[0,0,0,1]];

    Utmp=T(1:3,:)*pextend(Ubottle);

    triobj = [triobj,tribottle+size(Uobj,2)];
    Uobj=[Uobj,Utmp];
    instance_labels=[instance_labels,instance_cnt*ones(1,size(tribottle,2))];
end

    



figure(1);clf;trisurf(tri',mesh.shelves(1).U(1,:),mesh.shelves(1).U(2,:),mesh.shelves(1).U(3,:));axis equal;rotate3d on;hold on;
trisurf(triobj',Uobj(1,:),Uobj(2,:),Uobj(3,:));


pause;
         
         

for imindex = 1:48,imindex
    annotate_plot(settings,P_uncalib,imindex, tri, Utri, triobj, Uobj, instance_labels);
    pause
    close(imindex);
end
