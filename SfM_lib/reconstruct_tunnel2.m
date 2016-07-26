%Read settings
curr=pwd;
cd C:\Users\fredrik\Documents\MATLAB\tunnels\experiment3_data;
%cd C:\Users\fredrik\Documents\MATLAB\tunnels\experiment2;
settings = reconstr_setup2;
cd(curr);

%Perform pairwise matching for all pairs, (unless camera_graph is
%specified)

pairwise_matching(settings);

%%%Track points throughout the image collection. 

create_imdata(settings); 

%%%Compute all the relative orientations. Uses the 5 point algorithm for essential matrix.

rel_or_5points(settings);

%Performs rotation averaging using a RANSAC type procedure
%Does not perform very well for larger datasets due to the dimension of the
%search space. Will implement a new one.

settings.imnames = RANSAC_orientations(settings);

%%%Remove cameras that can't be estimated.

settings.imnames = Remove_disjoint_cameras(settings);

%%%Re-track points

create_imdata_2(settings); 

%%%Solve known rotation problem and run bundle.

settings.imnames = known_rotation_1(settings);

%%%Solve known rotation problem with updated rotations and run bundle.

settings.imnames = known_rotation_2(settings);

%Re-triangulates points with current estmated cameras
%Might increase number of inliers slightly.

triang_outl_reconst(settings);

%Remove points that are uncertain in the depth direction.
[U,P,u] = remove_uncertin_points(settings);

T=[P{1};[0,0,0,1]];
for ii=1:length(P);
    P{ii}=P{ii}*inv(T);
end
U=T*U;


save(strcat(settings.save_path,'result_data.mat'),'u','U','P','settings');

%Plot result
plot_result(settings,P,U,u);
figure(12);clf;
plot(motion(P),'b');axis equal;

if 0,
    nbr=length(P);    
    for ii=1:nbr-1,
        tmp=length(intersect(u_uncalib.index{ii},u_uncalib.index{ii+1}));
        if tmp<25,
            tmp,ii,error('Less than 25');
        end
    end
    for ii=1:nbr-2,
        c1=pflat(null(P_uncalib{ii}));
        c2=pflat(null(P_uncalib{ii+1}));
        c3=pflat(null(P_uncalib{ii+2}));
        v1 = c2(1:3)-c1(1:3);
        v2 = c3(1:3)-c2(1:3);v2 = v2 /norm(v1);
        v1 = v1 / norm(v1);
        if norm(v1-v2)>0.2,
            v1,v2,ii,
            error('Large deviation in translation');
        end
    end
        
    
    
    
end

