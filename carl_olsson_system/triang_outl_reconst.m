function triang_outl_reconst(settings);
%Recompute 3D-points using 2 view triangulation-RANSAC for fixed cameras.

save_path = settings.save_path;
visviews = settings.visviews;
pixtol2 = settings.pixtol2;
KK = settings.KK;

load(strcat(save_path,'str_mot.mat'));
load(strcat(save_path,'impoints4.mat'));

vis = zeros(1,u.pointnr);
for i=1:length(u.points)
    imvis = zeros(1,u.pointnr);
    imvis(u.index{i}) = 1;
    vis = vis + imvis;
end
vis0 = vis;
vis = vis >= visviews;
for i = 1:length(u.points)
    pp = NaN*ones(3,u.pointnr);
    pp(:,u.index{i}) = u.points{i};
    pp = pp(:,vis);
    if settings.storesift==1,
        ppsift = zeros(128,u.pointnr,'uint8');
        ppsift(:,u.index{i}) = u.sift{i};
        ppsift = ppsift(:,vis);
        u.sift{i} = ppsift(:,isfinite(pp(1,:)));
    end
    u.index{i} = find(isfinite(pp(1,:)));
    u.points{i} = pp(:,isfinite(pp(1,:)));
end
u.pointnr = sum(vis);

%For each par seeing common points.
maxinliers = zeros(1,u.pointnr);
best_cam1 = NaN(1,u.pointnr);
best_cam2 = NaN(1,u.pointnr);
best_U = NaN(4,u.pointnr);

nrcams = zeros(1,u.pointnr);
for i = 1:length(P);
    nrcams(u.index{i}) = nrcams(u.index{i}) + 1;
end

for i = 1:100;
    i
    %randomly select two cameras
    camnr1 = ceil(rand(size(nrcams)).*nrcams);
    camnr2 = ceil(rand(size(nrcams)).*(nrcams-1));
    camnr2(camnr2 >= camnr1) = camnr2(camnr2 >= camnr1) + 1;
    
    cams1 = zeros(1,u.pointnr);
    cams2 = zeros(1,u.pointnr);
    cumnrcams = zeros(1,u.pointnr);
    vis = false(1,u.pointnr);
    
    for j = 1:length(P);
        cumnrcams(u.index{j}) = cumnrcams(u.index{j}) + 1; 
        vis(u.index{j}) = true;
        cams1(camnr1 == cumnrcams & vis) = j;        
        cams2(camnr2 == cumnrcams & vis) = j;        
        vis(u.index{j}) = false;
    end    
    %Triangulate
    U = intsec2views_midpoint_mult_cam(P,u,cams1,cams2);
        
    p = NaN(3,u.pointnr);
    inliers = zeros(1,u.pointnr);
    vis = false(1,u.pointnr);
    %Count inliers
    for k=1:length(P);
        p(:,u.index{k}) = u.points{k};
        vis(u.index{k}) = true;
        err = sqrt(sum((p(:,vis)-pflat(P{k}*U(:,vis))).^2));
        depth = P{k}(3,:)*U(:,vis);
        inliers(vis) = inliers(vis) + (err < pixtol2/KK(1,1) & depth > 0);
        p(:,u.index{k}) = NaN;
        vis(u.index{k}) = false;
    end
    best_cam1((inliers >= maxinliers) ) = cams1(inliers >= maxinliers);
    best_cam2((inliers >= maxinliers) ) = cams2(inliers >= maxinliers);
    best_U(:,(inliers >= maxinliers) ) = U(:,(inliers >= maxinliers) );
    maxinliers((inliers >= maxinliers) ) = inliers((inliers >= maxinliers) );
end

imseq = cell(1,2);
for i = 1:length(best_cam1);
    if mod(i,100) == 0
        i
    end
    P1 = P{best_cam1(i)};
    P2 = P{best_cam2(i)};
    
    p1 = NaN(3,u.pointnr);
    p1(:,u.index{best_cam1(i)}) = u.points{best_cam1(i)};
    
    p2 = NaN(3,u.pointnr);
    p2(:,u.index{best_cam2(i)}) = u.points{best_cam2(i)};
    
    best_U(:,i) = intsec2views(P1,P2,p1(:,i),p2(:,i));
end
p = NaN(3,u.pointnr);
maxinliers = zeros(1,u.pointnr);
vis = false(1,u.pointnr);
%Count inliers
for k=1:length(P);
    p(:,u.index{k}) = u.points{k};
    vis(u.index{k}) = true;
    err = sqrt(sum((p(:,vis)-pflat(P{k}*U(:,vis))).^2));
    depth = P{k}(3,:)*U(:,vis);
    maxinliers(vis) = maxinliers(vis) + (err < pixtol2/KK(1,1) & depth > 0);
    p(:,u.index{k}) = NaN;
    vis(u.index{k}) = false;
end

%remove points with less inliers than the threshold
%vis = maxinliers >= visviews;
vis = maxinliers >= 2;
for i = 1:length(u.points)
    pp = NaN*ones(3,u.pointnr);
    pp(:,u.index{i}) = u.points{i};
    pp = pp(:,vis);
    if settings.storesift==1,
        ppsift = zeros(128,u.pointnr,'uint8');
        ppsift(:,u.index{i}) = u.sift{i};
        ppsift = ppsift(:,vis);
        u.sift{i} = ppsift(:,isfinite(pp(1,:)));
    end
    u.index{i} = find(isfinite(pp(1,:)));
    u.points{i} = pp(:,isfinite(pp(1,:)));
end
u.pointnr = sum(vis);
best_U = best_U(:,vis);

%Create new point-structure
for i = 1:length(u.points)
    pp = NaN*ones(3,u.pointnr);
    pp(:,u.index{i}) = u.points{i};
    vis = find(isfinite(pp(1,:)));
    err = sqrt(sum((pp(:,vis)-pflat(P{i}*best_U(:,vis))).^2));
    depth = P{i}(3,:)*best_U(:,vis);
    %Test for inliers
    inliers = err < pixtol2/KK(1,1) & depth > 0;
    %Discard outliers
    pp(:,vis(~inliers)) = NaN;
    %Test outliers as new tracks?
    if settings.storesift==1,
        ppsift = zeros(128,u.pointnr,'uint8');
        ppsift(:,u.index{i}) = u.sift{i};
        u.sift{i} = ppsift(:,isfinite(pp(1,:)));
    end
    u.index{i} = find(isfinite(pp(1,:)));
    u.points{i} = pp(:,isfinite(pp(1,:)));    
end

%Bundle
[U,P] = modbundle_sparse(best_U,P,u,20,0.01);

%okalib kamror
%P_uncalib = P;
%u_uncalib = u;
%for i = 1:length(P);
%    P_uncalib{i} = KK*P{i};
%    u_uncalib.points{i} = pflat(KK*u.points{i});
%    u_uncalib.sift{i} = u.sift{i};
%end
save(strcat(save_path,'str_mot2.mat'), 'U', 'P', 'u', 'imnames');

function y = pflat(x)
y = x./repmat(x(end,:),[size(x,1) 1]);

function y = pextend(x)
y = [x; ones(1,size(x,2))];

