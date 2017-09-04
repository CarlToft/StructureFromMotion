function imnames = known_rotation_1(settings)
save_path = settings.save_path;
visviews = settings.visviews;
pixtol1 = settings.pixtol1;
KK = settings.KK;
mininlnr = settings.mininlnr;

load(fullfile(save_path,'rotations2.mat'));
load(fullfile(save_path,'impoints3.mat'));

%keyboard;

%Ta bort punkter som inte syns i tillräckligt många kameror
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

%[U,P,s] = krot4_feas_LP_sparse(u,A,pixtol1./KK(1,1),0.1,100);
[U,P,s] = krot4_feas_LP_sparse(u,A,pixtol1./KK(1,1),0.01,1);
s0 = s;
outlnr = sum(s > 1e-5);
for i = 1:length(P);
    uu = NaN*ones(3,u.pointnr);
    uu(:,u.index{i}) = u.points{i};
    if settings.storesift==1,
        uusift = zeros(128,u.pointnr,'uint8');
        uusift(:,u.index{i}) = u.sift{i};
    end
    vis = find(isfinite(uu(1,:)));
    res = length(vis);
    outl = find(s(1:res) > 1e-5);
    uu(:,vis(outl)) = NaN;
    u.index{i} = isfinite(uu(1,:));
    u.points{i} = uu(:,u.index{i});
    if settings.storesift == 1,
        u.sift{i} = uusift(:,u.index{i});
    end
    s = s(res+1:end);
end

%Remove cameras with to few points
good_cams = true(size(P));
for i = 1:length(P)
    if sum(u.index{i},2) < mininlnr
        good_cams(i) = false;
    end
end
P = P(good_cams);
u.points = u.points(good_cams);
if settings.storesift==1,
    u.sift = u.sift(good_cams);
end
u.index = u.index(good_cams);
imnames = imnames(good_cams);

[U,P,lambda] = modbundle_sparse(U,P,u,20,0.01);

%Omskalning för visualize
T = diag([1 1 1 10]);

vis = 0;
for i=1:length(u.points)
    uu = NaN*ones(3,u.pointnr);
    uu(:,u.index{i}) = u.points{i};
    vis = vis + isfinite(uu(1,:));
end
vis = vis(1,:) >= 2;
for i = 1:length(u.points);
    uu = NaN*ones(3,u.pointnr);
    uu(:,u.index{i}) = u.points{i};
    uu = uu(:,vis);
    if settings.storesift==1,
        uusift = zeros(128,u.pointnr,'uint8');
        uusift(:,u.index{i}) = u.sift{i};
        uusift = uusift(:,vis);
        u.sift{i} = uusift(:,isfinite(uu(1,:)));
    end
    u.points{i} = uu(:,isfinite(uu(1,:)));
    u.index{i} = find(isfinite(uu(1,:)));
end
u.pointnr = sum(vis);

A = cell(size(P));
for i = 1:length(P);
    A{i} = P{i}(:,1:3);
end

%Spara structure rörelse
save(fullfile(save_path,'str_mot.mat'), 'U', 'P', 'u');
save(fullfile(save_path,'rotations3.mat'),'A');
load(fullfile(save_path,'impoints3.mat'));
u.points = u.points(good_cams);
if settings.storesift==1,
    u.sift = u.sift(good_cams);
end
u.index = u.index(good_cams);
imnames = imnames(good_cams);
save(fullfile(save_path,'impoints4.mat'),'u','imnames');



