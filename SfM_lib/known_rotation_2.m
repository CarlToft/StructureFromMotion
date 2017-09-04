function imnames = known_rotation_2(settings)
save_path = settings.save_path;
visviews = settings.visviews;
pixtol2 = settings.pixtol2;
KK = settings.KK;
mininlnr = settings.mininlnr;

load(fullfile(save_path,'rotations3.mat'));
load(fullfile(save_path,'impoints4.mat'));

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

%[U,P,s] = krot4_feas_LP_sparse(u,A,pixtol2./KK(1,1),0.01,100);
[U,P,s] = krot4_feas_LP_sparse(u,A,pixtol2./KK(1,1),0.01,1);
s0 = s;
outlnr = sum(s > 1e-5);
if settings.debug_match
    fprintf('KnownRotationII: Nr of Outliers: %3d.\n',outlnr);
end

for i = 1:length(P);
    uu = NaN*ones(3,u.pointnr);
    uu(:,u.index{i}) = u.points{i};
    vis = find(isfinite(uu(1,:)));
    res = length(vis);
    outl = find(s(1:res) > 1e-5);
    uu(:,vis(outl)) = NaN;
    if settings.storesift==1,
        uusift = zeros(128,u.pointnr,'uint8');
        uusift(:,u.index{i}) = u.sift{i};
        u.sift{i} = uusift(:,isfinite(uu(1,:)));
    end
    u.index{i} = find(isfinite(uu(1,:)));
    u.points{i} = uu(:,u.index{i});
    s = s(res+1:end);
end

[U,P] = modbundle_sparse(U,P,u,20,0.01);

A = cell(size(P));
for i = 1:length(P);
    A{i} = P{i}(:,1:3);
end

%P_uncalib = cell(0);
%u_uncalib = u;
%for i = 1:length(P);
%    P_uncalib{i} = KK*P{i};
%    u_uncalib.points{i} = pflat(KK*u.points{i});
%    u_uncalib.index{i} = u.index;
%    u_uncalib.sift{i} = u.sift{i};
%end

%save(fullfile(save_path,'str_mot.mat'), 'U', 'P', 'P_uncalib', 'u', 'u_uncalib','imnames');
save(fullfile(save_path,'str_mot.mat'), 'U', 'P', 'u', 'imnames');
save(fullfile(save_path,'rotations3.mat'),'A');

function y = pflat(x)
y = x./repmat(x(end,:),[size(x,1) 1]);

