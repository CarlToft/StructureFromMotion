function [U,P,u] = remove_uncertin_points(settings,U,P,u); 
save_path = settings.save_path;
uncertin_tol = settings.uncertin_tol;

if nargin<2,
    load(fullfile(save_path,'str_mot2.mat'));
end
d2f = zeros(size(U,2),9);

KK = settings.KK;
for i = 1:length(P)
    p = KK*P{i};
    r = p(:,1:3);
    t = p(:,4);
    uu = NaN*ones(3,u.pointnr);
    uu(:,u.index{i}) = KK*u.points{i};
    
    vis = isfinite(uu(1,:));
    
    a = uu(1,vis)'*r(3,:)-ones(size(uu(1,vis)'))*r(1,:);
    a0 = uu(1,vis)'*t(3,:)-t(1,:);
    
    b = uu(2,vis)'*r(3,:)-ones(size(uu(2,vis)'))*r(2,:);
    b0 = uu(2,vis)'*t(3,:)-t(2,:);
    
    c = ones(size(uu(1,vis)'))*r(3,:);
    c0 = ones(size(uu(1,vis)'))*t(3,:);
    
    x = sum(a.*U(1:3,vis)',2)+a0;
    y = sum(b.*U(1:3,vis)',2)+b0;
    z = sum(c.*U(1:3,vis)',2)+c0;
    
    %Compute the hessian
    h11 = 1./(z.^2);
    h12 = 0;
    h13 = -2*x./(z.^3);
    h21 = h12;
    h22 = 1./(z.^2);
    h23 = -2*y./(z.^3);
    h31 = h13;
    h32 = h23;
    h33 = 6*(x.^2+y.^2)./(z.^4);
    
    %columnstacked Hessian after change of coordinates.
    %term 1,1
    d2f(vis,1) = d2f(vis,1)     + a(:,1).*h11.*a(:,1)+a(:,1).*h12.*b(:,1) + a(:,1).*h13.*c(:,1) ...
        + b(:,1).*h21.*a(:,1)+b(:,1).*h22.*b(:,1) + b(:,1).*h23.*c(:,1) ...
        + c(:,1).*h31.*a(:,1)+c(:,1).*h32.*b(:,1) + c(:,1).*h33.*c(:,1);
    %term 2,1
    d2f(vis,2) = d2f(vis,2)     + a(:,2).*h11.*a(:,1)+a(:,2).*h12.*b(:,1) + a(:,2).*h13.*c(:,1) ...
        + b(:,2).*h21.*a(:,1)+b(:,2).*h22.*b(:,1) + b(:,2).*h23.*c(:,1) ...
        + c(:,2).*h31.*a(:,1)+c(:,2).*h32.*b(:,1) + c(:,2).*h33.*c(:,1);
    
    %term 3,1
    d2f(vis,3) = d2f(vis,3)     + a(:,3).*h11.*a(:,1)+a(:,3).*h12.*b(:,1) + a(:,3).*h13.*c(:,1) ...
        + b(:,3).*h21.*a(:,1)+b(:,3).*h22.*b(:,1) + b(:,3).*h23.*c(:,1) ...
        + c(:,3).*h31.*a(:,1)+c(:,3).*h32.*b(:,1) + c(:,3).*h33.*c(:,1);
    
    %term 1,2
    d2f(vis,4) = d2f(vis,4)     + a(:,1).*h11.*a(:,2)+a(:,1).*h12.*b(:,2) + a(:,1).*h13.*c(:,2) ...
        + b(:,1).*h21.*a(:,2)+b(:,1).*h22.*b(:,2) + b(:,1).*h23.*c(:,2) ...
        + c(:,1).*h31.*a(:,2)+c(:,1).*h32.*b(:,2) + c(:,1).*h33.*c(:,2);
    
    %term 2,2
    d2f(vis,5) = d2f(vis,5)     + a(:,2).*h11.*a(:,2)+a(:,2).*h12.*b(:,2) + a(:,2).*h13.*c(:,2) ...
        + b(:,2).*h21.*a(:,2)+b(:,2).*h22.*b(:,2) + b(:,2).*h23.*c(:,2) ...
        + c(:,2).*h31.*a(:,2)+c(:,2).*h32.*b(:,2) + c(:,2).*h33.*c(:,2);
    
    %term 3,2
    d2f(vis,6) = d2f(vis,6)     + a(:,3).*h11.*a(:,2)+a(:,3).*h12.*b(:,2) + a(:,3).*h13.*c(:,2) ...
        + b(:,3).*h21.*a(:,2)+b(:,3).*h22.*b(:,2) + b(:,3).*h23.*c(:,2) ...
        + c(:,3).*h31.*a(:,2)+c(:,3).*h32.*b(:,2) + c(:,3).*h33.*c(:,2);
    
    %term 1,3
    d2f(vis,7) = d2f(vis,7)     + a(:,1).*h11.*a(:,3)+a(:,1).*h12.*b(:,3) + a(:,1).*h13.*c(:,3) ...
        + b(:,1).*h21.*a(:,3)+b(:,1).*h22.*b(:,3) + b(:,1).*h23.*c(:,3) ...
        + c(:,1).*h31.*a(:,3)+c(:,1).*h32.*b(:,3) + c(:,1).*h33.*c(:,3);
    
    %term 2,3
    d2f(vis,8) = d2f(vis,8)     + a(:,2).*h11.*a(:,3)+a(:,2).*h12.*b(:,3) + a(:,2).*h13.*c(:,3) ...
        + b(:,2).*h21.*a(:,3)+b(:,2).*h22.*b(:,3) + b(:,2).*h23.*c(:,3) ...
        + c(:,2).*h31.*a(:,3)+c(:,2).*h32.*b(:,3) + c(:,2).*h33.*c(:,3);
    
    %term 3,3
    d2f(vis,9) = d2f(vis,9)     + a(:,3).*h11.*a(:,3)+a(:,3).*h12.*b(:,3) + a(:,3).*h13.*c(:,3) ...
        + b(:,3).*h21.*a(:,3)+b(:,3).*h22.*b(:,3) + b(:,3).*h23.*c(:,3) ...
        + c(:,3).*h31.*a(:,3)+c(:,3).*h32.*b(:,3) + c(:,3).*h33.*c(:,3);
    
    if 0
        vis = find(isfinite(uu(1,:)));
        point = vis(1);
        d = randn(3,1);
        fapprox = [];
        f0 = ((x(1)^2 + y(1)^2)/(z(1)^2));
        W = [a(1,:)' b(1,:)' c(1,:)'];
        gradf = W*2/z(1)^2*[x(1); y(1); (x(1)^2+y(1)^2)/(z(1)^2)];
        f = [];
        f2 = [];
        f3 = [];
        tint = -0.1:0.001:0.1;
        H = 1/(z(1)^2)*[1 0 -2*x(1)/z(1); 0 1 -2*y(1)/z(1); -2*x(1)/z(1) -2*y(1)/z(1) 6*(x(1)^2+y(1)^2)/(z(1)^2)];
        H = [a(1,:)' b(1,:)' c(1,:)']*H*[a(1,:); b(1,:); c(1,:)];
        for t = tint;
            res = 0;
            u2 = pflat(getpoints(imseq2{i}));
            f = [f sum((u2(:,point)- pflat(P{i}*(U(:,point)+t*[d;0]))).^2)];
            f2 = [f2 ((a(1,:)*(U(1:3,point)+t*d)+a0(1)).^2+(b(1,:)*(U(1:3,point)+t*d)+b0(1)).^2)./((c(1,:)*(U(1:3,point)+t*d)+c0(1)).^2)];
            f3 = [f3 f0+t^2*d'*H*d];
            fapprox = [fapprox f0+t*gradf'*d+t^2*d'*reshape(d2f(point,:),[3 3])*d];
        end
        plot(tint,[f2;f3;fapprox]);
    end
end
d2f = 2*d2f;

rmpoints = false(1,size(U,2));
for i = 1:length(P);
    
    uu = NaN*ones(3,u.pointnr);
    uu(:,u.index{i}) = KK*u.points{i};
    
    vis = isfinite(uu(1,:));
    PP = KK*P{i};
    c = -inv(PP(:,1:3))*PP(:,4);
    d = U(1:3,vis)-repmat(c,[1 sum(vis)]); %direction from the camera
    d = d./repmat(sqrt(d(1,:).^2+d(2,:).^2+d(3,:).^2),[3 1]);
    %rough change in function value d'*d2f*d
    deltaf =    d(1,:)'.*d2f(vis,1).*d(1,:)' + d(2,:)'.*d2f(vis,2).*d(1,:)' + d(3,:)'.*d2f(vis,3).*d(1,:)' +...
        d(1,:)'.*d2f(vis,4).*d(2,:)' + d(2,:)'.*d2f(vis,5).*d(2,:)' + d(3,:)'.*d2f(vis,6).*d(2,:)' +...
        d(1,:)'.*d2f(vis,7).*d(3,:)' + d(2,:)'.*d2f(vis,8).*d(3,:)' + d(3,:)'.*d2f(vis,9).*d(3,:)';
    
    vis = find(vis);
    rmpoints(vis(deltaf < uncertin_tol)) = true;
end
U = U(:,~rmpoints);

%for i = 1:size(u.points,2);
for i = 1:length(u.points);
    uu = NaN*ones(3,u.pointnr);
    uu(:,u.index{i}) = u.points{i};
    uu = uu(:,~rmpoints);
    if settings.storesift==1,
        usift = zeros(128,u.pointnr,'uint8');
        usift(:,u.index{i}) = u.sift{i};
        usift = usift(:,~rmpoints);
        u.sift{i} = usift(:,isfinite(uu(1,:)));
    end
    u.index{i} = find(isfinite(uu(1,:)));
    u.points{i} = uu(:,u.index{i});
end    
u.pointnr = sum(~rmpoints);
if nargin<2,
    save(fullfile(save_path,'red_str_mot.mat'),'u','U','P');
end
