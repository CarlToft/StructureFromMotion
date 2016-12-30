function imseq2 = undist_radial(imseq1,f,k)
%undistort all points in an image sequence
% INPUT:    imseq1 - cell of imagedata objects OR cell of 2xN matrices
%           f - vector of focal lengths
%           k - 2xn matrix of radial distortion parameter
% OUTPUT:   imseq2 - cell of imagedata objects
%
% MODEL: (x',y')=r(x,y)*(x,y) where r(x,y)=1+k(1)*(x^2+y^2)+k(2)*(x^2+y^2)^2

n=length(imseq1);
for ii=1:n,
    if isa(imseq1{ii},'imagedata'),
        tmp=pflat(getpoints(imseq1{ii}));tmp(3,:)=[];
    else
        tmp=imseq1{ii};
    end
    for jj=find(~isnan(tmp(1,:)));
        tmp(:,jj)=calibrate(tmp(:,jj),f(ii),k(:,ii));
    end
    imseq2{ii}=imagedata([],tmp);
end


function xc = calibrate(loc,f,k)

% loc is the point location, f is the focal length, k are the radial distortion parameters.

rdist = norm(loc/f);

pol = [k(2) 0 k(1) 0 1 -rdist];

rt = roots(pol);

rt = rt(abs(imag(rt))<eps);

rt = rt(rt > 0);

bb = abs(rt-rdist);

[m,g] = min(bb);

rt = rt(g);
if isempty(rt),
    scaling=1;
else
    scaling = rdist/rt;
end
xc = [loc/scaling/f];

