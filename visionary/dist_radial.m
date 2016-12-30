function imseq2 = dist_radial(imseq1,f,k)
%distort all points in an image sequence
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

rdist = norm(loc);
xc=(1+rdist^2*k(1)+rdist^4*k(2))*loc*f;

