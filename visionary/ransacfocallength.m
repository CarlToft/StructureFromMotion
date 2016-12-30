function [bestmot, ind, bestd] = ransacfocallength(imseq, iteration, threshold)
% [mot,ind]=ransacfocallength(imdata,iteration,threshold)
% Calculate the motion of two views using ransac based on the 6pts in 2
% views minimal case for cameras with a common unknown focallength.
% Input:
%   imseq - cell of IMAGEDATA
%   iteration - iterations in RANSAC (optional)
%   threshold - for inliers in relative magnitude (optional)
% Output:
%   mot - motion
%   ind - indices of used points
%   bestd - reprojection errors
% Default: iteration = 100, threshold = 0.01 (relative);

if nargin < 2,
  iteration = 100;
end;
if nargin < 3
   threshold = 0.01;
end

xx1 = pflat(getpoints(imseq{1}));
xx2 = pflat(getpoints(imseq{2}));

x1 = xx1(1:2,:);
x2 = xx2(1:2,:);
scale = mean(abs([x1(:); x2(:)]));
threshold = scale * threshold;

n1 = size(xx1,2);
consistent = 0;

for i = 1 : iteration,
  % draw a random permutation and select the first 6 points  
  per=randperm(n1);
  x1 = xx1(:, per(1 : 6));
  x2 = xx2(:, per(1 : 6));
  
  % solve mincase
  mm = focallength6pt({imagedata([], x1), imagedata([], x2)});
  
  for j = 1 : length(mm)
      mot = mm{j};
      % triangulate & reproject
      str = intsec2views(mot, imseq);
      im1 = project(str, mot, 1);
      y1 = getpoints(im1);
      r1 = sum((y1 - xx1).^2);
      im2 = project(str, mot, 2);
      y2 = getpoints(im2);
      r2 = sum((y2 - xx2).^2);
      res = r1 + r2;

      % find inliers
      index = find(res < threshold^2);

      if length(index) > consistent,
          consistent = length(index);
          ind = index;
          bestd = sum(res);
          bestmot = mot;
      end
  end
end