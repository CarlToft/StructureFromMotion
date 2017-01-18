% VISIONARY/DATASTRUCTURES & CLASSES
%
% There are four main classes in Visionary.
%
% 1. imagedata - an image and its features, i.e. 2D data
% 2. imagesequence - a set of imagedata objects
% 3. structure - 3D data of scene
% 4. motion - camera motion rel. the scene and calibration
%
% The folling methods exist in these classes. For more information,
% type HELP CLASS/METHOD.
%
% 1. IMAGEDATA
%
%  imagedata - create an IMAGEDATA object
%  plot - plot imagedata
%  plus - add all features of an imagedata object to another
%  addpoints
%  addlines
%  addconics
%  getfilename - return filename of image
%  getpoints - return points in homogeneous coordinates
%  getlines - return lines in stored format
%  gethomogeneouslines - return lines in homogeneous coordinates
%  getconics - return conics in dual coordinates
%  clearpoints - clear points
%  clearlines
%  clearconics
%  getimage - return image matrix
%  loadimage - load image into memory with imread
%  clearimage - clear current image in memory
%  size - return number of points,lines & conics
%  changecsystem - change coordinate system
%  getnormtrans - return a "good" image transformation to better numerics
%
% 2. IMAGESEQUENCE
%  Not implemented
%  getcommonfeatures - return features present in all images
%
% 3. STRUCTURE
%
%  structure - create a STRUCTURE object
%  plot - plot structure in 3D
%  plus - add all features of a structure object to another
%  addpoints
%  addlines
%  addquadrics
%  getpoints - return points in homogeneous coordinates
%  getlines
%  getpluckerlines - plucker coordinates of lines
%  getquadrics
%  clearpoints - clear points
%  clearlines
%  clearquadrics
%  size - return number of points,lines & quadrics
%  changecsystem - change coordinate system
%
% 4. MOTION
%
%  motion - create a MOTION object
%  plot - plot motion orbit
%  addcameras - add cameras
%  getcameras - return selected cameras
%  size - return number of cameras
%  changecsystem - change coordinate system
%  focalpoints - calculate camera positions from camera matrices.
help datastructures

