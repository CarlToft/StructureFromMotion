% VISIONARY - A computer vision toolbox.
%
% 0. Help on datastructures & classes
%   datastructures.m - Explanation of notations, classes etc.
%
% 1. Intersection or structure from motion is the problem
% of calculating scene structure when the camera motion is known
%   intseclinear - points, lines, quadrics with linear method
%   intsecpoints - points, linear method
%   intsec2views - optimal twoview point method (Hartley/Sturm)
%   intseclines - lines, linear method
%   intsecconics - conics/quadrics, linear method (three or more views)
%
% 2. Resection or absolute orientation or motion from structure
% is the problem of calculating camera position orientation
% and internal calibration when something is known about
% the structure in the scene.
%  a) Resection: Calibrated.
%   resec3pts - resection using 3 points (4 solutions)
% 
%  b) Resection: Uncalibrated.
%   reseclinear - resection, linear method using points
%   reseccam3 - given two cameras, calculate third camera using points
% 
% 3. Structure and motion or relative orientation is the
% problem of calculating both scene structure and
% the relative positions of the cameras at the same time.
%  a) Assuming known interior calibration:
%
%  b) Assuming unknown and possibly varying interior calibration.
%   smshape   - shape based factorisation method for points
%   smaffine  - affine factorisation for points, lines and conics
%   smaffineclos - affine reconstruction with closure constraints
%   smeucl8p  - euclidean reconstruction with 8 point algorithm
%   smproj8p  - projective reconstruction with 8 point algorithm
%   focallength6pt - calculates the motion of 2 views with an unknown common
%       focallength using the 6 point algorithm.
%   smfocallengths  - semicalibrated reconstruction with two unknown focallengths
%   ransac2views - calculates the motion of 2 views using RANSAC
%   ransac3views - calculates the motion of 3 views using RANSAC
%   ransacfocallength - calculates the motion of 2 views with an
%       unknown common focallength using RANSAC and the 6pt mincase.
%   sm6points - calculates structure & motion for 6 points in 3 views
%   calibrated_fivepoint - 5 point 10 solutions minimal case solver for 2 views
%
% 4. Bundle adjustment consists of routines for optimising both
% structure and motion with respect to image measurements.
%  Bundle:
%   bundleplc - bundle adjustment for points, lines and conics
%               bundle adjustment for patches coming soon!!!
%
% 5. Selfcalibration or autocalibration is the
% problem of determining the interior calibration of
% the camera assuming some knowledge, e.g. that
% the intrinsic parameters are constant.
%   flexlinear - Flexible calibration a la Pollefeys
%   flexnonlinear - Flexible calibration optimizing the absolute quadric
%   bundleplc - auto & flexible calibration in bundle adjustment
%   
%
% 6. Multilinear geometry.
%   trifocallinear    - Linear estimate of the trifocal tensor
%   tritop	      - trifocal to motion object
%   ptotri	      - motion object to trifocal
%   bitop	      - bifocal(fundamental matrix) to motion object
%   ptobi	      - motion object to bifocal
%   invessentialm  - Essential matrix to t, R1, R2
%
% Matrix generation and manipulation.
%   m2v      	      - transform symmetric matrix to vector form.
%   v2m               - transform vector to symmetric matrix form.
%   rq                - rq-factorisation
%
% Image Processing routines.
%   harris	      - Harris corner detector
%   harrispro         - Harris corner detector + "semilocal" autocorrelation
%   findpatch         - find patch using affine correlation
%   manline           - extract a line from an image
%   manconic          - extract a conic from an image
%   manpolygon        - extract a polygon from an image (lines+points)
%
% Construct simulated data.
%   randomscene       - create random scene and images
%
% Quaternions and rotations.
%   rot2quat          - Rotation matrix to unit quaternion
%   quat2rot          - Unit quaternion to rotation matrix
%   randrot           - create a random rotation matrix.
%
% Visualize
%   visualize_export  - export a 3D reconstruction to our visualizer
%   visualize_import  - import a 3D reconstruction from our visualizer
%   test_visualize    - test script for exporting to visualize
%
% Snavely's bundler.
%   bundler_import    - import Noah Snavely files to visionary
%   test_bundler      - test script for importing bundler files
%   imdist_radial     - distort image points with 2-parameter radial model
%   imundist_radial   - undistort image points with 2-parameter radial model
%
% General routines.
%   pflat             - normalise homogeneous coordinates (Z=1)
%   psphere           - normalise homogeneous coordinates (sphere)
%   pextend           - extend to homogeneous coordinates (add ones)
%   skew              - return skew matrix from a vector
%   coniccamera       - calculate 6x10 conic projection matrix
%   findproj          - calculate projective transformation using points
%   findaff           - calculate affine transformation using points
%   findeucl          - calculate euclidean transformation using points
%   project           - project structure with motion to imagedata
%   reproject         - compare projected structure with imagedata
%   rmspoints         - calculate RMS for reprojected points
%   rmslines          - calculate RMS for reprojected lines
%   plotline          - plot 2D/3D lines
%   plotconic         - plot 2D conics
%   plotquadric       - plot 3D quadrics
%   fitline           - fit a line to point set with covariances
%   fitconic          - fit a conic to point set with convariances
%
% Local routines.
%   ginput2           - alternative to Matlab's ginput
%
% Demo, see directory 'demo/'.
%
% Mathematical Imaging Group at the department of mathematics,
% Lund University, Lund, Sweden.
%



