# Structure from Motion

## Installing Dependencies

### CVX
Download standard bundle from http://cvxr.com/cvx/download/

### MOSEK
1. Download here https://mosek.com/resources/downloads. Academic license is available.
2. For Unix/Linux, start MATLAB with `LD_LIBRARY_PATH="<mosek_root>/8/tools/platform/linux64x86/bin" matlab`.
This is done in `matlab_start.sh`.

### VLFEAT
Install instructions: http://www.vlfeat.org/install-matlab.html

### Lowe´s SIFT implementation (optional)
Download and unpack http://www.cs.ubc.ca/~lowe/keypoints/siftDemoV4.zip.

### Calibration toolbox
Download and unpack http://www.vision.caltech.edu/bouguetj/calib_doc/.

## Running the example
The setup code assumes that the above libs have been placed in a common folder , 
below refered to as the library path.

The example is executed with `reconstruct_scene.m` and uses data from the folder `corner`.
Settings for the reconstruction are found in `corner/reconstr_setup_corner.m`.

For Unix/Linux (and possibly iOS?):
1. Start matlab with `./start_matlab.sh <your library path>`.

Otherwise:
1a. Go to repo root folder: `cd <your_path>/StructureFromMotion`
1b. Edit the library path in `setup.m` to match yours.

Last steps are the same for all:
2. Run `setup.m`
3. `cd SfM_lib`
4. Run `reconstruct_scene.m`
