# Structure from Motion

## Installing Dependencies

### CVX
Download standard bundle from http://cvxr.com/cvx/download/

### MOSEK
1. Download here https://mosek.com/downloads. Academic license is available.
2. For Unix/Linux, start MATLAB with `LD_LIBRARY_PATH="<mosek_root>/8/tools/platform/linux64x86/bin" matlab`.
This is done in `matlab_start.sh`.

### VLFEAT
Install instructions: http://www.vlfeat.org/install-matlab.html

### Lowe´s SIFT implementation (optional)
Download and unpack http://www.cs.ubc.ca/~lowe/keypoints/siftDemoV4.zip.

### Calibration toolbox
Download and unpack http://www.vision.caltech.edu/bouguetj/calib_doc/.

## Initiate session
The setup code assumes that the above libs have been placed in a common folder , 
below refered to as the library path.

**For Unix/Linux:**
1. Start matlab with `./start_matlab.sh <your library path>`.
2. Run `setup.m`

**Otherwise:**
1. Run `setup.m <your library path>`

## Running the example
The example is executed with `reconstruct_corner.m` and uses data from the folder `corner`.
Settings for the reconstruction are found in `corner/reconstr_setup.m`.
1. `cd SfM_lib`
2. Run `reconstruct_scene.m`

## Running on your images
1. Put your images in a folder together with your settings file named `reconstr_setup.m`.
Copy the corner example settings file to use as a template.
2. Make sure the settings file points correctly to your calibration data.
3. Run `reconstruct_scene <path to data>`
