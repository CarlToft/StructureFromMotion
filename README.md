# Structure from Motion

## Installing Dependencies

### CVX
Download standard bundle from http://cvxr.com/cvx/download/

### MOSEK
1. Download here https://mosek.com/resources/downloads. Academic license is available.
2. For Unix/Linux, start MATLAB with `LD_LIBRARY_PATH="<mosek_root>/8/tools/platform/linux64x86/bin" matlab`

### VLFEAT
Install instructions: http://www.vlfeat.org/install-matlab.html

### Lowe´s SIFT implementation
Download and unpack http://www.cs.ubc.ca/~lowe/keypoints/siftDemoV4.zip.

### Calibration toolbox
Download and unpack http://www.vision.caltech.edu/bouguetj/calib_doc/.

## Running the example
The example is executed with `reconstruct_scene.m` and uses data from the folder `corner`.
1. Go to repo root folder: `cd <your_path>/StructureFromMotion`
2. Edit the library path in `setup.m` to match yours.
3. Run `setup.m`
4. `cd SfM_lib`
5. Edit the VLFEAT, SIFT implementation paths in `corner/reconstr_setup_corner.m` to match yours.
6. Run `reconstruct_scene.m`
