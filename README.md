# Geometric Modeling Final Project: 2D Harmonic Coordinates
In this project, 2D mesh deformation based on harmonic coordinates was implemented and visualized.

## Compile & Run
To compile:
```
mkdir build
cd build
cmake ..
make
```
To run the program, there are three options:
```
./final_bin
```
This will run the algorithm based on default mesh and cage. The mesh and cage will be simple square.
```
./final_bin <mesh.off>
```
This will read the mesh and automatically build a cage based on the boundary of the mesh. Note that the mesh coordinates must be in 3D, where the `z` coordinate is considered `0` (not used). Also, the mesh must be a triangle mesh, and must have a boundary.
```
./final_bin <mesh.off> <cage.cage>
```
This will read the mesh and cage from file. The cage file should only contain the cage coordinates, in 3D, one vertex per row, where the `z` coordinate is considered `0` (not used). These vertices must be either ordered clockwise or counterclockwise. Please take a look at `square.cage`, which is an example of this format.

## Keyboard & Mouse Control
Following interactions are supported in the UI:
- Use mouse to click and drag any cage vertex, to change its location. The new mesh will be calculated automatically.
- Use mouse to click any cage vertex, and then use `W/A/S/D` to change its location. The new mesh will be calculated automatically.
- Press `U` to undo the cage vertex change (only once).
- Press `R` to reset the cage and mesh to their original positions. 
