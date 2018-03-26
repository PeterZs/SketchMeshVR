# Project Title

Master thesis on a sketch-based 3D modeling system for VR.

## Getting Started
This git repository contains both a VR (VR-SketchMesh) and non-VR (SketchMesh) version of the sketch-based 3D modeling software. SketchMesh depends on the included version of libigl and VR-SketchMesh depends on the specially adapted version of libigl called VR-Viewer plus the OculusSDK. 

### Prerequisites

The VR version is tested on Windows 10 with Microsoft Visual Studio 15 2017 (with Visual C++), the non-VR version is additionally also tested on OS X El Capitan with clang 8.0.0.

### Installing

Clone the repository including submodules (this contains an adapted version of libigl that is necessary) with 

```
git clone --recursive https://github.com/FloorVerhoeven/thesis.git
```

Create a build folder in the "SketchMesh" and/or "VR-SketchMesh" folder.
Use CMake to build the makefiles for your system. You can either do this via the CMake GUI (press configure and then generate; use the "Unix Makefiles" generator on macOS or the "Visual Studio 15 2017 Win64" one on Windows) or via the command line by navigating into the build folder and running:

```
cmake -DCMAKE_BUILD_TYPE=Release ../
```

Then for macOS, use the terminal to navigate to inside the build folder run the following to compile and run the executable (replace ./SketchMesh_bin with ./VR-SketchMesh_bin for the VR version:
```
make && ./SketchMesh_bin
```
This will launch the application.

## Usage of VR-SketchMesh
At all times, the control curves are overlayed over the mesh.

* For drawing a initial mesh, press the GRIP and TRIGGER buttons and draw the shape (while pressing the buttons down). Slower drawing results in more internal triangles but also a slower processing.
* For adding additional control curves on the existing mesh surface, press the TRIGGER button and draw the stroke (while keeping the button pressed). For strokes that start outside of the mesh but end inside it, only the points that were drawn on the mesh will be considered. Strokes that start and end off the mesh are curently not supported
* For pulling (deforming) the shape, you first have to toggle away from drawing mode by pressing the B button. Then you have to press the GRIP and TRIGGER buttons and click close to one of the vertices of the stroke that you want to pull on, and while holding the buttons down, pull to the desired position. If you clicked too far away from any of the selectable vertices, nothing will happen.
* For removing additional control curves, first toggle away from curve adding mode by pressing the thumbstick buttonn. Afterwards click the TRIGGER button while pointing to one of the vertices of the stroke that you want to remove. Its vertices will turn black and you need to click the TRIGGER button one of the vertices of the same curve again to confirm the removal.
* For cutting the mesh, press the GRIP button and draw a cutting stroke that starts and ends outside of the mesh. After drawing the cutting stroke, click the GRIP button in order to perform the cut. The mesh part on the left-hand side of the drawn stroke will be removed (when looking in the direction of drawing of the stroke).
* For making extrusions, toggle away from cutting mode by pressing the A button. Next you press and hold the GRIP button to draw an extrusion base stroke on the mesh. This stroke needs to fully contain at least one mesh vertex. Then press and hold the GRIP button again and draw the silhouette stroke of the extrusion.


## Usage of non-VR SketchMesh
At all times, the control curves are overlayed over the mesh.

* For drawing a initial mesh, press 'D' and draw the shape (while pressing the mouse down). Slower drawing results in more internal triangles but also a slower processing.
* For adding additional control curves on the existing mesh surface, press 'A' and draw the stroke (while pressing the mouse down). For strokes that start outside of the mesh but end inside it, only the points that were drawn on the mesh will be considered. Strokes that start and end off the mesh are curently not supported
* For pulling (deforming) the shape, press 'P' and click close to one of the vertices of the stroke that you want to pull on, and while holding the mouse button down, pull to the desired position. The menu on the left hand side allows for switching between smooth and sharp deformation. If you click too far away from one of the stroke vertices, the mouse movement will be interpreted as a navigation movement.
* For removing additional control curves, press 'R' and click close to one of the vertices of the stroke that you want to remove. Its vertices will turn black and you need to click one of the vertices of the same curve again to confirm the removal.
* For cutting the mesh, press 'C' and draw a cutting stroke that starts and ends outside of the mesh. After drawing the cutting stroke, click in order to perform the cut. The mesh part on the left-hand side of the drawn stroke will be removed (when looking in the direction of drawing of the stroke).
* For making extrusions, press 'E' and draw an extrusion base stroke on the mesh. This stroke needs to fully contain at least one mesh vertex. Then press 'N' to turn the mesh to get a side view in which you can draw a silhouette stroke. Press 'E' again and draw the silhouette stroke outside of the mesh.
* For navigating, press 'N' and drag.
