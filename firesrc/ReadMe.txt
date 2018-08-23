---------------------------------------------------
- Aline N, CIS563 Assignment 2, Fluid Sim
---------------------------------------------------

0. Building

This project requires Boost 1.37, Devil 1.6.7, Glut, and (optionally) renderman to build.  
Glut is included with the project.  The use of renderman can be disabled by removing the 
HAVE_PRMAN macro from the build settings.

I. Features

This assignment implements 3D smoke that lives in a closed box.  Features include:

 - Staggered grid (MACGrid creates several instances of GriddedData for each velocity direction)

 - Solids (Stored in mSolid in MACGrid.h)

 - Air, smoke density and temperature simulation (defined in MACGrid.h)
 -- Semi-lagrangian approach (See MACGrid::advectVelocity, MACGrid::advectTemperature, and MACGrid::advectDensity)
 -- Solve the poisson equation (See MACGridd::project)
 -- Add buoyancy forces (MacGrid::computeBouyancy)
 -- Add vorticity confinement forces (MACGrid::computeVorticityConfinement)

 - Rendering of smoke (MACGrid::draw)
 -- I render the smoke in sheets along either the X or Z directions, depending on the viewing direction
 -- I use a Cubic Spline (CubicGriddedData in GriddedData.h) for interpolating density and temperature

 - Use Conjugate Gradient algorithm to solve the system of linear equations (use this code ConjGrad.h (C++))
 -- Using a coordinate_matrix for A seemed to speed up the projection step from 1 fps to ~2/3 fps.
 -- I did try using the Jacobi preconditioner, but did not see a speedup.  Though cg_solve appears to need 
    a smaller number of iterations, the time savings were moot due to the overhead of using the preconditioner.

Sample output:

smoke_10x10x10.avi: Shows coarse smoke in a 10x10x10 grid
smoke_40x30x2.avi: Shows smoke in a think grid along with velocity, temperature, and pressure fields

II. UI Notes

Camera Controls: 
 - Alt-Left Button to rotate 
 - Alt-Middle Button to zoom
 - Reset camera with spacebar

Simulation Controls (available from right-click menu and keyboard):

 - Pause simulation ('=')
 - Play simulation  ('>')
 - Reset simulation ('<')
 - Toggle sources   ('s')
 - Toggle vorticity confinement
 - Toggle preconditioner
 - Hide solid
 - Make 40x30x2 Grid (Warning: slow)
 - Make 10x10x10 Grid (Warning: slow)

Display controls (available from right-click menu and keyboard):

 - Display smoke density ('1')
 - Display temperature   ('2')
 - Display pressure      ('3')
 - Toggle velocity field ('v')
 - Toggle vorticity forces ('f')
 - Toggle grid display     ('g')

III. Known Bugs/Limitations:

- Initializing the A and Jacobi preconditioner matrices using ublas's coordinate_matrix is very slow. 
  This effects startup and changing the grid size

- The smoke is not rendered correctly when viewed from the top or bottom because I did not write the
  code to support it.  Also, there are some sheet artifacts when viewing the smoke from the corners of
  the bounding box.

IV. References

- "Fluid Simulation: SIGGRAPH 2007 Course Notes" Robert Bridson, Matthias Muller-Fischer
- "Visual Simulation of Smoke", Ronald Fedkiw, Jos Stam, Henrik Wann Jensen
