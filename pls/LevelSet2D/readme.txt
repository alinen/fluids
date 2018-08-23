The file Container.h contains all the code neccessary to implement this example.
An object of class Container has fields that hold the initial grid, the level set, 
and the particles.  In the this example, 
the costructor of Container initializes the system using the MakeSphere function.
This function will initialize the values in the grid "init" to create an implicit surface 
representing Zalesak’s disk. A function like this is necessary to create the initial grid 
values representing the implicit surface. In the current version of this file, it is 
assumed that this function creates a signed distance function around the implicit function.
However, if this is not the case, the level set function "Reinitialize" has to be called
immediately after the "Initialize" function is called. The "Reinitialize" function
creates a signed distance function based on the interface location.

Once the system is initialized then calling Container::Update() from the appropriate function
(in our example glut_idle_cb()) will update the level set. In a more complex
example Update() will hold the simulator that computes the new velocities
of the grid points.
