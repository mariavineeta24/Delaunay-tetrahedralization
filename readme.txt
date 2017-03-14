Delaunay Tetrahedralization :

There is a lot of computation involved for constructing DT. Hence, if the computation is done in initializeGL(), then the program would run even for 1000 points(when SDF is used for point generation). Otherwise the code will tend to crash.

For the case of user interface, the computation is taking place during run time, hence if the density of points is too much, then it is likely to crash. But, using SDF for point generation, the DT was constructed for density = 500 and 1000 points, but would take a while.

Voronoi cant be computed without having Delaunay computed as both are dual to each other.


