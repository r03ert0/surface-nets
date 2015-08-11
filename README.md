# surface-nets
Isosurface extraction based on surface nets

Translated from javascript code by Mikola Lysenko, http://0fps.net/2012/07/12/smooth-voxel-terrain-part-2/

All the code is in one single file, without any dependencies. There are two functions, init_surfaceNets, and surfaceNets. The 'main' function shows how to use the code to transform a sphere in a volume into a mesh. The mesh format is:

* first line: number_of_vertices number_of_triangles
* x y z (for vertex 0)
* x y z (for vertex 1, etc)
* a b c (for triangle 0, indices are 0 based)
* a b c (for triangle 1, etc)

This is a render of the resulting mesh:

![sphere](https://cloud.githubusercontent.com/assets/2310732/9207777/f78b7822-406f-11e5-9a14-98044c06ad46.jpg)

There are some strange spikes on some cubes, but the result is overall ok.
