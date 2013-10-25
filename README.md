particleChase - scalable particle chasing

We provide here the most basic but high scalable version of a particle chaser for a given velocity field.  Due to the underlying amazing fast <a href="https://github.com/cburstedde/p4est">p4est</a> forest of octree structure - which proved to work in peta-scale (e.g. on Jaguar Cray XT5) - this work can easily be adapted to implement more complex algorithms like the Fast Multipole Method.</p>


installing
----------------------
To test it on your computer, several steps are needed.

<code>
$ git clone https://github.com/plattenschieber/p4est "my fork which has only a little change in the vtk output 
$ cd p4est
$ git clone https://github.com/plattenschieber/particleChase
$ ./bootstrap "please make sure that which version of (g)libtoolize you use
$ ./configure --enable-mpi "--enable-debug
$ make
$ make install
$ mpirun -np 4 particleChase/particleChase 
</code>
