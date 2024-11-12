# mandelbrot.c

**A C program calculates the mandelbrot set and outputs it graphically in a window.**

The user can zoom in on the window by selecting an area with the mouse. The time it took to compute the pixels is also rendered in the window. The mandelbrot computations are performed in parallel by 1 - n slave processes and the graphical output and user interaction are handled by the master process 0.

The benefit of the parallellization is clearly noticeable. Here are some test execution results when executed on maximum.cs.abo.fi: 2 processes (1 master and 1 slave), the computation of the initial window took 0.16 seconds. 20 processes (1 master and 19 slaves), the computation of the initial window took 0.06 seconds.

Compile the program with `mpicc -mpe=graphics mandelbrot.c -o mandelbrot -lm`

Run the program with `mpiexec -n x ./mandelbrot` where `x` is the number of processes.

![Screenshot of the mandelbrot.c program running.](https://github.com/markarja/mandelbrot/blob/main/mandelbrot.png)

**Author**: Markus Karjalainen (uid markarja, matnr. 29849)
**Version**: 2010-11-06
