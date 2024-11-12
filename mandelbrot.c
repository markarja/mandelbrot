#include <stdio.h>
#include <stdlib.h>
#define MPE_GRAPHICS
#include <mpi.h>
#include "mpe.h"

/**
 * A C program calculates the mandelbrot set
 * and outputs it graphically in a window. The 
 * user can zoom in on the window by selecting
 * an area with the mouse. The time it took to 
 * compute the pixels is also rendered in the 
 * window. The mandelbrot computations are 
 * performed in parallel by 1 - n slave processes 
 * and the graphical output and user interaction 
 * are handled by the master process 0. 
 * 
 * The benefit of the parallellization is 
 * clearly noticeable. Here are some test 
 * execution results when executed on 
 * maximum.cs.abo.fi:
 * 
 * 2 processes (1 master and 1 slave), 
 * the computation of the initial window took 
 * 0.16 seconds.
 * 
 * 20 processes (1 master and 19 slaves), 
 * the computation of the initial window took 
 * 0.06 seconds.
 *
 * Compile the program with 'mpicc -mpe=graphics mandelbrot.c -o mandelbrot -lm'
 *
 * Run the program with 'mpiexec -n x ./mandelbrot' 
 * where x is the number of processes.  
 * 
 * @author Markus Karjalainen (uid markarja, matnr. 29849)
 * @version 2010-11-06
 */

typedef struct {
  double real;
  double imag;
} complex_type;

static MPE_XGraph window;
static char *displayname = "";

int main(int argc, char *argv[]) {
  const int tag = 1;
  int END_OF_WINDOW = -1;
  int STOP_SIGNAL = -2;
  int RESIZE_WINDOW = -3;
  int id, np;
  int W = 800;
  int H = 800;
  MPE_Color colors[256];
  int row[W + 1];
  int xy_plane[W][H];
  int selected_area[4];
  int rownum;
  MPI_Status status;
  int x = 0, y = 0; 
  int x1 = 0, y1 = 0;
  int x2 = 0, y2 = 0;

  double real_min = -2.0;
  double imag_min = -2.0;
  double real_max = 2.0;
  double imag_max = 2.0;
  double scale_real = (real_max - real_min)/(double)W;
  double scale_imag = (imag_max - imag_min)/(double)H;
  double start = 0;
  double end = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);

  if(np < 2) {
    if(id == 0)
      printf("This program can only run on a number of processes greater than 1.\n");
    MPI_Finalize();
    exit(0);
  }

  if(id == 0) {
   
    //Master process

    MPE_Open_graphics(&window, MPI_COMM_WORLD, displayname, 0, 0, W, H + 15, 0);  
    MPE_Make_color_array(window, 256, colors);

    while(1) {
      int i = 0;
      int s = 1; //the current slave
      int rows_sent = 0;
      int rows_received = 0;
      start = MPI_Wtime();
      
      //Send a row number to all slaves
      for(i;i < np - 1;i++) {
	MPI_Send(&i, 1, MPI_INT, i + 1, tag, MPI_COMM_WORLD);
        rows_sent++;
      }

      do {
        //Receive a row from a slave
	MPI_Recv(&row, W + 1, MPI_INT, s, tag, MPI_COMM_WORLD, &status);
        rows_received++;
        //Collect the results 
	for(x = 0;x < W;x++) {
          xy_plane[x][row[H]] = row[x];
	}
        //If there are rows left to compute, send next row to slave s
        if(rows_sent < H) {
	  MPI_Send(&i, 1, MPI_INT, s, tag, MPI_COMM_WORLD);
          rows_sent++;
          i++;
        }
	s = s + 1;
        if(s == np) s = 1;
      } while(rows_received != rows_sent);

      //Send a end of window signal to all slaves
      for(i = 1;i < np;i++) {
	MPI_Send(&END_OF_WINDOW, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
      }

      end = MPI_Wtime();    
      
      //Output the graphics
      for(x = 0;x < W;x++) {
	for(y = 0;y < H;y++) {
	  MPE_Draw_point(window, y, x, xy_plane[x][y]);
        }
      }
    
      //Output the execution time
      char buffer[30];
      sprintf(buffer, "The computing of the window took %2.2f seconds.", end - start);
      MPE_Fill_rectangle(window, 0, H, H, 15, MPE_WHITE);
      MPE_Draw_string(window, 5, H + 10, MPE_BLACK, buffer);   
    
      MPE_Get_drag_region(window, MPE_BUTTON1, MPE_DRAG_SQUARE, &x1, &y1, &x2, &y2);
      //If the x-axis of the  selected are is less than 10, stop 
      if(abs(x1 - x2) < 10) {
	for(i = 1;i < np;i++) {
	  MPI_Send(&STOP_SIGNAL, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
	}
        break;
      } else {
        //Otherwise, send the selected area to the slaves
        selected_area[0] = x1;
        selected_area[1] = y1;
        selected_area[2] = x2;
        selected_area[3] = y2;
        for(i = 1;i < np;i++) {
          MPI_Send(&RESIZE_WINDOW, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
	  MPI_Send(&selected_area, 4, MPI_INT, i, tag, MPI_COMM_WORLD);
	}
      }
    }

    MPE_Close_graphics(&window);
    MPI_Finalize();
    exit(0);

  } else {
    
    //Slave process
    
    while(1) {
      while(1) {
        //Receive a row number, compute the row 
	//and send it back to the master process 
	MPI_Recv(&rownum, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
	if(rownum < 0) break;
        complex_type c;
	c.imag = 0.0;
	c.real = 0.0;
        c.real = real_min + ((double)rownum*scale_real);
	for(x = 0;x < W;x++) {
	  c.imag = imag_min + ((double)x*scale_imag);
	  row[x] = mandelbrot(c);
        }
	row[W] = rownum;
        MPI_Send(&row, H + 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
      }
      //Received a stop signal, break the loop
      if(rownum == STOP_SIGNAL) break;
      
      //Received a resize window signal, receive new x1,y1,x2 and y2
      //coordinates and re-calculate the min and max values and the 
      //scaling factors.
      if(rownum == RESIZE_WINDOW) {
        MPI_Recv(&selected_area, 4, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
        x1 = selected_area[0];
        y1 = selected_area[1];
        x2 = selected_area[2];
        y2 = selected_area[3];
        double tmp_real_min = real_min;
        double tmp_imag_min = imag_min;
        double tmp_real_max = real_max;
        double tmp_imag_max = imag_max;
        real_min = tmp_real_min + (tmp_real_max - tmp_real_min)*(double)x1/(double)W;
        imag_min = tmp_imag_min + (tmp_imag_max - tmp_imag_min)*(double)y1/(double)H;
        real_max = tmp_real_min + (tmp_real_max - tmp_real_min)*(double)x2/(double)W;
        imag_max = tmp_imag_min + (tmp_imag_max - tmp_imag_min)*(double)y2/(double)H;
        scale_real = (real_max - real_min)/(double)W;
        scale_imag = (imag_max - imag_min)/(double)H;
      }
    }
    MPI_Finalize();
    exit(0);

  }
}

/**
 * Function that calculates the mandelbrot
 * value for a given complex number c.
 */
int mandelbrot(complex_type c) {
  int count = 0;
  int max = 255;
  complex_type z;
  double len2, temp;
  z.real = 0.0, z.imag = 0.0;
  
  do {
    temp = z.real * z.real - z.imag * z.imag + c.real;
    z.imag = 2.0 * z.real * z.imag + c.imag;
    z.real = temp;  
    len2 = z.real * z.real + z.imag * z.imag;
    if(len2 > 4.0) break;
    count++;
  } while (count < max);

  return count;
}
