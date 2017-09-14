#include <iostream>
using namespace std;
#include <valarray>
#include <fstream>
#include <vector>
#include <math.h>
#include <cmath>
#include <string>

extern "C" {
	
	//////////////////////// AUTHOR(S): ////////////////////////
	// Arne Johannes Holmin
	//////////////////////// LANGUAGE: /////////////////////////
	// English
	/////////////////////////// LOG: ///////////////////////////
	// Start: 2014-02-20 - Clean version.
	////////////////////// DESCRIPTION: ////////////////////////
	// Median smoothing 3-D data along the last dimension.
	//////////////////////// VARIABLES: ////////////////////////
	// 'Y' is the vector of responce variables.
	// 'L1' is the length of the first dimension.
	// 'L2' is the length of the second dimension.
	// 'L3' is the length of the third dimension.
	// 'YSmooth' is the smoothed data to be returned.
	
	
	// A function that calculates the one-dimensional array index given the first (i) and the second (j) two-dimensional array index and the length of the first (Li) dimension of the array:
	int rotate3d_ind2d(int i, int j, int Li)
	{
		return i + j * Li;
	}
	
	// A function that calculates the one-dimensional array index given the first (i), the second (j), and third (k) three-dimensional array index and the length of the first (Li) and the second (Lj) dimension of the array:
	int rotate3d_ind3d(int i, int j, int k, int Li, int Lj)
	{
		return i + j * Li + k * Li * Lj;
	}
	
	// The main funciton for generating correlated vectors of autocorrelated pressure valuesf from sine waved:
	void rotate3d(double x[],  int *Npoints, int by[], int *Lby, double ang[], double A[], int *Nrot, int *paired, double out[])
	{
		// Initialize the position indices used in the rotations:
		// The rotation matrices:
		int ind_A0 = 0;
		int ind_A1 = 0;
		int ind_A2 = 0;
		int ind_A3 = 0;
		int ind_A4 = 0;
		int ind_A5 = 0;
		int ind_A6 = 0;
		int ind_A7 = 0;
		int ind_A8 = 0;
		// The input 3D-points:
		int ind_x0 = 0;
		int ind_x1 = 0;
		int ind_x2 = 0;
		// The output 3D-points:
		int ind_out0 = 0;
		int ind_out1 = 0;
		int ind_out2 = 0;
		
		// Declare the angle in each step in the for loop:
		double thisang = 0.0;
	
		// Declare the sine and cosine values:
		double thissin = 0.0;
		double thiscos = 0.0;
	
		// Declare temporary values of 'A':
		double temp1 = 0.0;
		double temp2 = 0.0;
		double temp3 = 0.0;
		double temp4 = 0.0;
		double temp5 = 0.0;
		double temp6 = 0.0;
		
		int numElementsInOut = 0;
		if(*paired == 1)
		{
			numElementsInOut = *Npoints * 3;
		}
		else
		{
			numElementsInOut = *Npoints * *Nrot * 3;
		}
		
		// Initiate the output:	
		for(int i = 0; i < numElementsInOut; i++)
		{
			out[i] = 0.0;
		}
		// Initiate the output:	
		//for(int i = 0; i < 9 * *Nrot; i++)
		//{
		//	A[i] = 0.0;
		//}
		
		// 1. Generate the rotation matrices:
		// Move through the rows of angles:
		for(int i3 = 0; i3 < *Nrot; i3++)
		{
			// Get the current posision in A, and add 1, 2, ..., 8:
			ind_A0 = rotate3d_ind2d(0,i3,9);
			ind_A1 = ind_A0 + 1;
			ind_A2 = ind_A0 + 2;
			ind_A3 = ind_A0 + 3;
			ind_A4 = ind_A0 + 4;
			ind_A5 = ind_A0 + 5;
			ind_A6 = ind_A0 + 6;
			ind_A7 = ind_A0 + 7;
			ind_A8 = ind_A0 + 8;
			
			// Move through the rotations:
			for(int i2 = 0; i2 < *Lby; i2++)
			{
				// Get the current angle:
				thisang = ang[rotate3d_ind2d(i3,i2,*Nrot)];
				thissin = sin(thisang);
				thiscos = cos(thisang);
				// x-rotation:
				if(by[i2] == 0)
				{
					temp1 = A[ind_A1] *  thiscos + A[ind_A2] * thissin;
					temp2 = A[ind_A1] * -thissin + A[ind_A2] * thiscos;
					temp3 = A[ind_A4] *  thiscos + A[ind_A5] * thissin;
					temp4 = A[ind_A4] * -thissin + A[ind_A5] * thiscos;
					temp5 = A[ind_A7] *  thiscos + A[ind_A8] * thissin;
					temp6 = A[ind_A7] * -thissin + A[ind_A8] * thiscos;
					A[ind_A1] = temp1;
					A[ind_A2] = temp2;
					A[ind_A4] = temp3;
					A[ind_A5] = temp4;
					A[ind_A7] = temp5;
					A[ind_A8] = temp6;
				}
				// y-rotation:
				else if(by[i2] == 1)
				{
					temp1 = A[ind_A0] * thiscos + A[ind_A2] * -thissin;
					temp2 = A[ind_A0] * thissin + A[ind_A2] *  thiscos;
					temp3 = A[ind_A3] * thiscos + A[ind_A5] * -thissin;
					temp4 = A[ind_A3] * thissin + A[ind_A5] *  thiscos;
					temp5 = A[ind_A6] * thiscos + A[ind_A8] * -thissin;
					temp6 = A[ind_A6] * thissin + A[ind_A8] *  thiscos;
					A[ind_A0] = temp1;
					A[ind_A2] = temp2;
					A[ind_A3] = temp3;
					A[ind_A5] = temp4;
					A[ind_A6] = temp5;
					A[ind_A8] = temp6;
				}
				// z-rotation:
				else
				{
					temp1 = A[ind_A0] *  thiscos + A[ind_A1] * thissin;
					temp2 = A[ind_A0] * -thissin + A[ind_A1] * thiscos;
					temp3 = A[ind_A3] *  thiscos + A[ind_A4] * thissin;
					temp4 = A[ind_A3] * -thissin + A[ind_A4] * thiscos;
					temp5 = A[ind_A6] *  thiscos + A[ind_A7] * thissin;
					temp6 = A[ind_A6] * -thissin + A[ind_A7] * thiscos;
					A[ind_A0] = temp1;
					A[ind_A1] = temp2;
					A[ind_A3] = temp3;
					A[ind_A4] = temp4;
					A[ind_A6] = temp5;
					A[ind_A7] = temp6;
				} // End of z rotation
			} // End of for i2
		} // End of for i3
		// 2. If paired, rotate each 3D-point by the corresponding rotation matrix:
		if(*paired == 1)
		{
			for(int i3 = 0; i3 < *Npoints; i3++)
			{
				// Get the current posision in 'x':
				ind_x0 = rotate3d_ind2d(i3,0,*Npoints);
				ind_x1 = ind_x0 + *Npoints;
				ind_x2 = ind_x1 + *Npoints;
				// Get the current posision in 'A', and add 1, 2, ..., 8:
				ind_A0 = rotate3d_ind2d(0,i3,9);
				// ultiply 'A' with 'x':
				out[ind_x0] = x[ind_x0] * A[ind_A0] + x[ind_x1] * A[ind_A0+3] + x[ind_x2] * A[ind_A0+6];
				out[ind_x1] = x[ind_x0] * A[ind_A0+1] + x[ind_x1] * A[ind_A0+4] + x[ind_x2] * A[ind_A0+7];
				out[ind_x2] = x[ind_x0] * A[ind_A0+2] + x[ind_x1] * A[ind_A0+5] + x[ind_x2] * A[ind_A0+8];
			}	
		}
		// 3. Else, rotate each 3D-point by the each rotation matrix:
		else
		{
			for(int i3 = 0; i3 < *Npoints; i3++)
			{
				for(int i2 = 0; i2 < *Nrot; i2++)
				{
					// Get the current posision in 'x':
					ind_x0 = rotate3d_ind2d(i3,0,*Npoints);
					ind_x1 = ind_x0 + *Npoints;
					ind_x2 = ind_x1 + *Npoints;
					// Get the current posision in 'A':
					ind_A0 = rotate3d_ind2d(0,i2,9);
					// Get the current posision in 'out':
					ind_out0 = rotate3d_ind3d(i3, 0, i2, *Npoints, 3);
					ind_out1 = rotate3d_ind3d(i3, 1, i2, *Npoints, 3);
					ind_out2 = rotate3d_ind3d(i3, 2, i2, *Npoints, 3);
					// ultiply 'A' with 'x':
					out[ind_out0] = x[ind_x0] * A[ind_A0] + x[ind_x1] * A[ind_A0+3] + x[ind_x2] * A[ind_A0+6];
					out[ind_out1] = x[ind_x0] * A[ind_A0+1] + x[ind_x1] * A[ind_A0+4] + x[ind_x2] * A[ind_A0+7];
					out[ind_out2] = x[ind_x0] * A[ind_A0+2] + x[ind_x1] * A[ind_A0+5] + x[ind_x2] * A[ind_A0+8];
				} // End of for i2
			} // End of for i3
		} // End of paired
	} // End of void
} // End of extern "C"
// 
