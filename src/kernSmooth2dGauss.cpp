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
	// Start: 2013-08-24 - Clean version.
	////////////////////// DESCRIPTION: ////////////////////////
	// Smoothing 2-D data with a Gaussian kernel.
	//////////////////////// VARIABLES: ////////////////////////
	// 'psx' is the vector of x-coordinates.
	// 'psy' is the vector of y-coordinates.
	// 'Y' is the vector of responce variables.
	// 'L1' is the length of the first dimension.
	// 'L2' is the length of the second dimension.
	// 'hx' is the bandwidth in the x-direction.
	// 'hy' is the bandwidth in the y-direction.
	// 'wx' is the maximum extent of the kernel in the x-direction.
	// 'wy' is the maximum extent of the kernel in the y-direction.
	// 'YSmooth' is the smoothed data to be returned.
	
	
	// A function that calculates the one-dimensional array index given the first (i) and the second (j) two-dimensional array index and the length of the first (Li) dimension of the array:
	int kernSmooth2dGauss_ind2d(int i, int j, int Li)
	{
		return i + j * Li;
	}
	
	// A function that returns the two-dimensional standard Gaussian distribution given distances and standard deviations D1 and sigma1 along the first dimension, and D2 and sigma2 along the second dimension. The value 1/sqrt((2*pi)^2) = 0.1591549431 is excluded due to division by the 'kernsum' in kernSmooth2dGauss():
	double kernSmooth2dGauss_GaussKern(double D1, double D2, double sigma1, double sigma2)
	{
		return exp( -0.5 * ( D1*D1 / (sigma1*sigma1) + D2*D2 / (sigma2*sigma2) ) ) / (sigma1 * sigma2);
	}
	
	// The main funciton for generating correlated vectors of autocorrelated pressure valuesf from sine waved:
	void kernSmooth2dGauss(double psx[], double psy[], double Y[], int *L1, int *L2, double *hx, double *hy, double *wx, double *wy, double YSmooth[])
	{
		// 'kernsum' is the sum of the kernel at each point of smoothing:
		double kernsum = 0.0;
		// 'kern' is the kernel value at the current voxel:
		double kern = 0.0;
		// 'NAs' is the number of NAs counted :
		int NAs = 0;
		
		// 'N' is the number data points:
		int N = *L1 * *L2;
		// 'distx' and 'disty' are the distances from the reference voxel to the current voxel:
		double distx = 0.0;
		double disty = 0.0;
		// 'pos1' and 'pos2' are the index positions in the array of a current compared to the reference voxel:
		int pos1 = 0;
		int pos2 = 0;
		// 'ind_reference' and 'ind_current' are the indexes of the reference voxel and the current voxel:
		int ind_reference = 0;
		int ind_current = 0;
		// Controls for the while loops:
		bool up1 = true;
		bool down1 = true;
		bool up2 = true;
		bool down2 = true;
		
		
		// Initiate the output:	
		for(int i = 0; i < N; i++)
		{
			YSmooth[i] = 0.0;
		}
		
		
		// Move through the voxels:
		for(int i2 = 0; i2 < *L2; i2++)
		{
			for(int i1 = 0; i1 < *L1; i1++)
			{
				// Get the index of the reference voxel:
				ind_reference = kernSmooth2dGauss_ind2d(i1,i2,*L1);
				// Reinitialize the sum of the kernel:
				kernsum = 0.0;
				
				// Up2:
				pos2 = 0;
				up2 = true;
				while(up2 && i2+pos2 < *L2) // while loop U2
				{
					// Initialize the kernel and the NA count:
					kern = 0.0;
					NAs = 0;
					
					// Up2, Up1:
					pos1 = 0;
					up1 = true;
					while(up1 && i1+pos1 < *L1) // while loop U2U1
					{
						// Get the distance to the current voxel:
						ind_current = kernSmooth2dGauss_ind2d(i1+pos1,i2+pos2,*L1);
						distx = abs(psx[ind_reference] - psx[ind_current]);
						disty = abs(psy[ind_reference] - psy[ind_current]);
						// Is the current voxel in range of the reference voxel?:
						if(distx < *wx && disty < *wy)
						{
							// Exclude the NaNs:
							if(Y[ind_current] == Y[ind_current])
							{
								kern = kernSmooth2dGauss_GaussKern(distx, disty, *hx, *hy);
								kernsum += kern;
								YSmooth[ind_reference] += Y[ind_current] * kern;
							}
							else
							{
								NAs += 1;
							}
						}
						else
						{
							up1 = false;
						}
						pos1 += 1;
					}
					
					// Up2, Down1:
					pos1 = -1;
					down1 = true;
					while(down1 && i1+pos1>=0) // while loop U2D1
					{
						// Get the distance to the current voxel:
						ind_current = kernSmooth2dGauss_ind2d(i1+pos1,i2+pos2,*L1);
						distx = abs(psx[ind_reference] - psx[ind_current]);
						disty = abs(psy[ind_reference] - psy[ind_current]);
						// Is the current voxel in range of the reference voxel?:
						if(distx < *wx && disty < *wy)
						{
							// Exclude the NaNs:
							if(Y[ind_current] == Y[ind_current])
							{
								kern = kernSmooth2dGauss_GaussKern(distx, disty, *hx, *hy);
								kernsum += kern;
								YSmooth[ind_reference] += Y[ind_current] * kern;
							}
							else
							{
								NAs += 1;
							}
						}
						else
						{
							down1 = false;
						}
						pos1 -= 1;
					}
					
					// If kern==0 and no NAs are registered along dimension 1, break the loop:
					if(kern == 0.0 && NAs == 0)
					{
						up2 = false;
					}
					pos2 += 1;
				}
				
				// Down2:
				pos2 = -1;
				down2 = true;
				while(down2 && i2+pos2 >= 0) // while loop D2
				{
					// Initialize the kernel and the NA count:
					kern = 0.0;
					NAs = 0;
					
					// Down2, Up1:
					pos1 = 0;
					up1 = true;
					while(up1 && i1+pos1 < *L1) // while loop D2U1
					{
						// Get the distance to the current voxel:
						ind_current = kernSmooth2dGauss_ind2d(i1+pos1,i2+pos2,*L1);
						distx = abs(psx[ind_reference] - psx[ind_current]);
						disty = abs(psy[ind_reference] - psy[ind_current]);
						// Is the current voxel in range of the reference voxel?:
						if(distx < *wx && disty < *wy)
						{
							// Exclude the NaNs:
							if(Y[ind_current] == Y[ind_current])
							{
								kern = kernSmooth2dGauss_GaussKern(distx, disty, *hx, *hy);
								kernsum += kern;
								YSmooth[ind_reference] += Y[ind_current] * kern;
							}
							else
							{
								NAs += 1;
							}
						}
						else
						{
							up1 = false;
						}
						pos1 += 1;
					}
					
					// Down2, Down1:
					pos1 = -1;
					down1 = true;
					while(down1 && i1+pos1>=0) // while loop D2D1
					{
						// Get the distance to the current voxel:
						ind_current = kernSmooth2dGauss_ind2d(i1+pos1,i2+pos2,*L1);
						distx = abs(psx[ind_reference] - psx[ind_current]);
						disty = abs(psy[ind_reference] - psy[ind_current]);
						// Is the current voxel in range of the reference voxel?:
						if(distx < *wx && disty < *wy)
						{
							// Exclude the NaNs:
							if(Y[ind_current] == Y[ind_current])
							{
								kern = kernSmooth2dGauss_GaussKern(distx, disty, *hx, *hy);
								kernsum += kern;
								YSmooth[ind_reference] += Y[ind_current] * kern;
							}
							else
							{
								NAs += 1;
							}
						}
						else
						{
							down1 = false;
						}
						pos1 -= 1;
					}
					// If kern==0 and no NAs are registered along dimension 1, break the loop:
					if(kern == 0.0 && NAs == 0)
					{
						down2 = false;
					}
					pos2 -= 1;
				}
			// Store the smoothed values:	
			YSmooth[ind_reference] = YSmooth[ind_reference] / kernsum;	
			} // End of for i1
		} // End of for i2
	} // End of void
} // End of extern "C"
