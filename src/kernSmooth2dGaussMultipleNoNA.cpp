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
	// Start: 2013-08-23 - Clean version.
	////////////////////// DESCRIPTION: ////////////////////////
	// Smoothing 3-D data with a Gaussian kernel in the first two dimensions simultaneously for all arrays along the third dimension.
	//////////////////////// VARIABLES: ////////////////////////
	// 'psx' is the vector of x-coordinates.
	// 'psy' is the vector of y-coordinates.
	// 'Y' is the vector of responce variables.
	// 'L1' is the length of the first dimension.
	// 'L2' is the length of the second dimension.
	// 'L3' is the length of the third dimension.
	// 'hx' is the bandwidth in the x-direction.
	// 'hy' is the bandwidth in the y-direction.
	// 'wx' is the maximum extent of the kernel in the x-direction.
	// 'wy' is the maximum extent of the kernel in the y-direction.
	// 'YSmooth' is the smoothed data to be returned.
	
	
	// A function that calculates the one-dimensional array index given the first (i) and the second (j) two-dimensional array index and the length of the first (Li) dimension of the array:
	int kernSmooth2dGaussMultipleNoNA_ind2d(int i, int j, int Li)
	{
		return i + j * Li;
	}
	
	// A function that calculates the one-dimensional array index given the first (i), the second (j), and third (k) three-dimensional array index and the length of the first (Li) and the second (Lj) dimension of the array:
	int kernSmooth2dGaussMultipleNoNA_ind3d(int i, int j, int k, int Li, int Lj)
	{
		return i + j * Li + k * Li * Lj;
	}
	
	// A function that returns the two-dimensional standard Gaussian distribution given distances and standard deviations D1 and sigma1 along the first dimension, and D2 and sigma2 along the second dimension. The value 1/sqrt((2*pi)^2) = 0.1591549431 is excluded due to division by the 'kernsum' in kernSmooth2dGauss():
	double kernSmooth2dGaussMultipleNoNA_GaussKern(double D1, double D2, double sigma1, double sigma2)
	{
		return exp( -0.5 * ( D1*D1 / (sigma1*sigma1) + D2*D2 / (sigma2*sigma2) ) ) / (sigma1 * sigma2);
	}
	
	// The main funciton for generating correlated vectors of autocorrelated pressure valuesf from sine waved:
	void kernSmooth2dGaussMultipleNoNA(double psx[], double psy[], double Y[], int *L1, int *L2, int *L3, double *hx, double *hy, double *wx, double *wy, double YSmooth[])
	{
		// 'kernsum' is the sum of the kernel at each point of smoothing:
		double kernsum = 0.0;
		// 'kern' is the kernel value at the current voxel:
		double kern = 0.0;
		// 'N' is the number data points:
		int N = *L1 * *L2 * *L3;
		// 'distx', 'disty' and 'distz' are the distances from the reference voxel to the current voxel:
		double distx = 0.0;
		double disty = 0.0;
		// 'pos1', 'pos2' and 'pos3' are the index positions in the array of a current compared to the reference voxel:
		int pos1 = 0;
		int pos2 = 0;
		// 'ind_referencePos', 'ind_referenceY', 'ind_currentPos' and 'ind_currentY' are the indexes of the reference voxel in the data vector and the position vectors, and the indexes of the current voxel in the data vector and the position vectors:
		int ind_referencePos = 0;
		int ind_referenceY[*L3];
		for(int i3 = 0; i3 < *L3; i3++)
		{
			ind_referenceY[i3] = 0;
		}
		int ind_currentPos = 0;
		int ind_currentY[*L3];
		for(int i3 = 0; i3 < *L3; i3++)
		{
			ind_currentY[i3] = 0;
		}
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
				// Get the index of the reference voxel in the position vectors:
				ind_referencePos = kernSmooth2dGaussMultipleNoNA_ind2d(i1,i2,*L1);
				// Get the indexes of the reference voxel in the data vector:
				for(int i3 = 0; i3 < *L3; i3++)
				{
					ind_referenceY[i3] = kernSmooth2dGaussMultipleNoNA_ind3d(i1,i2,i3,*L1,*L2);
				}
				// Reinitialize the sum of the kernel:
				kernsum = 0.0;
				
				// Sum the data weighted by the kernel values over the volume where the kernel is defined > 0:
				
				// Upwards along dimension 2:
				pos2 = 0;
				up2 = true;
				while(up2 && i2+pos2 < *L2)
				{
					kern = 0.0;
					
					// Upwards along dimension 1 (Up2, Up1):
					pos1 = 0;
					up1 = true;
					while(up1 && i1+pos1 < *L1)
					{
						// Get the index of the current voxel in the position vectors:
						ind_currentPos = kernSmooth2dGaussMultipleNoNA_ind2d(i1+pos1,i2+pos2,*L1);
						// Get the indexes of the reference voxel in the data vector:
						for(int i3 = 0; i3 < *L3; i3++)
						{
							ind_currentY[i3] = kernSmooth2dGaussMultipleNoNA_ind3d(i1+pos1,i2+pos2,i3,*L1,*L2);
						}
						
						// Get the distances to the reference voxel:
						distx = abs(psx[ind_referencePos] - psx[ind_currentPos]);
						disty = abs(psy[ind_referencePos] - psy[ind_currentPos]);
						
						if(distx < *wx && disty < *wy)
						{
							kern = kernSmooth2dGaussMultipleNoNA_GaussKern(distx, disty, *hx, *hy);
							kernsum += kern;
							for(int i3 = 0; i3 < *L3; i3++)
							{
								YSmooth[ind_referenceY[i3]] += Y[ind_currentY[i3]] * kern;
							}
						}
						else
						{
							up1 = false;
						}
						pos1 += 1;
					} // End of Upwards along dimension 1:
					
					// Downwards along dimension 1 (Up2, Down1):
					pos1 = -1;
					down1 = true;
					while(down1 && i1+pos1>=0)
					{
						// Get the index of the current voxel in the position vectors:
						ind_currentPos = kernSmooth2dGaussMultipleNoNA_ind2d(i1+pos1,i2+pos2,*L1);
						// Get the indexes of the reference voxel in the data vector:
						for(int i3 = 0; i3 < *L3; i3++)
						{
							ind_currentY[i3] = kernSmooth2dGaussMultipleNoNA_ind3d(i1+pos1,i2+pos2,i3,*L1,*L2);
						}
						
						distx = abs(psx[ind_referencePos] - psx[ind_currentPos]);
						disty = abs(psy[ind_referencePos] - psy[ind_currentPos]);
						if(distx < *wx && disty < *wy)
						{
							kern = kernSmooth2dGaussMultipleNoNA_GaussKern(distx, disty, *hx, *hy);
							kernsum += kern;
							for(int i3 = 0; i3 < *L3; i3++)
							{
								YSmooth[ind_referenceY[i3]] += Y[ind_currentY[i3]] * kern;
							}
						}
						else
						{
							down1 = false;
						}
						pos1 -= 1;
					} // End of Downwards along dimension 1:
					
					if(kern == 0.0)
					{
						up2 = false;
					}
					pos2 += 1;
				} // End of Upwards along dimension 2:
				
				// Downwards along dimension 2:
				pos2 = -1;
				down2 = true;
				while(down2 && i2+pos2 > 0)
				{
					kern = 0.0;
					
					// Upwards along dimension 1 (Down2, Up1):
					pos1 = 0;
					up1 = true;
					while(up1 && i1+pos1 < *L1)
					{
						// Get the index of the current voxel in the position vectors:
						ind_currentPos = kernSmooth2dGaussMultipleNoNA_ind2d(i1+pos1,i2+pos2,*L1);
						// Get the indexes of the reference voxel in the data vector:
						for(int i3 = 0; i3 < *L3; i3++)
						{
							ind_currentY[i3] = kernSmooth2dGaussMultipleNoNA_ind3d(i1+pos1,i2+pos2,i3,*L1,*L2);
						}
						
						distx = abs(psx[ind_referencePos] - psx[ind_currentPos]);
						disty = abs(psy[ind_referencePos] - psy[ind_currentPos]);
						if(distx < *wx && disty < *wy)
						{
							kern = kernSmooth2dGaussMultipleNoNA_GaussKern(distx, disty, *hx, *hy);
							kernsum += kern;
							for(int i3 = 0; i3 < *L3; i3++)
							{
								YSmooth[ind_referenceY[i3]] += Y[ind_currentY[i3]] * kern;
							}
						}
						else
						{
							up1 = false;
						}
						pos1 += 1;
					} // End of Upwards along dimension 1:
					
					// Downwards along dimension 1 (Down2, Down1):
					pos1 = -1;
					down1 = true;
					while(down1 && i1+pos1>=0)
					{
						// Get the index of the current voxel in the position vectors:
						ind_currentPos = kernSmooth2dGaussMultipleNoNA_ind2d(i1+pos1,i2+pos2,*L1);
						// Get the indexes of the reference voxel in the data vector:
						for(int i3 = 0; i3 < *L3; i3++)
						{
							ind_currentY[i3] = kernSmooth2dGaussMultipleNoNA_ind3d(i1+pos1,i2+pos2,i3,*L1,*L2);
						}
						
						distx = abs(psx[ind_referencePos] - psx[ind_currentPos]);
						disty = abs(psy[ind_referencePos] - psy[ind_currentPos]);
						if(distx < *wx && disty < *wy)
						{
							kern = kernSmooth2dGaussMultipleNoNA_GaussKern(distx, disty, *hx, *hy);
							kernsum += kern;
							for(int i3 = 0; i3 < *L3; i3++)
							{
								YSmooth[ind_referenceY[i3]] += Y[ind_currentY[i3]] * kern;
							}
						}
						else
						{
							down1 = false;
						}
						pos1 -= 1;
					} // End of Downwards along dimension 1:
					
					if(kern == 0.0)
					{
						down2 = false;
					}
					pos2 -= 1;
				} // End of Downwards along dimension 2:
				
				for(int i3 = 0; i3 < *L3; i3++)
				{
					YSmooth[ind_referenceY[i3]] = YSmooth[ind_referenceY[i3]] / kernsum;
				} // End of for i3
			} // End of for i1
		} // End of for i2
	} // End of void
} // End of extern "C"
