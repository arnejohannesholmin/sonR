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
	// Smoothing 2-D data with a Gaussian kernel in the first dimension simultaneously for all arrays along the second dimension.
	//////////////////////// VARIABLES: ////////////////////////
	// 'psx' is the vector of x-coordinates.
	// 'Y' is the vector of responce variables.
	// 'L1' is the length of the first dimension.
	// 'L2' is the length of the second dimension.
	// 'hx' is the bandwidth in the x-direction.
	// 'wx' is the maximum extent of the kernel in the x-direction.
	// 'YSmooth' is the smoothed data to be returned.
	
	
	// A function that calculates the one-dimensional array index given the first (i) and the second (j) two-dimensional array index and the length of the first (Li) dimension of the array:
	int kernSmooth1dGaussMultipleNoNA_ind2d(int i, int j, int Li)
	{
		return i + j * Li;
	}
	
	// A function that returns the two-dimensional standard Gaussian distribution given distances and standard deviations D1 and sigma1 along the first dimension, and D2 and sigma2 along the second dimension. The value 1/sqrt(2*pi) = 0.3989422804 is excluded due to division by the 'kernsum' in kernSmooth1dGauss():
	double kernSmooth1dGaussMultipleNoNA_GaussKern(double D1, double sigma1)
	{
		return exp( -0.5 * ( D1*D1 / (sigma1*sigma1) ) ) / sigma1;
	}
	
	// The main funciton for generating correlated vectors of autocorrelated pressure valuesf from sine waved:
	void kernSmooth1dGaussMultipleNoNA(double psx[], double Y[], int *L1, int *L2, double *hx, double *wx, double YSmooth[])
	{
		// 'kernsum' is the sum of the kernel at each point of smoothing:
		double kernsum = 0.0;
		// 'kern' is the kernel value at the current voxel:
		double kern = 0.0;
		// 'N' is the number data points:
		int N = *L1 * *L2;
		// 'distx', 'disty' and 'distz' are the distances from the reference voxel to the current voxel:
		double distx = 0.0;
		// 'pos1', 'pos2' and 'pos3' are the index positions in the array of a current compared to the reference voxel:
		int pos1 = 0;
		// 'ind_referencePos', 'ind_referenceY', 'ind_currentPos' and 'ind_currentY' are the indexes of the reference voxel in the data vector and the position vectors, and the indexes of the current voxel in the data vector and the position vectors:
		int ind_referencePos = 0;
		int ind_referenceY[*L2];
		for(int i2 = 0; i2 < *L2; i2++)
		{
			ind_referenceY[i2] = 0;
		}
		int ind_currentPos = 0;
		int ind_currentY[*L2];
		for(int i2 = 0; i2 < *L2; i2++)
		{
			ind_currentY[i2] = 0;
		}
		// Controls for the while loops:
		bool up1 = true;
		bool down1 = true;
		
		// Initiate the output:	
		for(int i = 0; i < N; i++)
		{
			YSmooth[i] = 0.0;
		}
		
		
		
		// Move through the voxels:
		for(int i1 = 0; i1 < *L1; i1++)
		{
			// Get the index of the reference voxel in the position vectors:
			ind_referencePos = i1;
			// Get the indexes of the reference voxel in the data vector:
			for(int i2 = 0; i2 < *L2; i2++)
			{
				ind_referenceY[i2] = kernSmooth1dGaussMultipleNoNA_ind2d(i1,i2,*L1);
			}
			// Reinitialize the sum of the kernel:
			kernsum = 0.0;
			
			// Upwards along dimension 1:
			pos1 = 0;
			up1 = true;
			while(up1 && i1+pos1 < *L1)
			{
				// Get the index of the current voxel in the position vectors:
				ind_currentPos = i1+pos1;
				// Get the indexes of the reference voxel in the data vector:
				for(int i2 = 0; i2 < *L2; i2++)
				{
					ind_currentY[i2] = kernSmooth1dGaussMultipleNoNA_ind2d(i1+pos1,i2,*L1);
				}
				
				// Get the distances to the reference voxel:
				distx = abs(psx[ind_referencePos] - psx[ind_currentPos]);
				
				if(distx < *wx)
				{
					kern = kernSmooth1dGaussMultipleNoNA_GaussKern(distx, *hx);
					kernsum += kern;
					for(int i2 = 0; i2 < *L2; i2++)
					{
						YSmooth[ind_referenceY[i2]] += Y[ind_currentY[i2]] * kern;
					}
				}
				else
				{
					up1 = false;
				}
				pos1 += 1;
			} // End of Upwards along dimension 1:
			
			// Downwards along dimension 1:
			pos1 = -1;
			down1 = true;
			while(down1 && i1+pos1>=0)
			{
				// Get the index of the current voxel in the position vectors:
				ind_currentPos = i1+pos1;
				// Get the indexes of the reference voxel in the data vector:
				for(int i2 = 0; i2 < *L2; i2++)
				{
					ind_currentY[i2] = kernSmooth1dGaussMultipleNoNA_ind2d(i1+pos1,i2,*L1);
				}
				
				distx = abs(psx[ind_referencePos] - psx[ind_currentPos]);
				if(distx < *wx)
				{
					kern = kernSmooth1dGaussMultipleNoNA_GaussKern(distx, *hx);
					kernsum += kern;
					for(int i2 = 0; i2 < *L2; i2++)
					{
						YSmooth[ind_referenceY[i2]] += Y[ind_currentY[i2]] * kern;
					}
				}
				else
				{
					down1 = false;
				}
				pos1 -= 1;
			} // End of Downwards along dimension 1:
			
			for(int i2 = 0; i2 < *L2; i2++)
			{
				YSmooth[ind_referenceY[i2]] = YSmooth[ind_referenceY[i2]] / kernsum;
			} // End of for i2
		} // End of for i1
	} // End of void
} // End of extern "C"
