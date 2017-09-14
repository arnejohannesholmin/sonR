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
	// Smoothing 1-D data with a Gaussian kernel.
	//////////////////////// VARIABLES: ////////////////////////
	// 'psx' is the vector of x-coordinates.
	// 'Y' is the vector of responce variables.
	// 'L1' is the length of the first dimension.
	// 'hx' is the bandwidth in the x-direction.
	// 'wx' is the maximum extent of the kernel in the x-direction.
	// 'YSmooth' is the smoothed data to be returned.
	
	
	// A function that returns the two-dimensional standard Gaussian distribution given distances and standard deviations D1 and sigma1 along the first dimension, and D2 and sigma2 along the second dimension. The value 1/sqrt(2*pi) = 0.3989422804 is excluded due to division by the 'kernsum' in kernSmooth1dGauss():
	double kernSmooth1dGauss_GaussKern(double D1, double sigma1)
	{
		return exp( -0.5 * ( D1*D1 / (sigma1*sigma1) ) ) / sigma1;
	}
	
	// The main funciton for generating correlated vectors of autocorrelated pressure valuesf from sine waved:
	void kernSmooth1dGauss(double psx[], double Y[], int *L1, double *hx, double *wx, double YSmooth[])
	{
		// 'kernsum' is the sum of the kernel at each point of smoothing:
		double kernsum = 0.0;
		// 'kern' is the kernel value at the current voxel:
		double kern = 0.0;
		// 'N' is the number data points:
		int N = *L1;
		// 'distx' is the distances from the reference voxel to the current voxel:
		double distx = 0.0;
		// 'pos1' is the index positions in the array of a current compared to the reference voxel:
		int pos1 = 0;
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
			// Reinitialize the sum of the kernel:
			kernsum = 0.0;
			
			pos1 = 0;
			up1 = true;
			while(up1 && i1+pos1 < *L1)
			{
				distx = abs(psx[i1] - psx[i1+pos1]);
				
				if(distx < *wx)
				{
					// Exclude the NaNs:
					if(Y[i1+pos1] == Y[i1+pos1])
					{
						kern = kernSmooth1dGauss_GaussKern(distx, *hx);
						kernsum += kern;
						YSmooth[i1] += Y[i1+pos1] * kern;
					}
				}
				else
				{
					up1 = false;
				}
				pos1 += 1;
			}
			
			pos1 = -1;
			down1 = true;
			while(down1 && i1+pos1>=0)
			{
				distx = abs(psx[i1] - psx[i1+pos1]);
				if(distx < *wx)
				{
					// Exclude the NaNs:
					if(Y[i1+pos1] == Y[i1+pos1])
					{
						kern = kernSmooth1dGauss_GaussKern(distx, *hx);
						kernsum += kern;
						YSmooth[i1] += Y[i1+pos1] * kern;
					}
				}
				else
				{
					down1 = false;
				}
				pos1 -= 1;
			}
			
		YSmooth[i1] = YSmooth[i1] / kernsum;	
		} // End of for i1
	} // End of void
} // End of extern "C"
