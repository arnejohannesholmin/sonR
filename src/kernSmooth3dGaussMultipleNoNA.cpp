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
	// Start: 2012-05-03 - Clean version.
	////////////////////// DESCRIPTION: ////////////////////////
	// Smoothing 3-D data with a Gaussian kernel. Note that four dimensional arrays are used as input, but the smoothing is done at each position in the fourth dimenstion separately. The fourth diemension is included in order to save CPU time.
	//////////////////////// VARIABLES: ////////////////////////
	// 'psx' is the vector of x-coordinates.
	// 'psy' is the vector of y-coordinates.
	// 'psz' is the vector of z-coordinates.
	// 'Y' is the vector of responce variables.
	// 'L1' is the length of the first dimension.
	// 'L2' is the length of the second dimension.
	// 'L3' is the length of the third dimension.
	// 'L4' is the length of the fourth dimension.
	// 'hx' is the bandwidth in the x-direction.
	// 'hy' is the bandwidth in the y-direction.
	// 'hz' is the bandwidth in the z-direction.
	// 'wx' is the maximum extent of the kernel in the x-direction.
	// 'wy' is the maximum extent of the kernel in the y-direction.
	// 'wz' is the maximum extent of the kernel in the z-direction.
	// 'YSmooth' is the smoothed data to be returned.
	
	
	// A function that calculates the one-dimensional array index given the first (i) and the second (j) two-dimensional array index and the length of the first (Li) dimension of the array:
	int kernSmooth3dGaussMultipleNoNA_ind2d(int i, int j, int Li)
	{
		return i + j * Li;
	}
	
	// A function that calculates the one-dimensional array index given the first (i), the second (j), and third (k) three-dimensional array index and the length of the first (Li) and the second (Lj) dimension of the array:
	int kernSmooth3dGaussMultipleNoNA_ind3d(int i, int j, int k, int Li, int Lj)
	{
		return i + j * Li + k * Li * Lj;
	}
	
	// A function that calculates the one-dimensional array index given the first (i), the second (j), the third (k), and the fourth (p) four-dimensional array index and the length of the first (Li), the second (Lj), and the third (Lk) dimension of the array:
	int kernSmooth3dGaussMultipleNoNA_ind4d(int i, int j, int k, int p, int Li, int Lj, int Lk)
	{
		return i + j * Li + k * Li * Lj + p * Li * Lj * Lk;
	}
	
	// A function that calculates the one-dimensional array index given the first (i), the second (j), and third (k) three-dimensional array index and the length of the first (Li) and the second (Lj) dimension of the array:
	double kernSmooth3dGaussMultipleNoNA_GaussKern(double D1, double D2, double D3, double sigma1, double sigma2, double sigma3)
	{
		return 0.1591549 / (sigma1 * sigma2 * sigma3) * exp( -0.5 * ( D1*D1 / (sigma1*sigma1) + D2*D2 / (sigma2*sigma2) + D3*D3 / (sigma3*sigma3) ) );
	}
	
	// The main funciton for generating correlated vectors of autocorrelated pressure valuesf from sine waved:
	//void kernSmooth3dGauss(double *psx, double *psy, double *psz, double *Y, int *L1, int *L2, int *L3, double *hx, double *hy, double *hz, double *wx, double *wy, double *wz, double *YSmooth)
	//{
		
	void kernSmooth3dGaussMultipleNoNA(double psx[], double psy[], double psz[], double Y[], int *L1, int *L2, int *L3,  int *L4, double *hx, double *hy, double *hz, double *wx, double *wy, double *wz, double YSmooth[])
	{
		// 'kernsum' is the sum of the kernel at each point of smoothing:
		double kernsum = 0.0;
		// 'kern' is the kernel value at the current voxel:
		double kern = 0.0;
		// 'N' is the number data points:
		int N = *L1 * *L2 * *L3 * *L4;
		// 'distx', 'disty' and 'distz' are the distances from the reference voxel to the current voxel:
		double distx = 0.0;
		double disty = 0.0;
		double distz = 0.0;
		// 'pos1', 'pos2' and 'pos3' are the index positions in the array of a current compared to the reference voxel:
		int pos1 = 0;
		int pos2 = 0;
		int pos3 = 0;
		// 'ind_referencePos', 'ind_referenceY', 'ind_currentPos' and 'ind_currentY' are the indexes of the reference voxel in the data vector and the position vectors, and the indexes of the current voxel in the data vector and the position vectors:
		int ind_referencePos = 0;
		std::vector<int>  ind_referenceY(*L4);
		
		for(int i4 = 0; i4 < *L4; i4++)
		{
			ind_referenceY[i4] = 0;
		}
		int ind_currentPos = 0;
		std::vector<int>  ind_currentY(*L4);
		for(int i4 = 0; i4 < *L4; i4++)
		{
			ind_currentY[i4] = 0;
		}
		// Controls for the while loops:
		bool up1 = true;
		bool down1 = true;
		bool up2 = true;
		bool down2 = true;
		bool up3 = true;
		bool down3 = true;
		
		// Initiate the output:	
		for(int i = 0; i < N; i++)
		{
			YSmooth[i] = 0.0;
		}
		
		
		
		// Move through the voxels:
		for(int i3 = 0; i3 < *L3; i3++)
		{
			for(int i2 = 0; i2 < *L2; i2++)
			{
				for(int i1 = 0; i1 < *L1; i1++)
				{
					// Get the index of the reference voxel in the position vectors:
					ind_referencePos = kernSmooth3dGaussMultipleNoNA_ind3d(i1,i2,i3,*L1,*L2);
					// Get the indexes of the reference voxel in the data vector:
					for(int i4 = 0; i4 < *L4; i4++)
					{
						ind_referenceY[i4] = kernSmooth3dGaussMultipleNoNA_ind4d(i1,i2,i3,i4,*L1,*L2,*L3);
					}
					// Reinitialize the sum of the kernel:
					kernsum = 0.0;
					
					// Sum the data weighted by the kernel values over the volume where the kernel is defined > 0:
					pos3 = 0;
					up3 = true;
					while(up3 && i3+pos3 < *L3)
					{
						kern = 0.0;
						
						pos2 = 0;
						up2 = true;
						while(up2 && i2+pos2 < *L2)
						{
							kern = 0.0;
							
							pos1 = 0;
							up1 = true;
							while(up1 && i1+pos1 < *L1)
							{
								// Get the index of the current voxel in the position vectors:
								ind_currentPos = kernSmooth3dGaussMultipleNoNA_ind3d(i1+pos1,i2+pos2,i3+pos3,*L1,*L2);
								// Get the indexes of the reference voxel in the data vector:
								for(int i4 = 0; i4 < *L4; i4++)
								{
									ind_currentY[i4] = kernSmooth3dGaussMultipleNoNA_ind4d(i1+pos1,i2+pos2,i3+pos3,i4,*L1,*L2,*L3);
								}
								
								// Get the distances to the reference voxel:
								distx = abs(psx[ind_referencePos] - psx[ind_currentPos]);
								disty = abs(psy[ind_referencePos] - psy[ind_currentPos]);
								distz = abs(psz[ind_referencePos] - psz[ind_currentPos]);
								
								if(distx < *wx && disty < *wy && distz < *wz)
								{
									kern = kernSmooth3dGaussMultipleNoNA_GaussKern(distx, disty, distz, *hx, *hy, *hz);
									kernsum += kern;
									for(int i4 = 0; i4 < *L4; i4++)
									{
										YSmooth[ind_referenceY[i4]] += Y[ind_currentY[i4]] * kern;
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
								// Get the index of the current voxel in the position vectors:
								ind_currentPos = kernSmooth3dGaussMultipleNoNA_ind3d(i1+pos1,i2+pos2,i3+pos3,*L1,*L2);
								// Get the indexes of the reference voxel in the data vector:
								for(int i4 = 0; i4 < *L4; i4++)
								{
									ind_currentY[i4] = kernSmooth3dGaussMultipleNoNA_ind4d(i1+pos1,i2+pos2,i3+pos3,i4,*L1,*L2,*L3);
								}
								
								distx = abs(psx[ind_referencePos] - psx[ind_currentPos]);
								disty = abs(psy[ind_referencePos] - psy[ind_currentPos]);
								distz = abs(psz[ind_referencePos] - psz[ind_currentPos]);
								if(distx < *wx && disty < *wy && distz < *wz)
								{
									kern = kernSmooth3dGaussMultipleNoNA_GaussKern(distx, disty, distz, *hx, *hy, *hz);
									kernsum += kern;
									for(int i4 = 0; i4 < *L4; i4++)
									{
										YSmooth[ind_referenceY[i4]] += Y[ind_currentY[i4]] * kern;
									}
								}
								else
								{
									down1 = false;
								}
								pos1 -= 1;
							}
							
							if(kern == 0.0)
							{
								up2 = false;
							}
							pos2 += 1;
						}
						
						pos2 = -1;
						down2 = true;
						while(down2 && i2+pos2 > 0)
						{
							kern = 0.0;
							
							pos1 = 0;
							up1 = true;
							while(up1 && i1+pos1 < *L1)
							{
								// Get the index of the current voxel in the position vectors:
								ind_currentPos = kernSmooth3dGaussMultipleNoNA_ind3d(i1+pos1,i2+pos2,i3+pos3,*L1,*L2);
								// Get the indexes of the reference voxel in the data vector:
								for(int i4 = 0; i4 < *L4; i4++)
								{
									ind_currentY[i4] = kernSmooth3dGaussMultipleNoNA_ind4d(i1+pos1,i2+pos2,i3+pos3,i4,*L1,*L2,*L3);
								}
								
								distx = abs(psx[ind_referencePos] - psx[ind_currentPos]);
								disty = abs(psy[ind_referencePos] - psy[ind_currentPos]);
								distz = abs(psz[ind_referencePos] - psz[ind_currentPos]);
								if(distx < *wx && disty < *wy && distz < *wz)
								{
									kern = kernSmooth3dGaussMultipleNoNA_GaussKern(distx, disty, distz, *hx, *hy, *hz);
									kernsum += kern;
									for(int i4 = 0; i4 < *L4; i4++)
									{
										YSmooth[ind_referenceY[i4]] += Y[ind_currentY[i4]] * kern;
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
								// Get the index of the current voxel in the position vectors:
								ind_currentPos = kernSmooth3dGaussMultipleNoNA_ind3d(i1+pos1,i2+pos2,i3+pos3,*L1,*L2);
								// Get the indexes of the reference voxel in the data vector:
								for(int i4 = 0; i4 < *L4; i4++)
								{
									ind_currentY[i4] = kernSmooth3dGaussMultipleNoNA_ind4d(i1+pos1,i2+pos2,i3+pos3,i4,*L1,*L2,*L3);
								}
								
								distx = abs(psx[ind_referencePos] - psx[ind_currentPos]);
								disty = abs(psy[ind_referencePos] - psy[ind_currentPos]);
								distz = abs(psz[ind_referencePos] - psz[ind_currentPos]);
								if(distx < *wx && disty < *wy && distz < *wz)
								{
									kern = kernSmooth3dGaussMultipleNoNA_GaussKern(distx, disty, distz, *hx, *hy, *hz);
									kernsum += kern;
									for(int i4 = 0; i4 < *L4; i4++)
									{
										YSmooth[ind_referenceY[i4]] += Y[ind_currentY[i4]] * kern;
									}
								}
								else
								{
									down1 = false;
								}
								pos1 -= 1;
							}
							
							if(kern == 0.0)
							{
								down2 = false;
							}
							pos2 -= 1;
						}
						
						if(kern == 0.0)
						{
							up3 = false;
						}
						pos3 += 1;
					}
					
					
					pos3 = -1;
					down3 = true;
					while(down3 && i3+pos3 >= 0)
					{
						kern = 0.0;
						
						pos2 = 0;
						up2 = true;
						while(up2 && i2+pos2 < *L2)
						{
							kern = 0.0;
							
							pos1 = 0;
							up1 = true;
							while(up1 && i1+pos1 < *L1)
							{
								// Get the index of the current voxel in the position vectors:
								ind_currentPos = kernSmooth3dGaussMultipleNoNA_ind3d(i1+pos1,i2+pos2,i3+pos3,*L1,*L2);
								// Get the indexes of the reference voxel in the data vector:
								for(int i4 = 0; i4 < *L4; i4++)
								{
									ind_currentY[i4] = kernSmooth3dGaussMultipleNoNA_ind4d(i1+pos1,i2+pos2,i3+pos3,i4,*L1,*L2,*L3);
								}
								
								distx = abs(psx[ind_referencePos] - psx[ind_currentPos]);
								disty = abs(psy[ind_referencePos] - psy[ind_currentPos]);
								distz = abs(psz[ind_referencePos] - psz[ind_currentPos]);
								if(distx < *wx && disty < *wy && distz < *wz)
								{
									kern = kernSmooth3dGaussMultipleNoNA_GaussKern(distx, disty, distz, *hx, *hy, *hz);
									kernsum += kern;
									for(int i4 = 0; i4 < *L4; i4++)
									{
										YSmooth[ind_referenceY[i4]] += Y[ind_currentY[i4]] * kern;
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
								// Get the index of the current voxel in the position vectors:
								ind_currentPos = kernSmooth3dGaussMultipleNoNA_ind3d(i1+pos1,i2+pos2,i3+pos3,*L1,*L2);
								// Get the indexes of the reference voxel in the data vector:
								for(int i4 = 0; i4 < *L4; i4++)
								{
									ind_currentY[i4] = kernSmooth3dGaussMultipleNoNA_ind4d(i1+pos1,i2+pos2,i3+pos3,i4,*L1,*L2,*L3);
								}
								
								distx = abs(psx[ind_referencePos] - psx[ind_currentPos]);
								disty = abs(psy[ind_referencePos] - psy[ind_currentPos]);
								distz = abs(psz[ind_referencePos] - psz[ind_currentPos]);
								if(distx < *wx && disty < *wy && distz < *wz)
								{
									kern = kernSmooth3dGaussMultipleNoNA_GaussKern(distx, disty, distz, *hx, *hy, *hz);
									kernsum += kern;
									for(int i4 = 0; i4 < *L4; i4++)
									{
										YSmooth[ind_referenceY[i4]] += Y[ind_currentY[i4]] * kern;
									}
								}
								else
								{
									down1 = false;
								}
								pos1 -= 1;
							}
							
							if(kern == 0.0)
							{
								up2 = false;
							}
							pos2 += 1;
						}
						
						pos2 = -1;
						down2 = true;
						while(down2 && i2+pos2 > 0)
						{
							kern = 0.0;
							
							pos1 = 0;
							up1 = true;
							while(up1 && i1+pos1 < *L1)
							{
								// Get the index of the current voxel in the position vectors:
								ind_currentPos = kernSmooth3dGaussMultipleNoNA_ind3d(i1+pos1,i2+pos2,i3+pos3,*L1,*L2);
								// Get the indexes of the reference voxel in the data vector:
								for(int i4 = 0; i4 < *L4; i4++)
								{
									ind_currentY[i4] = kernSmooth3dGaussMultipleNoNA_ind4d(i1+pos1,i2+pos2,i3+pos3,i4,*L1,*L2,*L3);
								}
								
								distx = abs(psx[ind_referencePos] - psx[ind_currentPos]);
								disty = abs(psy[ind_referencePos] - psy[ind_currentPos]);
								distz = abs(psz[ind_referencePos] - psz[ind_currentPos]);
								if(distx < *wx && disty < *wy && distz < *wz)
								{
									kern = kernSmooth3dGaussMultipleNoNA_GaussKern(distx, disty, distz, *hx, *hy, *hz);
									kernsum += kern;
									for(int i4 = 0; i4 < *L4; i4++)
									{
										YSmooth[ind_referenceY[i4]] += Y[ind_currentY[i4]] * kern;
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
								// Get the index of the current voxel in the position vectors:
								ind_currentPos = kernSmooth3dGaussMultipleNoNA_ind3d(i1+pos1,i2+pos2,i3+pos3,*L1,*L2);
								// Get the indexes of the reference voxel in the data vector:
								for(int i4 = 0; i4 < *L4; i4++)
								{
									ind_currentY[i4] = kernSmooth3dGaussMultipleNoNA_ind4d(i1+pos1,i2+pos2,i3+pos3,i4,*L1,*L2,*L3);
								}
								
								distx = abs(psx[ind_referencePos] - psx[ind_currentPos]);
								disty = abs(psy[ind_referencePos] - psy[ind_currentPos]);
								distz = abs(psz[ind_referencePos] - psz[ind_currentPos]);
								if(distx < *wx && disty < *wy && distz < *wz)
								{
									kern = kernSmooth3dGaussMultipleNoNA_GaussKern(distx, disty, distz, *hx, *hy, *hz);
									kernsum += kern;
									for(int i4 = 0; i4 < *L4; i4++)
									{
										YSmooth[ind_referenceY[i4]] += Y[ind_currentY[i4]] * kern;
									}
								}
								else
								{
									down1 = false;
								}
								pos1 -= 1;
							}
							
							if(kern == 0.0)
							{
								down2 = false;
							}
							pos2 -= 1;
						}
						
						if(kern == 0.0)
						{
							down3 = false;
						}
						pos3 -= 1;
					} // End of while 3
					for(int i4 = 0; i4 < *L4; i4++)
					{
						YSmooth[ind_referenceY[i4]] = YSmooth[ind_referenceY[i4]] / kernsum;
					} // End of for i4
				} // End of for i1
			} // End of for i2
		} // End of for i3
	} // End of void
} // End of extern "C"
