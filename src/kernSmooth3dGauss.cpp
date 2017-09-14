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
	// Smoothing 3-D data with a Gaussian kernel.
	//////////////////////// VARIABLES: ////////////////////////
	// 'psx' is the vector of x-coordinates.
	// 'psy' is the vector of y-coordinates.
	// 'psz' is the vector of z-coordinates.
	// 'Y' is the vector of responce variables.
	// 'L1' is the length of the first dimension.
	// 'L2' is the length of the second dimension.
	// 'L3' is the length of the third dimension.
	// 'hx' is the bandwidth in the x-direction.
	// 'hy' is the bandwidth in the y-direction.
	// 'hz' is the bandwidth in the z-direction.
	// 'wx' is the maximum extent of the kernel in the x-direction.
	// 'wy' is the maximum extent of the kernel in the y-direction.
	// 'wz' is the maximum extent of the kernel in the z-direction.
	// 'YSmooth' is the smoothed data to be returned.
	
	
	// A function that calculates the one-dimensional array index given the first (i) and the second (j) two-dimensional array index and the length of the first (Li) dimension of the array:
	int kernSmooth3dGauss_ind2d(int i, int j, int Li)
	{
		return i + j * Li;
	}
	
	// A function that calculates the one-dimensional array index given the first (i), the second (j), and third (k) three-dimensional array index and the length of the first (Li) and the second (Lj) dimension of the array:
	int kernSmooth3dGauss_ind3d(int i, int j, int k, int Li, int Lj)
	{
		return i + j * Li + k * Li * Lj;
	}
	
	// A function that calculates the one-dimensional array index given the first (i), the second (j), and third (k) three-dimensional array index and the length of the first (Li) and the second (Lj) dimension of the array:
	double kernSmooth3dGauss_GaussKern(double D1, double D2, double D3, double sigma1, double sigma2, double sigma3)
	{
		return 0.1591549 / (sigma1 * sigma2 * sigma3) * exp( -0.5 * ( D1*D1 / (sigma1*sigma1) + D2*D2 / (sigma2*sigma2) + D3*D3 / (sigma3*sigma3) ) );
	}
	
	// The main funciton for generating correlated vectors of autocorrelated pressure valuesf from sine waved:
	//void kernSmooth3dGauss(double *psx, double *psy, double *psz, double *Y, int *L1, int *L2, int *L3, double *hx, double *hy, double *hz, double *wx, double *wy, double *wz, double *YSmooth)
	//{
		
	void kernSmooth3dGauss(double psx[], double psy[], double psz[], double Y[], int *L1, int *L2, int *L3, double *hx, double *hy, double *hz, double *wx, double *wy, double *wz, double YSmooth[])
	{
		// 'kernsum' is the sum of the kernel at each point of smoothing:
		double kernsum = 0.0;
		// 'kern' is the kernel value at the current voxel:
		double kern = 0.0;
		// 'NAs' and 'NAs2' are the number of NAs counted along the first dimension, and along the first two dimensions, respectively:
		int NAs = 0;
		int NAs2 = 0;
		
		// 'N' is the number data points:
		int N = *L1 * *L2 * *L3;
		// 'distx', 'disty' and 'distz' are the distances from the reference voxel to the current voxel:
		double distx = 0.0;
		double disty = 0.0;
		double distz = 0.0;
		// 'pos1', 'pos2' and 'pos3' are the index positions in the array of a current compared to the reference voxel:
		int pos1 = 0;
		int pos2 = 0;
		int pos3 = 0;
		// 'ind_reference' and 'ind_current' are the indexes of the reference voxel and the current voxel:
		int ind_reference = 0;
		int ind_current = 0;
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
					// Get the index of the reference voxel:
					ind_reference = kernSmooth3dGauss_ind3d(i1,i2,i3,*L1,*L2);
					// Reinitialize the sum of the kernel:
					kernsum = 0.0;
					
					// Up3:
					pos3 = 0;
					up3 = true;
					while(up3 && i3+pos3 < *L3)
					{
						// Initialize the kernel and the NA count:
						kern = 0.0;
						NAs2 = 0;
						
						// Up3, Up2:
						pos2 = 0;
						up2 = true;
						while(up2 && i2+pos2 < *L2)
						{
							// Initialize the kernel:
							kern = 0.0;
							NAs = 0;
							
							// Up3, Up2, Up1:
							pos1 = 0;
							up1 = true;
							while(up1 && i1+pos1 < *L1)
							{
								// Get the distance to the current voxel:
								ind_current = kernSmooth3dGauss_ind3d(i1+pos1,i2+pos2,i3+pos3,*L1,*L2);
								distx = abs(psx[ind_reference] - psx[ind_current]);
								disty = abs(psy[ind_reference] - psy[ind_current]);
								distz = abs(psz[ind_reference] - psz[ind_current]);
								// Is the current voxel in range of the reference voxel?:
								if(distx < *wx && disty < *wy && distz < *wz)
								{
									// Exclude the NaNs:
									if(Y[ind_current] == Y[ind_current])
									{
										kern = kernSmooth3dGauss_GaussKern(distx, disty, distz, *hx, *hy, *hz);
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
							
							// Up3, Up2, Down1:
							pos1 = -1;
							down1 = true;
							while(down1 && i1+pos1>=0)
							{
								// Get the distance to the current voxel:
								ind_current = kernSmooth3dGauss_ind3d(i1+pos1,i2+pos2,i3+pos3,*L1,*L2);
								distx = abs(psx[ind_reference] - psx[ind_current]);
								disty = abs(psy[ind_reference] - psy[ind_current]);
								distz = abs(psz[ind_reference] - psz[ind_current]);
								// Is the current voxel in range of the reference voxel?:
								if(distx < *wx && disty < *wy && distz < *wz)
								{
									// Exclude the NaNs:
									if(Y[ind_current] == Y[ind_current])
									{
										kern = kernSmooth3dGauss_GaussKern(distx, disty, distz, *hx, *hy, *hz);
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
							NAs2 += NAs;
						}
						
						// Up3, Down2:
						pos2 = -1;
						down2 = true;
						while(down2 && i2+pos2 > 0)
						{
							// Initialize the kernel and the NA count:
							kern = 0.0;
							NAs = 0;
							
							// Up3, Down2, Up1:
							pos1 = 0;
							up1 = true;
							while(up1 && i1+pos1 < *L1)
							{
								// Get the distance to the current voxel:
								ind_current = kernSmooth3dGauss_ind3d(i1+pos1,i2+pos2,i3+pos3,*L1,*L2);
								distx = abs(psx[ind_reference] - psx[ind_current]);
								disty = abs(psy[ind_reference] - psy[ind_current]);
								distz = abs(psz[ind_reference] - psz[ind_current]);
								// Is the current voxel in range of the reference voxel?:
								if(distx < *wx && disty < *wy && distz < *wz)
								{
									// Exclude the NaNs:
									if(Y[ind_current] == Y[ind_current])
									{
										kern = kernSmooth3dGauss_GaussKern(distx, disty, distz, *hx, *hy, *hz);
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
							
							// Up3, Down2, Down1:
							pos1 = -1;
							down1 = true;
							while(down1 && i1+pos1>=0)
							{
								// Get the distance to the current voxel:
								ind_current = kernSmooth3dGauss_ind3d(i1+pos1,i2+pos2,i3+pos3,*L1,*L2);
								distx = abs(psx[ind_reference] - psx[ind_current]);
								disty = abs(psy[ind_reference] - psy[ind_current]);
								distz = abs(psz[ind_reference] - psz[ind_current]);
								// Is the current voxel in range of the reference voxel?:
								if(distx < *wx && disty < *wy && distz < *wz)
								{
									// Exclude the NaNs:
									if(Y[ind_current] == Y[ind_current])
									{
										kern = kernSmooth3dGauss_GaussKern(distx, disty, distz, *hx, *hy, *hz);
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
							NAs2 += NAs;
						}
						// If kern==0 and no NAs are registered along dimension 1 and 2, break the loop:
						if(kern == 0.0 && NAs2 == 0)
						{
							up3 = false;
						}
						pos3 += 1;
					}
					
					// Down3:
					pos3 = -1;
					down3 = true;
					while(down3 && i3+pos3 >= 0)
					{
						// Initialize the kernel and the NA count:
						kern = 0.0;
						NAs2 = 0;
						
						// Down3, Up2:
						pos2 = 0;
						up2 = true;
						while(up2 && i2+pos2 < *L2)
						{
							// Initialize the kernel:
							kern = 0.0;
							NAs = 0;
							
							// Down3, Up2, Up1:
							pos1 = 0;
							up1 = true;
							while(up1 && i1+pos1 < *L1)
							{
								// Get the distance to the current voxel:
								ind_current = kernSmooth3dGauss_ind3d(i1+pos1,i2+pos2,i3+pos3,*L1,*L2);
								distx = abs(psx[ind_reference] - psx[ind_current]);
								disty = abs(psy[ind_reference] - psy[ind_current]);
								distz = abs(psz[ind_reference] - psz[ind_current]);
								// Is the current voxel in range of the reference voxel?:
								if(distx < *wx && disty < *wy && distz < *wz)
								{
									// Exclude the NaNs:
									if(Y[ind_current] == Y[ind_current])
									{
										kern = kernSmooth3dGauss_GaussKern(distx, disty, distz, *hx, *hy, *hz);
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
							
							// Down3, Up2, Down1:
							pos1 = -1;
							down1 = true;
							while(down1 && i1+pos1>=0)
							{
								// Get the distance to the current voxel:
								ind_current = kernSmooth3dGauss_ind3d(i1+pos1,i2+pos2,i3+pos3,*L1,*L2);
								distx = abs(psx[ind_reference] - psx[ind_current]);
								disty = abs(psy[ind_reference] - psy[ind_current]);
								distz = abs(psz[ind_reference] - psz[ind_current]);
								// Is the current voxel in range of the reference voxel?:
								if(distx < *wx && disty < *wy && distz < *wz)
								{
									// Exclude the NaNs:
									if(Y[ind_current] == Y[ind_current])
									{
										kern = kernSmooth3dGauss_GaussKern(distx, disty, distz, *hx, *hy, *hz);
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
							NAs2 += NAs;
						}
						
						// Down3, Down2:
						pos2 = -1;
						down2 = true;
						while(down2 && i2+pos2 > 0)
						{
							// Initialize the kernel:
							kern = 0.0;
							NAs = 0;
							
							// Down3, Down2, Up1:
							pos1 = 0;
							up1 = true;
							while(up1 && i1+pos1 < *L1)
							{
								// Get the distance to the current voxel:
								ind_current = kernSmooth3dGauss_ind3d(i1+pos1,i2+pos2,i3+pos3,*L1,*L2);
								distx = abs(psx[ind_reference] - psx[ind_current]);
								disty = abs(psy[ind_reference] - psy[ind_current]);
								distz = abs(psz[ind_reference] - psz[ind_current]);
								// Is the current voxel in range of the reference voxel?:
								if(distx < *wx && disty < *wy && distz < *wz)
								{
									// Exclude the NaNs:
									if(Y[ind_current] == Y[ind_current])
									{
										kern = kernSmooth3dGauss_GaussKern(distx, disty, distz, *hx, *hy, *hz);
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
							
							// Down3, Down2, Down1:
							pos1 = -1;
							down1 = true;
							while(down1 && i1+pos1>=0)
							{
								// Get the distance to the current voxel:
								ind_current = kernSmooth3dGauss_ind3d(i1+pos1,i2+pos2,i3+pos3,*L1,*L2);
								distx = abs(psx[ind_reference] - psx[ind_current]);
								disty = abs(psy[ind_reference] - psy[ind_current]);
								distz = abs(psz[ind_reference] - psz[ind_current]);
								// Is the current voxel in range of the reference voxel?:
								if(distx < *wx && disty < *wy && distz < *wz)
								{
									// Exclude the NaNs:
									if(Y[ind_current] == Y[ind_current])
									{
										kern = kernSmooth3dGauss_GaussKern(distx, disty, distz, *hx, *hy, *hz);
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
							NAs2 += NAs;
						}
						// If kern==0 and no NAs are registered along dimension 1 and 2, break the loop:
						if(kern == 0.0)
						{
							down3 = false;
						}
						pos3 -= 1;
					} // End of while 3
					YSmooth[ind_reference] = YSmooth[ind_reference] / kernsum;	
				} // End of for i1
			} // End of for i2
		} // End of for i3
	} // End of void
} // End of extern "C"
