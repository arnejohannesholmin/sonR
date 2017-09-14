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
	int median3d_ind2d(int i, int j, int Li)
	{
		return i + j * Li;
	}
	
	// A function that calculates the one-dimensional array index given the first (i), the second (j), and third (k) three-dimensional array index and the length of the first (Li) and the second (Lj) dimension of the array:
	int median3d_ind3d(int i, int j, int k, int Li, int Lj)
	{
		return i + j * Li + k * Li * Lj;
	}
	
	// A function that calculates the median of a sorted array, given the length and number of NAs:
	double median3d_getMedianOfSorted(double x[], int lengthOfArray, int NumberOfNAs)
	{
		double out = 0.0;
		int NumberOfNotNAs = lengthOfArray - NumberOfNAs;
		
		if(lengthOfArray == NumberOfNAs)
		{
			out = 0.0 / 0.0;
		}
		else if(NumberOfNotNAs % 2 == 0)
		{
			NumberOfNotNAs = NumberOfNotNAs / 2;
			out = (x[NumberOfNotNAs-1] + x[NumberOfNotNAs]) / 2;
		}
		else
		{
			NumberOfNotNAs = (NumberOfNotNAs - 1) / 2;
			out = x[NumberOfNotNAs];	
		}
		return out;
	}
	
	// The main funciton for generating correlated vectors of autocorrelated pressure valuesf from sine waved:
	void median3d(double Y[],  int *L1, int *L2, int *L3, int *margin, double YSmooth[], double sorted[])
	{
		// Initialize the position indices used in the smoothing:
		int ind_smooth = 0;
		int ind_current = 0;
		// The number of NAs:
		int NAs = 0;
		
		// 1. If smoothing along the first dimension:
		if(*margin==1)
		{
			// Initiate the 'sorted':	
			for(int i = 0; i < *L1; i++)
			{
				sorted[i] = 0.0;
			}
			// Initiate the output:	
			for(int i = 0; i < *L2 * *L3; i++)
			{
				YSmooth[i] = 0.0;
			}
			
			// Move through the voxels:
			for(int i3 = 0; i3 < *L3; i3++)
			{
				for(int i2 = 0; i2 < *L2; i2++)
				{
					// Reset the number of NAs:
					NAs = 0;
					// Get the values to sort:
					ind_smooth = median3d_ind2d(i2,i3,*L2);
					for(int i1 = 0; i1 < *L1; i1++)
					{
						ind_current = median3d_ind3d(i1,i2,i3,*L1,*L2);
						if(Y[ind_current] != Y[ind_current])
						{
							NAs += 1;
							//sorted[i1] = 1.0 / 0;
							sorted[i1] = 999999999999999999;
						}
						else
						{
							sorted[i1] = Y[ind_current];
						}
					}
					// Sort:
					std::sort(sorted, sorted + *L1);
					// Get the median:
					YSmooth[ind_smooth] = median3d_getMedianOfSorted(sorted, *L1, NAs);
				} // End of for i2
			} // End of for i3
		} // End of *margin=1
		
		
		
		// 2. If smoothing along the second dimension:
		else if(*margin==2)
		{
			// Initiate the 'sorted':	
			for(int i = 0; i < *L2; i++)
			{
				sorted[i] = 0.0;
			}
			// Initiate the output:	
			for(int i = 0; i < *L1 * *L3; i++)
			{
				YSmooth[i] = 0.0;
			}
			
			// Move through the voxels:
			for(int i3 = 0; i3 < *L3; i3++)
			{
				for(int i1 = 0; i1 < *L1; i1++)
				{
					// Reset the number of NAs:
					NAs = 0;
					// Get the values to sort:
					ind_smooth = median3d_ind2d(i1,i3,*L1);
					for(int i2 = 0; i2 < *L2; i2++)
					{
						ind_current = median3d_ind3d(i1,i2,i3,*L1,*L2);
						if(Y[ind_current] != Y[ind_current])
						{
							NAs += 1;
							//sorted[i2] = 1.0 / 0;
							sorted[i2] = 999999999999999999;
						}
						else
						{
							sorted[i2] = Y[ind_current];
						}
					}
					// Sort:
					std::sort(sorted, sorted + *L2);
					// Get the median:
					YSmooth[ind_smooth] = median3d_getMedianOfSorted(sorted, *L2, NAs);
				} // End of for i1
			} // End of for i3
		} // End of *margin=2
		
		
		
		// 3. If smoothing along the third dimension:
		else
		{
			// Initiate 'sorted':	
			for(int i = 0; i < *L3; i++)
			{
				sorted[i] = 0.0;
			}
			// Initiate the output:	
			for(int i = 0; i < *L1 * *L2; i++)
			{
				YSmooth[i] = 0.0;
			}
			
			// Move through the voxels:
			for(int i2 = 0; i2 < *L2; i2++)
			{
				for(int i1 = 0; i1 < *L1; i1++)
				{
					// Reset the number of NAs:
					NAs = 0;
					// Get the values to sort:
					ind_smooth = median3d_ind2d(i1,i2,*L1);
					for(int i3 = 0; i3 < *L3; i3++)
					{
						ind_current = median3d_ind3d(i1,i2,i3,*L1,*L2);
						if(Y[ind_current] != Y[ind_current])
						{
							NAs += 1;
							//sorted[i3] = 1.0 / 0;
							sorted[i3] = 999999999999999999;
						}
						else
						{
							sorted[i3] = Y[ind_current];
						}
					}
					// Sort:
					std::sort(sorted, sorted + *L3);
					// Get the median:
					YSmooth[ind_smooth] = median3d_getMedianOfSorted(sorted, *L3, NAs);
				} // End of for i1
			} // End of for i2
		} // End of *margin=3
	//delete [] sorted;
	//delete ind_smooth;
	//delete ind_current;
	//delete NAs;
	} // End of void
} // End of extern "C"
// 
