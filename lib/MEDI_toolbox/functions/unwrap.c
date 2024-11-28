/*=================================================================
 * 3d unwrapping algorithm
 * Similar to R. Cusack, N. Papadakis NeuroImage 16, 754 - 764 (2002)
 * Use intensity guided region growing algorithm;
 * Starts from a seed point, grow(unwrap) the edge. 
 * Always grow the edge point with the greatest intensity. 
 * Algorithm is slow
 *=================================================================*/
/* $Revision: 1.5.6.2 $ */

// 2013-10-10 LJQ 路径按照ntensity guided region growing algorithm
//                每次pixel unwrap时的比较点为第一次作为邻居的点

//2013-10-10 LJQ　Try to Get the path with the quality map 
//                 Magnetic Resonance in Medicine 62:1085-1090 (2009)   
//                 文献中要求second-order partial derivatives 
//                 好像不能work, 待仔细确认     

#include "unwrap.h"


int main(int argc, char **argv)
{
    int dim[3];
	int numel;
	double *phase, *mag;
	char ifname[256];				/* input filename		*/
	char ofname[256];				/* output filename		*/
	point * MRimage;
	int max_ind = 0;					/* the index of the element with largest intensity */
	int i;
	double noise_level = 0;
	double noise_ratio = 5e-3;
	
	FILE *fid_in, *fid_out;					/* input file ID	*/
	
	//2010-08-16 LJQ modifify the codes for float data
	float *fTempData;
    
	//2010-07-09 LJQ add Segment
	/*short *seg;*/

	if (argc >= 6)
	{
		dim[0] = atoi(argv[3]);
		dim[1] = atoi(argv[4]);
		dim[2] = atoi(argv[5]);
		numel =  dim[0]*dim[1]*dim[2];
		mag = (double *) calloc( numel, sizeof(double));
		phase = (double *) calloc( numel, sizeof(double));
		

		//2010-08-16
		fTempData = (float *) calloc( numel, sizeof(float));
		fid_in = fopen(argv[1], "rb");		/* read in the magnitude map		*/
		if (fid_in)
		{
			fread(fTempData, sizeof(float), numel, fid_in);	
			fclose(fid_in);
		}
		else
		{
			printf(" %s doesn't exist\n", ifname);
			return -1;
		}

		for(i=0; i<numel; i++)
			mag[i] = (double)fTempData[i];
		

		
		fid_in = fopen(argv[2], "rb");		/* read in the phase map */
		if (fid_in)
		{
			fread(fTempData, sizeof(float), numel, fid_in);		
			fclose(fid_in);
		}
		else
		{
			printf(" %s doesn't exist\n", ifname);
			return -1;
		}
		
		for(i=0; i<numel; i++)
			phase[i] = (double)fTempData[i];

		//2010-07-09 LJQ Segment
		/*
		phase = (short *) calloc( numel, sizeof(short));
		fid_in = fopen("seg.img, "rb");		
		if (fid_in)
		{
			fread(seg, sizeof(short), numel, fid_in);	
			fclose(fid_in);
		}
		else
		{
			printf(" %s doesn't exist\n", "segment file");
			return -1;
		}
        */


		//if (argc == 7) noise_ratio = atof(argv[6]);
		if (argc == 7) noise_level = atof(argv[6]);
		//2010-03-07 LJQ modified the code because noise_level is obtained in "step1_read_files.m"
		
		//2010-03-07 LJQ
		printf("noise_level = % 10.5f\n",noise_level);
	}
	else
	{
		printf("syntax: unwrap magnitude_map phase_map dimension1 dimension2 dimension3\n");
		return -1;
	}

	
	MRimage = assemble_data(mag, phase, dim, &max_ind);

	//noise_level = MRimage[max_ind].mag * noise_ratio;
	//2010-03-07 LJQ modified the code because noise_level is obtained in "step1_read_files.m"

	unwrap(MRimage, dim[2]/2*dim[1]*dim[0] + dim[1]/2*dim[0] + dim[0]/2,  dim, noise_level);	/* choose the middle of the FOV as the seeding point	*/
	
	for (i= 0; i<numel; i++)
		phase[i] = MRimage[i].phase;


	for(i=0; i<numel; i++)
		fTempData[i] = (float)phase[i];

	sprintf(ofname, "%s.unwrapped", argv[2]);				/* generate the output filename	*/
	fid_out = fopen(ofname, "wb+");
		fwrite(fTempData, sizeof(float), numel, fid_out);
	fclose(fid_out);

	for (i= 0; i<numel; i++)
		fTempData[i] =(float) MRimage[i].phase_quality;

	sprintf(ofname, "%s.phasequality", argv[2]);				/* generate the output filename	*/
	fid_out = fopen(ofname, "wb+");
		fwrite(fTempData, sizeof(float), numel, fid_out);
	fclose(fid_out);
		

	free(MRimage); MRimage = NULL;
	free(mag); mag = NULL;
	free(phase); phase = NULL;
	free(fTempData); fTempData = NULL;

	return 0;

}

static void swap(int *a, int *b)
{
int c;
	c = *a;
	*a = *b;
	*b = c;
}

static point* assemble_data( double *mag, double *phase, int *dim, int *max_ind)
{
int ind, numel;

point * ret_image;

double tempmean_0,tempmean_1,tempmean_2;

	numel = dim[0]*dim[1]*dim[2];
	ret_image = (point *) calloc(numel, sizeof(point));

	for (ind=0; ind< numel; ind++)
	{
    	ret_image[ind].pos[2]=(int)((ind)/(dim[1]*dim[0]));
    	ret_image[ind].pos[1]=(int)((ind - ret_image[ind].pos[2] *dim[1]*dim[0])/dim[0]);
    	ret_image[ind].pos[0]=(int)(ind - ret_image[ind].pos[1]*dim[0] - ret_image[ind].pos[2]*dim[1]*dim[0]);
		ret_image[ind].ind = ind;
		ret_image[ind].mag = mag[ind];
		ret_image[ind].phase = phase[ind];
		ret_image[ind].planned = 0;
		ret_image[ind].pre = 0;
		ret_image[ind].phase_derivative[0] =0;
		ret_image[ind].phase_derivative[1] =0;
		ret_image[ind].phase_derivative[2] =0;
		ret_image[ind].phase_derivative_2nd[0] =0;
		ret_image[ind].phase_derivative_2nd[1] =0;
		ret_image[ind].phase_derivative_2nd[2] =0;
		ret_image[ind].phase_quality = 100.0;
		if (mag[ind] > mag[*max_ind]) *max_ind = ind;
	}	
	
	//calculate the first order partial phase derivations
	for (ind=0; ind< numel; ind++)
	{
		if(ret_image[ind].pos[0] < (dim[0]-1))
		{
			ret_image[ind].phase_derivative[0] = ret_image[ind+1].phase - ret_image[ind].phase; 
		}

		if(ret_image[ind].pos[1] < (dim[1]-1))
		{
			ret_image[ind].phase_derivative[1] = ret_image[ind+dim[0]].phase - ret_image[ind].phase; 
		}

		if(ret_image[ind].pos[2] < (dim[2]-1))
		{
			ret_image[ind].phase_derivative[2] = ret_image[ind+dim[0]*dim[1]].phase - ret_image[ind].phase; 
		}
	}

	//calculate the second order partial phase derivations
	for (ind=0; ind< numel; ind++)
	{
		if(ret_image[ind].pos[0] < (dim[0]-1))
		{
			ret_image[ind].phase_derivative_2nd[0] = ret_image[ind+1].phase_derivative[0] - ret_image[ind].phase_derivative[0]; 
		}

		if(ret_image[ind].pos[1] < (dim[1]-1))
		{
			ret_image[ind].phase_derivative_2nd[1] = ret_image[ind+dim[0]].phase_derivative[1] - ret_image[ind].phase_derivative[1]; 
		}

		if(ret_image[ind].pos[2] < (dim[2]-1))
		{
			ret_image[ind].phase_derivative_2nd[2] = ret_image[ind+dim[0]*dim[1]].phase_derivative[2] - ret_image[ind].phase_derivative[2]; 
		}
	}

	//calculate the phase quality
	for (ind=0; ind< numel; ind++)
	{
		if((ret_image[ind].pos[0] > 0) && (ret_image[ind].pos[0] < (dim[0]-1)) && (ret_image[ind].pos[1] > 0) && (ret_image[ind].pos[1] < (dim[1]-1))&& (ret_image[ind].pos[2] > 0) && (ret_image[ind].pos[2] < (dim[2]-1)))
		{
			 tempmean_0 = ( ret_image[ind-1].phase_derivative_2nd[0] + ret_image[ind].phase_derivative_2nd[0] + ret_image[ind+1].phase_derivative_2nd[0] )/3.0;
			 tempmean_1 = ( ret_image[ind-dim[0]].phase_derivative_2nd[1] + ret_image[ind].phase_derivative_2nd[1] + ret_image[ind+dim[0]].phase_derivative_2nd[1] )/3.0;
			 tempmean_2 = ( ret_image[ind-dim[0]*dim[1]].phase_derivative_2nd[2] + ret_image[ind].phase_derivative_2nd[2] + ret_image[ind+dim[0]*dim[1]].phase_derivative_2nd[2] )/3.0;
			 ret_image[ind].phase_quality =
				  sqrt(pow(ret_image[ind-1].phase_derivative_2nd[0] -tempmean_0,2) + pow(ret_image[ind].phase_derivative_2nd[0] -tempmean_0,2) + pow(ret_image[ind+1].phase_derivative_2nd[0] -tempmean_0,2))
				+ sqrt(pow(ret_image[ind-dim[0]].phase_derivative_2nd[1] -tempmean_1,2) + pow(ret_image[ind].phase_derivative_2nd[1] -tempmean_1,2) + pow(ret_image[ind+dim[0]].phase_derivative_2nd[1] -tempmean_1,2))
			 	+ sqrt(pow(ret_image[ind-dim[0]*dim[1]].phase_derivative_2nd[2] -tempmean_2,2) + pow(ret_image[ind].phase_derivative_2nd[2] -tempmean_2,2) + pow(ret_image[ind+dim[0]*dim[1]].phase_derivative_2nd[2] -tempmean_2,2));
		}
	}


	return ret_image;

}






static void get_neighbors(int *nei,int *pos,int *dim)
{
    
    int x,y,z;


	x = pos[0];
	y = pos[1];
	z = pos[2];

	if(x>0)
        nei[0]=x-1+y*dim[0]+z*dim[0]*dim[1];
    else
        nei[0]=-1;
    
    if (x<(dim[0]-1))
        nei[1]=x+1+y*dim[0]+z*dim[0]*dim[1];
    else
        nei[1]=-1;

    if (y>0)
        nei[2]=x+(y-1)*dim[0]+z*dim[0]*dim[1];
    else
        nei[2]=-1;
    
    if (y<(dim[1]-1))
        nei[3]=x+(1+y)*dim[0]+z*dim[0]*dim[1];
    else
        nei[3]=-1;

    if (z>0)
        nei[4]=x+y*dim[0]+(z-1)*dim[0]*dim[1];
    else
        nei[4]=-1;
    
    if (z<(dim[2]-1))
        nei[5]=x+y*dim[0]+(z+1)*dim[0]*dim[1];
    else
        nei[5]=-1;

}






static double round(double num)
{
    if (num>=0)
    {
        if ((num-floor(num))>=0.5)
            return floor(num)+1;
        else
            return floor(num);
    }
    else
    {
        if ((num-floor(num))<=0.5)
            return floor(num);
        else 
            return floor(num)+1;
    }
}



static int * heap_create(int seed, int * dim)
{
int numel;
int *heap;
	numel = dim[0]*dim[1]*dim[2];
	heap = (int *) calloc(numel, sizeof(int));
	heap[0] = seed;
	return heap;
}




static void heap_add(point *ori_image, int * heap, int index, int * length)
{
int ind, parent, goon;
	
	heap[*length] = index;
	goon = 1;
	
	ind = *length;
	while ((ind>0) && (goon))
	{
		parent = (ind-1)/2;
	
		if (ori_image[heap[parent]].mag < ori_image[heap[ind]].mag) 
		{
			swap( &heap[parent], &heap[ind]);
			ind = parent;
		}
		else
		{
			goon = 0;
		}
	}
	*length = *length + 1;
}


static int heap_pop( point * ori_image, int * heap, int * length)
{
int left, right, ind, goon, child;
int ret_val;	
	*length = *length - 1;
	ret_val = heap[0];
	heap[0] = heap[*length];
	heap[*length] = 0;
	
	ind = 0;
	goon = 1;
	while ((ind*2+1<= *length - 1) && (goon))
	{
		left = ind * 2 +1;
		right = ind*2 +2;
		if ( ori_image[heap[left]].mag > ori_image[heap[right]].mag )
			child = left;
		else
			child = right;
		if (child > (*length - 1)) child = left;
		
		if (ori_image[heap[ind]].mag < ori_image[heap[child]].mag)
		{
			swap( &heap[ind], &heap[child]);
			ind = child;
		}
		else
		{
			goon = 0;
		}
	}
	return ret_val;
}

static void heap_add_phasequality(point *ori_image, int * heap, int index, int * length)
{
int ind, parent, goon;
	
	heap[*length] = index;
	goon = 1;
	
	ind = *length;
	while ((ind>0) && (goon))
	{
		parent = (ind-1)/2;
	
		if (ori_image[heap[parent]].phase_quality > ori_image[heap[ind]].phase_quality) 
		{
			swap( &heap[parent], &heap[ind]);
			ind = parent;
		}
		else
		{
			goon = 0;
		}
	}
	*length = *length + 1;
}


static int heap_pop_phasequality( point * ori_image, int * heap, int * length)
{
int left, right, ind, goon, child;
int ret_val;	
	*length = *length - 1;
	ret_val = heap[0];
	heap[0] = heap[*length];
	heap[*length] = 0;
	
	ind = 0;
	goon = 1;
	while ((ind*2+1<= *length - 1) && (goon))
	{
		left = ind * 2 +1;
		right = ind*2 +2;
		if ( ori_image[heap[left]].phase_quality < ori_image[heap[right]].phase_quality )
			child = left;
		else
			child = right;
		if (child > (*length - 1)) child = left;
		
		if (ori_image[heap[ind]].phase_quality > ori_image[heap[child]].phase_quality)
		{
			swap( &heap[ind], &heap[child]);
			ind = child;
		}
		else
		{
			goon = 0;
		}
	}
	return ret_val;
}

void unwrap(point *ori_image, int seed, int *dim, double noise)
{
int nei[6];
int curr_ind, counter;
int i;
int * heap;
int heap_length = 0;

/*FILE *fp;*/

	counter = 0;
	curr_ind = seed;

	/* beginning of the searching algorithm	*/
    
	heap = heap_create(seed, dim);
	heap_length++;
	ori_image[heap[0]].pre = -1;
	ori_image[heap[0]].planned = 1;

	/*fp = fopen("position.txt","w"); */

	do
	{
       
		/*
		if(counter >3&&counter<1000)
		{
	    fprintf(fp,"counter= %d\n",counter);
		fprintf(fp,"heap_length= %d\n",heap_length);
		fprintf(fp,"position  = %10d  %10d %10d\n",ori_image[ori_image[ ori_image[ori_image[heap[0]].pre].pre ].pre].pos[0],ori_image[ori_image[ ori_image[ori_image[heap[0]].pre].pre ].pre].pos[1],ori_image[ori_image[ ori_image[ori_image[heap[0]].pre].pre ].pre].pos[2]);
		fprintf(fp,"position  = %10d  %10d %10d\n",ori_image[ ori_image[ori_image[heap[0]].pre].pre ].pos[0],ori_image[ ori_image[ori_image[heap[0]].pre].pre ].pos[1],ori_image[ ori_image[ori_image[heap[0]].pre].pre ].pos[2]);
		fprintf(fp,"position  = %10d  %10d %10d\n",ori_image[ori_image[heap[0]].pre].pos[0],ori_image[ori_image[heap[0]].pre].pos[1],ori_image[ori_image[heap[0]].pre].pos[2]);
		fprintf(fp,"position  = %10d  %10d %10d\n",ori_image[heap[0]].pos[0],ori_image[heap[0]].pos[1],ori_image[heap[0]].pos[2]);
		
		}
		*/
		
		counter ++; 
/*		
		if (ori_image[heap[0]].phase < 10)
			ori_image[heap[0]].phase = 10;
		else 
			ori_image[heap[0]].phase = ori_image[heap[0]].phase + 10;
*/		
		
	
		if (ori_image[heap[0]].pre >=0) 
			ori_image[heap[0]].phase = ori_image[heap[0]].phase - round( ( ori_image[heap[0]].phase - ori_image[ ori_image[heap[0]].pre ].phase )/(2*PI) )   *2*PI;
	
		
        /*2011-12-12 LJQ Replace the above codes with the following codes, 为了避免个别奇异点的影响,采用连续三个点判断********/
		
		/*
		if (ori_image[heap[0]].pre >=0 && ori_image[ori_image[heap[0]].pre].pre >=0 && ori_image[ori_image[ori_image[heap[0]].pre].pre].pre >=0) 
		{
			if(fabs(round( ( ori_image[heap[0]].phase - ori_image[ ori_image[heap[0]].pre ].phase )/(2*PI) ))>0 && fabs(round( ( ori_image[heap[0]].phase - ori_image[ ori_image[ori_image[heap[0]].pre].pre ].phase )/(2*PI) ))>0 && fabs(round( ( ori_image[heap[0]].phase -ori_image[ori_image[ ori_image[ori_image[heap[0]].pre].pre ].pre].phase )/(2*PI) ))>0)
			{
				ori_image[heap[0]].phase = ori_image[heap[0]].phase - round( ( ori_image[heap[0]].phase - ori_image[ ori_image[heap[0]].pre ].phase )/(2*PI) )   *2*PI;
			}
			else
			{
				ori_image[heap[0]].phase = ori_image[heap[0]].phase;
			}

		}*/
		
		get_neighbors(nei, ori_image[heap[0]].pos, dim);

		//curr_ind = heap_pop( ori_image, heap, &heap_length);
		curr_ind = heap_pop_phasequality( ori_image, heap, &heap_length);

		for (i = 0; i<6; i++)
			if ( (nei[i]>=0) && (0 == ori_image[nei[i]].planned ) && (ori_image[nei[i]].mag > noise) )
			{
				
					ori_image[nei[i]].pre = curr_ind;
					//heap_add(ori_image, heap, nei[i], &heap_length);
					heap_add_phasequality(ori_image, heap, nei[i], &heap_length);
				
					ori_image[nei[i]].planned = 1;

									
			}
			

		
	} while (heap_length>0 );

	
	printf("counter = %d\n", counter);
	free(heap); heap = NULL;

}
