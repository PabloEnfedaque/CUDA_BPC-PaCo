#include <time.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"



void CUDA_BPC(int, DATATYPE* , int , int, int*);

#include <stdlib.h>
#include <string.h>



char* concat(const char *s1, const char *s2)
{
    char *result = (char*)malloc(strlen(s1)+strlen(s2)+1);//+1 for the zero-terminator
    strcpy(result, s1);
    strcat(result, s2);
    return result;
}


void load_CBlock(int* samples, char* fileName, int Input_xSize, int Input_ySize, int CBlock_coordinates, int CBlock_xSize, int CBlock_ySize){
	
	FILE* file = fopen(fileName, "rb");
	
	if(file != NULL){
		char* line = NULL;
		int c;
		size_t n = 0;
		getline(&line, &n, file);

		fscanf(file, "*%s %d *%d ", &c);				
		
		for(int i=0;i<CBlock_ySize;i++){
			for(int j=0;j<CBlock_xSize;j++){
		
				fscanf(file, "%d ", &c);
				samples[j+(i*Input_xSize)+CBlock_coordinates] = c;
			}
			
		}
			
	}else printf("\n		ERROR - Wrong input file %s \n\n", fileName);
	fclose(file);
	
}


void load_CStream(int* samples, char* fileName, int CStream_coordinates, int CBlock_xSize, int CBlock_ySize){
	
	FILE* file = fopen(fileName, "rb");
	
	if(file != NULL){			
		
		for(int i=0;i<CBlock_ySize*CBlock_xSize;i++){
				
			int tmp1 = fgetc(file) << 8;
			int tmp2 = fgetc(file);
			samples[i+CStream_coordinates] = tmp1 | tmp2;
					
		}	
			
	}else printf("\n		ERROR - Wrong input file %s \n\n", fileName);
	fclose(file);
	
}


void load_Data(int* input, int* output, int dwt_lvls, int Input_xSize, int Input_ySize, int NCBlocks, int CBlock_xSize, int CBlock_ySize, int CBlockDataLoad, int CStreamDataLoad, const char* folderCB, const char* folderCS){	
	
		
	char file[50];
	int aux = 0;
	int CB_index=0;
	int	subb_xsize = 0;
	int	subb_ysize = 0;

	int	xcoordinate = 0;
	int	ycoordinate = 0;

	subb_xsize = Input_xSize >> dwt_lvls;
	subb_ysize = Input_ySize >> dwt_lvls;

	for(int j=0; j < subb_ysize; j = j + CBlock_ySize)
	{
		for(int i=0; i < subb_xsize; i = i + CBlock_xSize)
		{

			if(CStreamDataLoad==1){
				sprintf(file, concat(folderCS,"%d.bin"), CB_index);

				aux= (((xcoordinate+i)/CBlock_xSize) + (((ycoordinate+j)/CBlock_ySize))*(Input_xSize/CBlock_xSize))*(CBlock_xSize*CBlock_ySize);
				load_CStream(output, file, aux, CBlock_xSize, CBlock_ySize);
			}

			if(CBlockDataLoad==1){

				sprintf(file, concat(folderCB,"%d.txt"), CB_index);
				load_CBlock(input, file, Input_xSize, Input_ySize, xcoordinate + (ycoordinate*Input_xSize) + i + (j*Input_xSize), CBlock_xSize, CBlock_ySize);				
			}

			CB_index++;
			
		}
	}


	for(int decomposition_lvl= dwt_lvls; decomposition_lvl>0 ; decomposition_lvl--)
	{
		
		subb_xsize = Input_xSize >> decomposition_lvl;
		subb_ysize = Input_ySize >> decomposition_lvl;

		xcoordinate = subb_xsize;
		ycoordinate = 0;

		for(int j=0; j < subb_ysize; j = j + CBlock_ySize)
		{
			for(int i=0; i < subb_xsize; i = i + CBlock_xSize)
			{
				if(CStreamDataLoad==1){
					sprintf(file, concat(folderCS,"%d.bin"), CB_index);
					aux= (((xcoordinate+i)/CBlock_xSize) + (((ycoordinate+j)/CBlock_ySize))*(Input_xSize/CBlock_xSize))*(CBlock_xSize*CBlock_ySize);
					load_CStream(output, file, aux, CBlock_xSize, CBlock_ySize);
				}
				if(CBlockDataLoad==1){
					sprintf(file, concat(folderCB,"%d.txt"), CB_index);
					load_CBlock(input, file, Input_xSize, Input_ySize, xcoordinate+(ycoordinate*Input_xSize) + i + (j*Input_xSize), CBlock_xSize, CBlock_ySize);
				}
				CB_index++;
				
			}
		}
		
		xcoordinate = 0;
		ycoordinate = subb_ysize;

		for(int j=0; j < subb_ysize; j = j + CBlock_ySize)
		{
			for(int i=0; i < subb_xsize; i = i + CBlock_xSize)
			{
				if(CStreamDataLoad==1){
					sprintf(file, concat(folderCS,"%d.bin"), CB_index);
					aux= (((xcoordinate+i)/CBlock_xSize) + (((ycoordinate+j)/CBlock_ySize))*(Input_xSize/CBlock_xSize))*(CBlock_xSize*CBlock_ySize);
					load_CStream(output, file, aux, CBlock_xSize, CBlock_ySize);
				}
				if(CBlockDataLoad==1){
					sprintf(file, concat(folderCB,"%d.txt"), CB_index);
					load_CBlock(input, file, Input_xSize, Input_ySize, xcoordinate+(ycoordinate*Input_xSize) + i + (j*Input_xSize), CBlock_xSize, CBlock_ySize);
				}
				CB_index++;
				
			}
		}

		xcoordinate = subb_xsize;
		ycoordinate = subb_ysize;

		for(int j=0; j < subb_ysize; j = j + CBlock_ySize)
		{
			for(int i=0; i < subb_xsize; i = i + CBlock_xSize)
			{
				if(CStreamDataLoad==1){
					sprintf(file, concat(folderCS,"%d.bin"), CB_index);
					aux= (((xcoordinate+i)/CBlock_xSize) + (((ycoordinate+j)/CBlock_ySize))*(Input_xSize/CBlock_xSize))*(CBlock_xSize*CBlock_ySize);
					load_CStream(output, file, aux, CBlock_xSize, CBlock_ySize);
				}
				if(CBlockDataLoad==1){
					sprintf(file, concat(folderCB,"%d.txt"), CB_index);
					load_CBlock(input, file, Input_xSize, Input_ySize, xcoordinate+(ycoordinate*Input_xSize) + i + (j*Input_xSize), CBlock_xSize, CBlock_ySize);
				}
				CB_index++;
			}
		}
	}

}


void check_Image(DATATYPE* input1, DATATYPE* input2, int Input_xSize, int Input_ySize, long int *correct_flag, long int *wrong_count, int *correct_CB){



	for(int i = 0; i< (Input_xSize*Input_ySize) ; i++)
	{	

		if(((DATATYPE)input1[i]) != ((int)input2[i])){

			if(*correct_flag==-1){
					*correct_flag = i;
					*correct_CB = 0;
					*wrong_count= (*wrong_count) +1 ;
			
				}
			else ++(*wrong_count);
		}	

	}
	
}


void Print_Results(int xSize,int ySize,long int correct_flag, long int wrong_count, int correct_CB, int levels, int direction){

	if(direction==0){
		if(correct_flag==-1)	printf("TEST	ENCODER	input size: %dx%d	data has %d DWT levels	>>>	OK\n", xSize, ySize, levels);
		else printf("TEST	ENCODER	input size: %dx%d	data has %d DWT levels	>>>	ERROR ...... first error position:  CodeBlock nº: %d row %ld, column %ld //  number mismatches: %ld (out of %d)\n", xSize, ySize, levels, correct_CB, correct_flag/xSize, correct_flag%xSize, wrong_count, xSize*ySize);

	}
	else{
		if(correct_flag==-1)	printf("TEST	DECODER	input size: %dx%d	data has %d DWT levels	>>>	OK\n", xSize, ySize, levels);
		else printf("TEST	DECODER	input size: %dx%d	data has %d DWT levels	>>>	ERROR ...... first error position:  CodeBlock nº: %d row %ld, column %ld //  number mismatches: %ld (out of %d)\n", xSize, ySize, levels, correct_CB, correct_flag/xSize, correct_flag%xSize, wrong_count, xSize*ySize);

	}
}

void load_LUTs(int** LUT, const char* folder)
{

	int* LUT_h;

	char fileName[50]; 

	FILE* file; 

	int c = -1;
	int max_bitplane = -1;

	sprintf(fileName, concat(folder,"ref.txt"));
	file = fopen(fileName, "rb");	

	if(file != NULL){

		while (fscanf(file, "%*d %*d %d : %*d ", &c)!=EOF){
			if( c > max_bitplane ) max_bitplane = c;
		}
			
	}else printf("\n		ERROR - Wrong input file %s \n\n", fileName);

	sprintf(fileName, concat(folder,"sig.txt"));
	file = fopen(fileName, "rb");	

	if(file != NULL){

		while (fscanf(file, "%*d %*d %d : %*d %*d %*d %*d %*d %*d %*d %*d %*d ", &c)!=EOF){
			if( c > max_bitplane ) max_bitplane = c;
		}
		
	}else printf("\n		ERROR - Wrong input file %s \n\n", fileName);

	sprintf(fileName, concat(folder,"sign.txt"));
	file = fopen(fileName, "rb");	

	if(file != NULL){

		while (fscanf(file, "%*d %*d %d : %*d %*d %*d %*d ", &c)!=EOF){
			if( c > max_bitplane ) max_bitplane = c;
		}
		
	}else printf("\n		ERROR - Wrong input file %s \n\n", fileName);

	max_bitplane++;

	int refinement_N = 	((LUT_N_SUBBANDS) * (max_bitplane) * (N_CONTEXT_REFINEMENT) * (DWT_LEVELS)) + (max_bitplane*N_CONTEXT_REFINEMENT);
	int significance_N = 	((LUT_N_SUBBANDS) * (max_bitplane) * (N_CONTEXT_SIGNIFICANCE) * (DWT_LEVELS)) + (max_bitplane*N_CONTEXT_SIGNIFICANCE);  
	int sign_N = 		((LUT_N_SUBBANDS) * (max_bitplane) * (N_CONTEXT_SIGN) * (DWT_LEVELS)) + (max_bitplane*N_CONTEXT_SIGN);    

	#if 	CPASS_3 == 0

		LUT_h = (int*)malloc((refinement_N + significance_N + sign_N) * sizeof(int));
	
	#else

		LUT_h = (int*)malloc((refinement_N + (significance_N*2) + (sign_N*2)) * sizeof(int));

	#endif

	int i = 0;
	int bitplane = 0;
	int prev_bitplane = -1;
	int c_aux[] = {0,0,0,0,0,0,0,0,0};

	sprintf(fileName, concat(folder,"ref.txt"));
	//fileName = "../images/stats/ref.txt";
	file = fopen(fileName, "rb");	

	if(file != NULL){


		while (fscanf(file, "%*d %*d %d : %d ", &bitplane, &c)!=EOF){
			if(bitplane <= prev_bitplane){				
				for(int z = 0; z < (max_bitplane - prev_bitplane)*N_CONTEXT_REFINEMENT; z++) LUT_h[i + (prev_bitplane*N_CONTEXT_REFINEMENT) + z + 1] = 64;
				i += max_bitplane*N_CONTEXT_REFINEMENT;
			} 
			prev_bitplane = bitplane;
			LUT_h[i + (bitplane*N_CONTEXT_REFINEMENT)] = c;
			
		}
			
	}else printf("\n		ERROR - Wrong input file %s \n\n", fileName);


	i=refinement_N;
	bitplane = 0;
	prev_bitplane = -1;

	sprintf(fileName, concat(folder,"sig.txt"));
	//fileName = "../images/stats/sig.txt";
	file = fopen(fileName, "rb");	

	if(file != NULL){

		while (fscanf(file, "%*d %*d %d : %d %d %d %d %d %d %d %d %d ", &bitplane, &(c_aux[0]), &(c_aux[1]), &(c_aux[2]), &(c_aux[3]), &(c_aux[4]), &(c_aux[5]), &(c_aux[6]), &(c_aux[7]),&(c_aux[8]))!=EOF){
			if(bitplane <= prev_bitplane){				
				for(int z = 0; z < (max_bitplane - prev_bitplane)*N_CONTEXT_SIGNIFICANCE; z++) LUT_h[i + (prev_bitplane*N_CONTEXT_SIGNIFICANCE) + z + 1] = 64;
				i += (max_bitplane*N_CONTEXT_SIGNIFICANCE);
			} 
			prev_bitplane = bitplane;
			LUT_h[i + (bitplane*N_CONTEXT_SIGNIFICANCE)] = c_aux[0];
			LUT_h[i + (bitplane*N_CONTEXT_SIGNIFICANCE) + 1] = c_aux[1];
			LUT_h[i + (bitplane*N_CONTEXT_SIGNIFICANCE) + 2] = c_aux[2];
			LUT_h[i + (bitplane*N_CONTEXT_SIGNIFICANCE) + 3] = c_aux[3];
			LUT_h[i + (bitplane*N_CONTEXT_SIGNIFICANCE) + 4] = c_aux[4];
			LUT_h[i + (bitplane*N_CONTEXT_SIGNIFICANCE) + 5] = c_aux[5];
			LUT_h[i + (bitplane*N_CONTEXT_SIGNIFICANCE) + 6] = c_aux[6];
			LUT_h[i + (bitplane*N_CONTEXT_SIGNIFICANCE) + 7] = c_aux[7];
			LUT_h[i + (bitplane*N_CONTEXT_SIGNIFICANCE) + 8] = c_aux[8];
			
	
		}

		
	}else printf("\n		ERROR - Wrong input file %s \n\n", fileName);

	i=refinement_N + significance_N;
	bitplane = 0;
	prev_bitplane = -1;

	sprintf(fileName, concat(folder,"sign.txt"));
	//fileName = "../images/stats/sign.txt";
	file = fopen(fileName, "rb");	

	if(file != NULL){

		while (fscanf(file, "%*d %*d %d : %d %d %d %d ", &bitplane, &(c_aux[0]), &(c_aux[1]), &(c_aux[2]), &(c_aux[3]))!=EOF){
			if(bitplane <= prev_bitplane){				
				for(int z = 0; z < (max_bitplane - prev_bitplane)*N_CONTEXT_SIGN; z++) LUT_h[i + (prev_bitplane*N_CONTEXT_SIGN) + z + 1] = 64;
				i += (max_bitplane*N_CONTEXT_SIGN);
			} 
			prev_bitplane = bitplane;
			LUT_h[i + (bitplane*N_CONTEXT_SIGN)] = c_aux[0];
			LUT_h[i + (bitplane*N_CONTEXT_SIGN) + 1] = c_aux[1];
			LUT_h[i + (bitplane*N_CONTEXT_SIGN) + 2] = c_aux[2];
			LUT_h[i + (bitplane*N_CONTEXT_SIGN) + 3] = c_aux[3];
				
		}

		while (fscanf(file, "%*d %*d %d : %*d %*d %*d %*d ", &c)!=EOF){
			if( c > max_bitplane ) max_bitplane = c;
		}
		
	}else printf("\n		ERROR - Wrong input file %s \n\n", fileName);


	#if 	CPASS_3 == 1

		i= refinement_N + significance_N + sign_N;
		bitplane = 0;
		prev_bitplane = -1;

		sprintf(fileName, concat(folder,"cp_sig.txt"));
		//fileName = "../images/stats/sig.txt";
		file = fopen(fileName, "rb");	

		if(file != NULL){

			while (fscanf(file, "%*d %*d %d : %d %d %d %d %d %d %d %d %d ", &bitplane, &(c_aux[0]), &(c_aux[1]), &(c_aux[2]), &(c_aux[3]), &(c_aux[4]), &(c_aux[5]), &(c_aux[6]), &(c_aux[7]),&(c_aux[8]))!=EOF){
				if(bitplane <= prev_bitplane){				
					for(int z = 0; z < (max_bitplane - prev_bitplane)*N_CONTEXT_SIGNIFICANCE; z++) LUT_h[i + (prev_bitplane*N_CONTEXT_SIGNIFICANCE) + z + 1] = 64;
					i += (max_bitplane*N_CONTEXT_SIGNIFICANCE);
				} 
				prev_bitplane = bitplane;
				LUT_h[i + (bitplane*N_CONTEXT_SIGNIFICANCE)] = c_aux[0];
				LUT_h[i + (bitplane*N_CONTEXT_SIGNIFICANCE) + 1] = c_aux[1];
				LUT_h[i + (bitplane*N_CONTEXT_SIGNIFICANCE) + 2] = c_aux[2];
				LUT_h[i + (bitplane*N_CONTEXT_SIGNIFICANCE) + 3] = c_aux[3];
				LUT_h[i + (bitplane*N_CONTEXT_SIGNIFICANCE) + 4] = c_aux[4];
				LUT_h[i + (bitplane*N_CONTEXT_SIGNIFICANCE) + 5] = c_aux[5];
				LUT_h[i + (bitplane*N_CONTEXT_SIGNIFICANCE) + 6] = c_aux[6];
				LUT_h[i + (bitplane*N_CONTEXT_SIGNIFICANCE) + 7] = c_aux[7];
				LUT_h[i + (bitplane*N_CONTEXT_SIGNIFICANCE) + 8] = c_aux[8];
				
			}
		
		}else printf("\n		ERROR - Wrong input file %s \n\n", fileName);


		i= refinement_N + significance_N + sign_N + significance_N;
		bitplane = 0;
		prev_bitplane = -1;

		sprintf(fileName, concat(folder,"cp_sign.txt"));
		file = fopen(fileName, "rb");	

		if(file != NULL){

			while (fscanf(file, "%*d %*d %d : %d %d %d %d ", &bitplane, &(c_aux[0]), &(c_aux[1]), &(c_aux[2]), &(c_aux[3]))!=EOF){
				if(bitplane <= prev_bitplane){				
					for(int z = 0; z < (max_bitplane - prev_bitplane)*N_CONTEXT_SIGN; z++) LUT_h[i + (prev_bitplane*N_CONTEXT_SIGN) + z + 1] = 64;
					i += (max_bitplane*N_CONTEXT_SIGN);
				} 
				prev_bitplane = bitplane;
				LUT_h[i + (bitplane*N_CONTEXT_SIGN)] = c_aux[0];
				LUT_h[i + (bitplane*N_CONTEXT_SIGN) + 1] = c_aux[1];
				LUT_h[i + (bitplane*N_CONTEXT_SIGN) + 2] = c_aux[2];
				LUT_h[i + (bitplane*N_CONTEXT_SIGN) + 3] = c_aux[3];
				
			}

			while (fscanf(file, "%*d %*d %d : %*d %*d %*d %*d ", &c)!=EOF){
				if( c > max_bitplane ) max_bitplane = c;
			}
		
		}else printf("\n		ERROR - Wrong input file %s \n\n", fileName);	


	#endif


	#if 	CPASS_3 == 0

		HANDLE_ERROR(cudaMalloc ((void**) &(*LUT), (refinement_N + significance_N + sign_N) * sizeof(int)));

		HANDLE_ERROR(cudaMemcpy(*LUT, LUT_h, (refinement_N + significance_N + sign_N) * sizeof(int) , cudaMemcpyHostToDevice));
	
	#else

		HANDLE_ERROR(cudaMalloc ((void**) &(*LUT), (refinement_N + (significance_N*2) + (sign_N*2)) * sizeof(int)));

		HANDLE_ERROR(cudaMemcpy(*LUT, LUT_h, (refinement_N + (significance_N*2) + (sign_N*2)) * sizeof(int) , cudaMemcpyHostToDevice));

	#endif

	free(LUT_h);

}



int main(int argc, char** argv)
{

	int ySize, xSize, CBlock_xSize, CBlock_ySize;
	int levels;
	DATATYPE *input, *original_samples;
	long int correct_flag= -1;
	long int wrong_count= 0;
	int correct_CB = -1;
	const char *folderCB, *folderLUT;

	srand(time(NULL));

	printf("\nPerforming code test...\n\n");		

	printf("------------------------------------------------------------------------------------\n");	

		printf("			TEST BPC\n");	

	printf("------------------------------------------------------------------------------------\n");
	printf("*** Comparing coder versus reference results ***\n");
	printf("------------------------------------------------------------------------------------\n\n");
	printf("\n");

	int* LUT;

#if 	CPASS_3 == 0	

	printf("+++++++++++++++++++++++++++++++\n");
	printf("+++ 2 Coding passes encoder +++\n");
	printf("+++++++++++++++++++++++++++++++\n\n");
	printf("\n");

	folderLUT = "LUTs/2_coding_pass/";

#else

	printf("+++++++++++++++++++++++++++++++\n");
	printf("+++ 3 Coding passes encoder +++\n");
	printf("+++++++++++++++++++++++++++++++\n\n");
	printf("\n");

	folderLUT = "LUTs/3_coding_pass/";

#endif

	xSize = 2048;
	ySize = 2048;

	levels = DWT_LEVELS;

	folderCB = "sample_image/";	

	CBlock_xSize = CBLOCK_WIDTH; 
	CBlock_ySize = CBLOCK_LENGTH;

	HANDLE_ERROR(cudaHostAlloc((void**)&input, (xSize) * (ySize) * sizeof(DATATYPE), cudaHostAllocDefault));

	original_samples = (DATATYPE*)malloc((xSize) * (ySize) * sizeof(DATATYPE));
			
	load_Data(input, 0, levels, xSize, ySize, (xSize/CBlock_xSize)*(ySize/CBlock_ySize), CBlock_xSize, CBlock_ySize, 1, 0, folderCB, 0);

	memcpy(original_samples, input, (xSize) * (ySize) * sizeof(DATATYPE));

	load_LUTs(&LUT, folderLUT);

	CUDA_BPC(CODE, input , xSize , ySize, LUT);

//Decode ---------------------------------------------------------------------

	load_LUTs(&LUT, folderLUT);

	CUDA_BPC(DECODE, input , xSize , ySize, LUT);

	correct_flag=-1;
	wrong_count=0;
	correct_CB= -1;

	check_Image(input, original_samples, xSize, ySize, &correct_flag, &wrong_count, &correct_CB);

	Print_Results(xSize, ySize, correct_flag, wrong_count, correct_CB, levels, 1);

	HANDLE_ERROR(cudaFreeHost(input));
	free(original_samples);
	printf("\n");


	printf("\n");
	printf("------------------------------------------------------------------------------------\n\n");

	HANDLE_ERROR(cudaDeviceReset());

	printf("\n");

		
	return(0);
}
