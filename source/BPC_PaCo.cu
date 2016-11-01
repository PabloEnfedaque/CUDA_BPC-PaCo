#include <math.h>
#include "common.h"

//#define THREAD_DEBUG_ID -1

__device__ int debug = 0;
__device__ int debug_thread = 10;


// -------------------------------------------------------------------------
//DEVICE FUNCTIONS ---------------------------------------------------------
// -------------------------------------------------------------------------

inline __device__ void Read_Coefficients(DATATYPE* Input, int* TData, int TCoordinate, int Input_XSize)
{

	int AUX_TCoordinate = TCoordinate;
	int negative = 0;
	

	for(int i=0;i<CBLOCK_LENGTH;i++)
	{
		TData[i*NELEMENTS_THREAD_X]= Input[AUX_TCoordinate];

		negative = (TData[i*NELEMENTS_THREAD_X]<0)? 1:0;
		TData[i*NELEMENTS_THREAD_X] = (abs(TData[i*NELEMENTS_THREAD_X])<<1) + negative;		


		TData[(i*NELEMENTS_THREAD_X)+1]= Input[AUX_TCoordinate+1];

		negative = (TData[(i*NELEMENTS_THREAD_X)+1]<0)?  1:0;
		TData[(i*NELEMENTS_THREAD_X)+1] = (abs(TData[(i*NELEMENTS_THREAD_X)+1])<<1) + negative;

#if 	CPASS_3 == 1

		//set CP flag to 1 at start
		TData[i*NELEMENTS_THREAD_X]|=(1<<30);
		TData[(i*NELEMENTS_THREAD_X)+1]|=(1<<30);
#endif

		AUX_TCoordinate += Input_XSize;
	}
	
}

inline __device__ void Write_Coefficients(DATATYPE* Output, unsigned int* TData, int TCoordinate, int Input_XSize)
{

	int AUX_TCoordinate = TCoordinate;

	for(int i=0;i<CBLOCK_LENGTH;i++)
	{

		Output[AUX_TCoordinate] = (TData[i*NELEMENTS_THREAD_X]&0x1FFFFFFF)>>1;
		if((TData[i*NELEMENTS_THREAD_X]&1)==1) Output[AUX_TCoordinate] = -Output[AUX_TCoordinate];

		Output[AUX_TCoordinate + 1] = (TData[(i*NELEMENTS_THREAD_X)+1]&0x1FFFFFFF)>>1;
		if((TData[(i*NELEMENTS_THREAD_X)+1]&1)==1) Output[AUX_TCoordinate + 1] = -Output[AUX_TCoordinate + 1];

		AUX_TCoordinate += Input_XSize;
	}
	
}

inline __device__ void Initialize_Coefficients(int* TData)
{

	for(int i=0;i<CBLOCK_LENGTH;i++)
	{
		TData[i*NELEMENTS_THREAD_X]= 0;	
		TData[(i*NELEMENTS_THREAD_X)+1] = 0;

#if 	CPASS_3 == 1

		//set CP flag to 1 at start
		TData[i*NELEMENTS_THREAD_X]|=(1<<30);
		TData[(i*NELEMENTS_THREAD_X)+1]|=(1<<30);
#endif
	}
	
}


inline __device__ void Find_Subband(int* CB_decomposition_level, int* CB_subband, int TCoordinate_X, int TCoordinate_Y, int Input_XSize, int Input_YSize)
{

	int Input_decomposition_levels = DWT_LEVELS;
	int aux = 1;

	while((*CB_decomposition_level)==-1){

		if(aux==Input_decomposition_levels+1){
			*CB_decomposition_level = Input_decomposition_levels;
			*CB_subband = 0;						
		}
		else if((TCoordinate_X >=(Input_XSize>>aux))||(TCoordinate_Y>=(Input_YSize>>aux)))
		{

			*CB_decomposition_level = aux-1;

			if(TCoordinate_X>=(Input_XSize>>aux)) 
				if(TCoordinate_Y>=(Input_YSize>>aux)) 
					*CB_subband = 2;
				else				
					*CB_subband = 0;
			else 
				*CB_subband = 1;
		}
		else aux++;		
	}
}


inline __device__ void Find_MSB(unsigned int* TData, int* MSB, volatile int* shared_buffer)
{
	int max = 0;
	for(int i=0;i<CBLOCK_LENGTH;i++)
		max|=((TData[(i*2)]>>1)|(TData[(i*2)+1]>>1));

	#if 	FERMI == 1

		shared_buffer[((threadIdx.x%THREADBLOCK_SIZE)+1)] = max;
		max|= shared_buffer[((threadIdx.x%THREADBLOCK_SIZE)+1)+1];

		shared_buffer[((threadIdx.x%THREADBLOCK_SIZE)+1)] = max;
		max|= shared_buffer[((threadIdx.x%THREADBLOCK_SIZE)+1)+2];

		shared_buffer[((threadIdx.x%THREADBLOCK_SIZE)+1)] = max;
		max|= shared_buffer[((threadIdx.x%THREADBLOCK_SIZE)+1)+4];

		shared_buffer[((threadIdx.x%THREADBLOCK_SIZE)+1)] = max;
		max|= shared_buffer[((threadIdx.x%THREADBLOCK_SIZE)+1)+8];

		shared_buffer[((threadIdx.x%THREADBLOCK_SIZE)+1)] = max;
		max|= shared_buffer[((threadIdx.x%THREADBLOCK_SIZE)+1)+16];

		shared_buffer[((threadIdx.x%THREADBLOCK_SIZE)+1)] = max;
		max= shared_buffer[((threadIdx.x%THREADBLOCK_SIZE)+1)-(threadIdx.x%WARPSIZE)];
		
	#else

		max|=__shfl_down(max, 1);
		max|=__shfl_down(max, 2);
		max|=__shfl_down(max, 4);
		max|=__shfl_down(max, 8);
		max|=__shfl_down(max, 16);

		max =__shfl(max, 0);		

	#endif

	#if 	CPASS_3 == 1

		//set CP flag to 0 to compute the MSB
		max&=0xDFFFFFFF;

	#endif
		
	*MSB = 32-__ffs(__brev(max)); //Reverse bits and get least significant position

}


inline __device__ int Compute_Context(		
							unsigned int u_l_coeff, 	unsigned int u_coeff, 	unsigned int u_r_coeff, 
							unsigned int l_coeff, 		 			unsigned int r_coeff, 
							unsigned int b_l_coeff, 	unsigned int b_coeff, 	unsigned int b_r_coeff)
{

	return(	(u_r_coeff>>31) + (u_coeff>>31) + (u_l_coeff>>31) + (r_coeff>>31) + 
			(l_coeff>>31) + (b_r_coeff>>31) + (b_coeff>>31) + (b_l_coeff>>31) );	
}


//CONTEXT SIGN FUNCTIONS		-------------------------------------------------------------------------------------------------


inline __device__ unsigned int Compute_Sign_Context(int horizontal, int vertical)
{
	int context;

//Sign contexts are from -4..-1,0,1..4, but they are represented with the sign as first least significant bit

	if(horizontal == 0){
		if(vertical == 0){
			context = 0; 
		}else if(vertical > 0){
			context = 2;
		}else{
			context = 3;
		}
	}else if(horizontal > 0){
		if(vertical == 0){
			context = 4;
		}else if(vertical > 0){
			context = 6;
		}else{
			context = 0;
		}
	}else{
		if(vertical == 0){
			context = 5;
		}else if(vertical > 0){
			context = 1;
		}else{
			context = 7;
		}
	}

	return(context);
}


inline __device__ unsigned int Compute_Sign_Context(						
								int u_coeff, 		
						int l_coeff, 		 	int r_coeff, 
							 	int b_coeff)
{	

	int left = (l_coeff>>31)==0?0:((l_coeff&1)==1?-1:1);
	int right = (r_coeff>>31)==0?0:((r_coeff&1)==1?-1:1);
	int upper = (u_coeff>>31)==0?0:((u_coeff&1)==1?-1:1);
	int bottom = (b_coeff>>31)==0?0:((b_coeff&1)==1?-1:1);

	return( Compute_Sign_Context( left + right, upper + bottom));
}

inline __device__ void Initialize_LUT_pointers(int* LUT_refinement_pointer, int* LUT_significance_pointer , int* LUT_sign_pointer, int CB_decomposition_level, int CB_subband, int MSB)
{

	int LUT_pointer_aux = 0;	

	*LUT_refinement_pointer =  	(CB_decomposition_level * LUT_N_SUBBANDS * LUT_N_BITPLANES *  N_CONTEXT_REFINEMENT) + 
					(CB_subband * LUT_N_BITPLANES *  N_CONTEXT_REFINEMENT) + 
					(MSB*N_CONTEXT_REFINEMENT);	

	LUT_pointer_aux = ((LUT_N_SUBBANDS) * (LUT_N_BITPLANES) * (N_CONTEXT_REFINEMENT) * (DWT_LEVELS)) + (LUT_N_BITPLANES*N_CONTEXT_REFINEMENT);

	*LUT_significance_pointer = 	(CB_decomposition_level * LUT_N_SUBBANDS * LUT_N_BITPLANES *  N_CONTEXT_SIGNIFICANCE) + 
					(CB_subband * LUT_N_BITPLANES *  N_CONTEXT_SIGNIFICANCE) + 
					(MSB*N_CONTEXT_SIGNIFICANCE) + 
					LUT_pointer_aux;
	
	LUT_pointer_aux += ((LUT_N_SUBBANDS) * (LUT_N_BITPLANES) * (N_CONTEXT_SIGNIFICANCE) * (DWT_LEVELS)) + (LUT_N_BITPLANES*N_CONTEXT_SIGNIFICANCE);
	
	*LUT_sign_pointer =  (CB_decomposition_level * LUT_N_SUBBANDS * LUT_N_BITPLANES *  N_CONTEXT_SIGN) + 
				(CB_subband * LUT_N_BITPLANES *  N_CONTEXT_SIGN) + 
				(MSB*N_CONTEXT_SIGN) +
				LUT_pointer_aux;
}


inline __device__ void Update_LUT_pointers(int* LUT_refinement_pointer, int* LUT_significance_pointer , int* LUT_sign_pointer)
{
	*LUT_significance_pointer -=  N_CONTEXT_SIGNIFICANCE;	
	*LUT_sign_pointer -=  N_CONTEXT_SIGN;	
	*LUT_refinement_pointer -=  N_CONTEXT_REFINEMENT;
}


inline __device__ void Arithmetic_encoder(unsigned int Symbol , VOLATILE int* Code_Stream, int Code_Stream_pointer, unsigned int* Reserved_codeword, unsigned int* AC_Interval_Lower, unsigned int* AC_Interval_Size, int Shared_pointer, volatile int* Code_Stream_Shared, int Probability)
{
	
	unsigned int aux = 0;

	//printf("%d ", (threadIdx.x & 0x1f));

	if((*AC_Interval_Size)==0){	
		
		aux= __ballot(1);

		*AC_Interval_Lower	=	0;
		*AC_Interval_Size	=	(1<<CODEWORD_SIZE)-1;
		
		*Reserved_codeword 	= 	__popc(aux<<(WARPSIZE-(threadIdx.x & 0x1f))) + Code_Stream_Shared[ Shared_pointer ]+ Code_Stream_pointer;
		if((aux>>(threadIdx.x & 0x1f))==1)	Code_Stream_Shared[ Shared_pointer ] += __popc(aux);

				
	}

	aux = (((*AC_Interval_Size)*Probability)>>(MULT_PRECISION)) + Symbol;

	if((Symbol==0))		*AC_Interval_Size = aux;
	else{	
		
		*AC_Interval_Size -= aux;
		*AC_Interval_Lower += aux;
	}		

	if((*AC_Interval_Size)==0){
		
		Code_Stream[*Reserved_codeword] = *AC_Interval_Lower;

	}
}

inline __device__ void Arithmetic_decoder(unsigned int* Symbol , VOLATILE int* Code_Stream, int Code_Stream_pointer, unsigned int* Codeword, unsigned int* AC_Interval_Lower, unsigned int* AC_Interval_Size, int Shared_pointer, volatile int* Code_Stream_Shared, int Probability)
{
	
	unsigned int aux = 0;
	unsigned int aux2 = 0;


	if((*AC_Interval_Size)==0){	
		
		aux= __ballot(1);

		*AC_Interval_Lower	=	0;
		*AC_Interval_Size	=	(1<<CODEWORD_SIZE)-1;
		
		aux2  	= 	__popc(aux<<(WARPSIZE-(threadIdx.x & 0x1f))) + Code_Stream_Shared[ Shared_pointer ] + Code_Stream_pointer;

		*Codeword = Code_Stream[ aux2 ];

		if((aux>>(threadIdx.x & 0x1f))==1)	Code_Stream_Shared[ Shared_pointer ] += __popc(aux);	

		//if((threadIdx.x==THREAD_DEBUG_ID)&&(blockIdx.x==0)) printf("Reserved CW: %d ", aux2);	
		
	}

	aux =  (((*AC_Interval_Size)*Probability)>>(MULT_PRECISION)) + 1;
	aux2 = *AC_Interval_Lower + aux;

	if(*Codeword >= aux2)
	{
		*AC_Interval_Size -= aux;
		*AC_Interval_Lower = aux2;

		*Symbol = 1;
	}
	else{

		*AC_Interval_Size = aux - 1;

		*Symbol = 0;
			
	}	
 
}


inline __device__ void Share_Left_Borders(unsigned int* TData, int Coeff_index, volatile int* shared_buffer, int* u_auxiliar, int* m_auxiliar, int* b_auxiliar){

#if 	FERMI == 1

		shared_buffer[((threadIdx.x%THREADBLOCK_SIZE)+1)] = TData[Coeff_index-1];
		*u_auxiliar= shared_buffer[((threadIdx.x%THREADBLOCK_SIZE)+1)-1];

		shared_buffer[((threadIdx.x%THREADBLOCK_SIZE)+1)] = TData[Coeff_index+1];
		*m_auxiliar= shared_buffer[((threadIdx.x%THREADBLOCK_SIZE)+1)-1];

		shared_buffer[((threadIdx.x%THREADBLOCK_SIZE)+1)] = TData[Coeff_index+3];
		*b_auxiliar= shared_buffer[((threadIdx.x%THREADBLOCK_SIZE)+1)-1];

	#else		

		*u_auxiliar=	__shfl_up(TData[Coeff_index-1], 1);
		*m_auxiliar=	__shfl_up(TData[Coeff_index+1], 1);
		*b_auxiliar=	__shfl_up(TData[Coeff_index+3], 1);		


	#endif

}

inline __device__ void Share_Right_Borders(unsigned int* TData, int Coeff_index, volatile int* shared_buffer, int* u_auxiliar, int* m_auxiliar, int* b_auxiliar){

#if 	FERMI == 1
			
		shared_buffer[((threadIdx.x%THREADBLOCK_SIZE)+1)] = TData[Coeff_index-3];
		*u_auxiliar= shared_buffer[((threadIdx.x%THREADBLOCK_SIZE)+1)+1];

		shared_buffer[((threadIdx.x%THREADBLOCK_SIZE)+1)] = TData[Coeff_index-1];
		*m_auxiliar= shared_buffer[((threadIdx.x%THREADBLOCK_SIZE)+1)+1];

		shared_buffer[((threadIdx.x%THREADBLOCK_SIZE)+1)] = TData[Coeff_index+1];
		*b_auxiliar= shared_buffer[((threadIdx.x%THREADBLOCK_SIZE)+1)+1];

	#else

		*u_auxiliar=	__shfl_down(TData[Coeff_index-3], 1);
		*m_auxiliar=	__shfl_down(TData[Coeff_index-1], 1);
		*b_auxiliar=	__shfl_down(TData[Coeff_index+1], 1);		


	#endif

}

inline __device__ void Correct_CB_borders(unsigned int* TData1, unsigned int* TData4, unsigned int* TData6, unsigned int* TData3, unsigned int* TData5, unsigned int* TData8, int direction)
{

	if(direction==1){						
		if((threadIdx.x%32)==0)
		{
			*TData1 = 0;
			*TData4 = 0;
			*TData6 = 0;
		}		 	
	}
	else{
		if((threadIdx.x%32)==31)
		{
			*TData3 = 0;
			*TData5 = 0;
			*TData8 = 0;
		}				
	}
	
}

inline __device__ void SPP_Encoder(unsigned int* TData, int Coeff_index, int Bit_Plane, VOLATILE REG_DATATYPE* Code_Stream, int Code_Stream_pointer, unsigned int* Reserved_codeword, unsigned int* AC_Interval_Lower, unsigned int* AC_Interval_Size, int Shared_pointer, volatile int* Code_Stream_Shared, int direction, int LUT_significance_pointer, int LUT_sign_pointer, volatile int* shared_buffer, int const* __restrict__ LUT, unsigned int TData1, unsigned int TData2, unsigned int TData3, unsigned int TData4, unsigned int TData5, unsigned int TData6, unsigned int TData7, unsigned int TData8)
{	

	unsigned int Context = 0;

	if(!(TData[Coeff_index]>>31))
	{		

		Correct_CB_borders(&TData1, &TData4, &TData6, &TData3, &TData5, &TData8, direction);

#if 	CPASS_3 == 1
	
		if((TData1>>31) || (TData2>>31) || (TData3>>31) || (TData4>>31) || (TData5>>31) || (TData6>>31) || (TData7>>31) || (TData8>>31))
		{
		
#endif					
				
			Context = Compute_Context(	TData1, TData2, TData3, TData4, TData5, TData6, TData7, TData8);

			Arithmetic_encoder((TData[Coeff_index]>>(Bit_Plane+1))&1, Code_Stream, Code_Stream_pointer, Reserved_codeword, AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, LUT[LUT_significance_pointer+Context]);

			
			//if((threadIdx.x==THREAD_DEBUG_ID)&&(blockIdx.x==0)) printf("SSP    -- Coeff: %d Symbol: %d Reserved CW: %d AC_L: %d AC_S: %d Context: %d Prob: %d\n", Coeff_index, (TData[Coeff_index]>>(Bit_Plane+1))&1, *Reserved_codeword, *AC_Interval_Lower, *AC_Interval_Size, Context, LUT[LUT_significance_pointer+Context]);

			if(((TData[Coeff_index]>>(Bit_Plane+1))&1) == 1)
			{ 
			
				TData[Coeff_index]|=(1<<31);

				Context = Compute_Sign_Context( TData2, TData4, TData5, TData7);
	
				Arithmetic_encoder(	((TData[Coeff_index]&1)==(Context&1))?0:1,	Code_Stream, Code_Stream_pointer, Reserved_codeword, AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, LUT[LUT_sign_pointer+(Context>>1)]);

				//if((threadIdx.x==THREAD_DEBUG_ID)&&(blockIdx.x==0)) printf("SPP SIGN    -- Coeff: %d Symbol: %d Reserved CW: %d AC_L: %d AC_S: %d \n", Coeff_index, ((TData[Coeff_index]&1)==(Context&1))?0:1, *Reserved_codeword, *AC_Interval_Lower, *AC_Interval_Size);

				
			}

#if 	CPASS_3 == 1	

		}else TData[Coeff_index]|=(1<<30);

#endif	
	}
}

inline __device__ void SPP_Decoder(unsigned int* TData, int Coeff_index, int Mask, VOLATILE REG_DATATYPE* Code_Stream, int Code_Stream_pointer, unsigned int* Reserved_codeword, unsigned int* AC_Interval_Lower, unsigned int* AC_Interval_Size, int Shared_pointer, volatile int* Code_Stream_Shared, int direction, int LUT_significance_pointer, int LUT_sign_pointer, volatile int* shared_buffer, int const* __restrict__ LUT, unsigned int TData1, unsigned int TData2, unsigned int TData3, unsigned int TData4, unsigned int TData5, unsigned int TData6, unsigned int TData7, unsigned int TData8)
{	

	unsigned int Context = 0;
	unsigned int Symbol = 0;

	if(!(TData[Coeff_index]>>31))
	{		

		Correct_CB_borders(&TData1, &TData4, &TData6, &TData3, &TData5, &TData8, direction);

#if 	CPASS_3 == 1
	
		if((TData1>>31) || (TData2>>31) || (TData3>>31) || (TData4>>31) || (TData5>>31) || (TData6>>31) || (TData7>>31) || (TData8>>31))
		{
		
#endif					
				
			Context = Compute_Context( TData1, TData2, TData3, TData4, TData5, TData6, TData7, TData8);

			Arithmetic_decoder(&Symbol, Code_Stream, Code_Stream_pointer, Reserved_codeword, AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, LUT[LUT_significance_pointer+Context]);

			//if((threadIdx.x==THREAD_DEBUG_ID)&&(blockIdx.x==0)) printf("SSP    -- Coeff: %d Symbol: %d CW: %d AC_L: %d AC_S: %d Context: %d Prob: %d\n", Coeff_index, Symbol, *Reserved_codeword, *AC_Interval_Lower, *AC_Interval_Size, Context, LUT[LUT_significance_pointer+Context]);

			if(Symbol == 1)
			{ 

				TData[Coeff_index]|= Mask;
			
				TData[Coeff_index]|=(1<<31);			

				Context = Compute_Sign_Context( TData2, TData4, TData5, TData7);						
	
				Arithmetic_decoder( &Symbol, Code_Stream, Code_Stream_pointer, Reserved_codeword, AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, LUT[LUT_sign_pointer+(Context>>1)]);

				Symbol=((Symbol&1)==(Context&1))?0:1;	

				//if((threadIdx.x==THREAD_DEBUG_ID)&&(blockIdx.x==0)) printf("SPP SIGN    -- Coeff: %d Symbol: %d CW: %d AC_L: %d AC_S: %d \n", Coeff_index, Symbol, *Reserved_codeword, *AC_Interval_Lower, *AC_Interval_Size);

				TData[Coeff_index]|= (Symbol&1);			

			}

#if 	CPASS_3 == 1	

		}else TData[Coeff_index]|=(1<<30);

#endif	
	}
}


inline __device__ void CP_Encoder(unsigned int* TData, int Coeff_index, int Bit_Plane, VOLATILE REG_DATATYPE* Code_Stream, int Code_Stream_pointer, unsigned int* Reserved_codeword, unsigned int* AC_Interval_Lower, unsigned int* AC_Interval_Size, int Shared_pointer, volatile int* Code_Stream_Shared, int direction, int LUT_significance_pointer, int LUT_sign_pointer, volatile int* shared_buffer, int const* __restrict__ LUT, unsigned int TData1, unsigned int TData2, unsigned int TData3, unsigned int TData4, unsigned int TData5, unsigned int TData6, unsigned int TData7, unsigned int TData8)
{	

	unsigned int Context = 0;

	if((TData[Coeff_index]>>30)&1)
	{		

		Correct_CB_borders(&TData1, &TData4, &TData6, &TData3, &TData5, &TData8, direction);					
			
		Context = Compute_Context( TData1, TData2, TData3, TData4, TData5, TData6, TData7, TData8);


		Arithmetic_encoder((TData[Coeff_index]>>(Bit_Plane+1))&1, Code_Stream, Code_Stream_pointer, Reserved_codeword, AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, LUT[LUT_significance_pointer+Context]);

		//if((threadIdx.x==THREAD_DEBUG_ID)&&(blockIdx.x==0)) printf("CP    -- Coeff: %d Symbol: %d Reserved CW: %d AC_L: %d AC_S: %d \n", Coeff_index, (TData[Coeff_index]>>(Bit_Plane+1))&1, *Reserved_codeword, *AC_Interval_Lower, *AC_Interval_Size);


		TData[Coeff_index]&=0xBFFFFFFF;

		if(((TData[Coeff_index]>>(Bit_Plane+1))&1) == 1)
		{ 
	
			TData[Coeff_index]|=(1<<31);
			TData[Coeff_index]|=(1<<29);			

			Context = Compute_Sign_Context( TData2, TData4, TData5, TData7);

			Arithmetic_encoder(	((TData[Coeff_index]&1)==(Context&1))?0:1,	Code_Stream, Code_Stream_pointer, Reserved_codeword, AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, LUT[LUT_sign_pointer+(Context>>1)]);

			//if((threadIdx.x==THREAD_DEBUG_ID)&&(blockIdx.x==0)) printf("CP SIGN    -- Coeff: %d Symbol: %d Reserved CW: %d AC_L: %d AC_S: %d \n", Coeff_index, ((TData[Coeff_index]&1)==(Context&1))?0:1, *Reserved_codeword, *AC_Interval_Lower, *AC_Interval_Size);

		}
	}
}

inline __device__ void CP_Decoder(unsigned int* TData, int Coeff_index, int Mask, VOLATILE REG_DATATYPE* Code_Stream, int Code_Stream_pointer, unsigned int* Reserved_codeword, unsigned int* AC_Interval_Lower, unsigned int* AC_Interval_Size, int Shared_pointer, volatile int* Code_Stream_Shared, int direction, int LUT_significance_pointer, int LUT_sign_pointer, volatile int* shared_buffer, int const* __restrict__ LUT, unsigned int TData1, unsigned int TData2, unsigned int TData3, unsigned int TData4, unsigned int TData5, unsigned int TData6, unsigned int TData7, unsigned int TData8)
{	

	unsigned int Context = 0;
	unsigned int Symbol = 0;

	if((TData[Coeff_index]>>30)&1)
	{		

		Correct_CB_borders(&TData1, &TData4, &TData6, &TData3, &TData5, &TData8, direction);					
			
		Context = Compute_Context(	TData1, TData2, TData3, TData4, TData5, TData6, TData7, TData8);

		Arithmetic_decoder(&Symbol, Code_Stream, Code_Stream_pointer, Reserved_codeword, AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, LUT[LUT_significance_pointer+Context]);

		//if((threadIdx.x==THREAD_DEBUG_ID)&&(blockIdx.x==0)) printf("CP    -- Coeff: %d Symbol: %d CW: %d AC_L: %d AC_S: %d \n", Coeff_index, Symbol, *Reserved_codeword, *AC_Interval_Lower, *AC_Interval_Size);

		TData[Coeff_index]&=0xBFFFFFFF;

		if(Symbol == 1)
		{ 

			TData[Coeff_index]|= Mask;	

			TData[Coeff_index]|=(1<<31);
			TData[Coeff_index]|=(1<<29);			

			Context = Compute_Sign_Context( TData2, TData4, TData5, TData7);

			Arithmetic_decoder(	&Symbol, Code_Stream, Code_Stream_pointer, Reserved_codeword, AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, LUT[LUT_sign_pointer+(Context>>1)]);

			Symbol=((Symbol&1)==(Context&1))?0:1;	

			//if((threadIdx.x==THREAD_DEBUG_ID)&&(blockIdx.x==0)) printf("CP SIGN    -- Coeff: %d Symbol: %d CW: %d AC_L: %d AC_S: %d \n", Coeff_index, Symbol, *Reserved_codeword, *AC_Interval_Lower, *AC_Interval_Size);


			TData[Coeff_index]|= (Symbol&1);
			

		}			


	}
}



inline __device__ void MRP_Encoder(unsigned int* TData, int Coeff_index, int Bit_Plane, VOLATILE REG_DATATYPE* Code_Stream, int Code_Stream_pointer, unsigned int* Reserved_codeword, unsigned int* AC_Interval_Lower, unsigned int* AC_Interval_Size, int  Shared_pointer, volatile int* Code_Stream_Shared, int Probability)
{	

	if((TData[Coeff_index]>>29)&1)
	{
			
			Arithmetic_encoder((TData[Coeff_index]>>(Bit_Plane+1))&1, Code_Stream, Code_Stream_pointer, Reserved_codeword, AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, Probability);

			//if((threadIdx.x==THREAD_DEBUG_ID)&&(blockIdx.x==0)) printf("REF    -- Coeff: %d Symbol: %d Reserved CW: %d AC_L: %d AC_S: %d \n", Coeff_index, (TData[Coeff_index]>>(Bit_Plane+1))&1, *Reserved_codeword, *AC_Interval_Lower, *AC_Interval_Size);

	}
	else if(TData[Coeff_index]>>31) TData[Coeff_index]|=(1<<29);
 }


inline __device__ void MRP_Decoder(unsigned int* TData, int Coeff_index, int Mask, VOLATILE REG_DATATYPE* Code_Stream, int Code_Stream_pointer, unsigned int* Reserved_codeword, unsigned int* AC_Interval_Lower, unsigned int* AC_Interval_Size, int  Shared_pointer, volatile int* Code_Stream_Shared, int Probability, int Bit_Plane)
{	
	
	unsigned int Symbol = 0;

	if((TData[Coeff_index]>>29)&1)
	{
			

			Arithmetic_decoder(&Symbol, Code_Stream, Code_Stream_pointer, Reserved_codeword, AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, Probability);

			//if((threadIdx.x==THREAD_DEBUG_ID)&&(blockIdx.x==0)) printf("REF    -- Coeff: %d Symbol: %d CW: %d AC_L: %d AC_S: %d \n", Coeff_index, Symbol, *Reserved_codeword, *AC_Interval_Lower, *AC_Interval_Size);

			//Delete previous aproximate bit
			TData[Coeff_index]&= (~Mask);
			//Write new bit and aproximation
			TData[Coeff_index]|= (Mask&(((Symbol<<1)+1)<<(Bit_Plane)));

		
	}
	else if(TData[Coeff_index]>>31) TData[Coeff_index]|=(1<<29);
 }



inline __device__ void SPP_Encoder_launcher(unsigned int* TData, int Bit_Plane, VOLATILE REG_DATATYPE* Code_Stream, int Code_Stream_pointer, unsigned int* Reserved_codeword, unsigned int* AC_Interval_Lower, unsigned int* AC_Interval_Size, int Shared_pointer, volatile int* Code_Stream_Shared, int LUT_significance_pointer, int LUT_sign_pointer, volatile int* shared_buffer, int const* __restrict__ LUT)
{	

	int i = 0;

	int u_auxiliar=0;
	int m_auxiliar=0;
	int b_auxiliar=0;

	int Coeff_index = i*2;

	Share_Left_Borders(TData, Coeff_index, shared_buffer, &u_auxiliar, &m_auxiliar, &b_auxiliar);

	SPP_Encoder(		TData, Coeff_index, Bit_Plane, Code_Stream, Code_Stream_pointer, Reserved_codeword, 
				AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, 1, LUT_significance_pointer, LUT_sign_pointer, shared_buffer, LUT, 
				0,		0, 			0,   
				m_auxiliar,				TData[Coeff_index+1],
				b_auxiliar,	TData[Coeff_index+2],	TData[Coeff_index+3]	);

	Coeff_index = (i*2)+1;

	Share_Right_Borders(TData, Coeff_index, shared_buffer, &u_auxiliar, &m_auxiliar, &b_auxiliar);

	SPP_Encoder(		TData, Coeff_index, Bit_Plane, Code_Stream, Code_Stream_pointer, Reserved_codeword, 						
				AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, 2, LUT_significance_pointer, LUT_sign_pointer, shared_buffer, LUT,
				0,			0, 			0,      
				TData[Coeff_index-1],				m_auxiliar,
				TData[Coeff_index+1],	TData[Coeff_index+2],	b_auxiliar	);

	for(i=1;i<(CBLOCK_LENGTH -1);i++)
	{

		Coeff_index = i*2;
	
		Share_Left_Borders(TData, Coeff_index,shared_buffer, &u_auxiliar, &m_auxiliar, &b_auxiliar);		

		SPP_Encoder(	TData, Coeff_index, Bit_Plane, Code_Stream, Code_Stream_pointer, Reserved_codeword, 
				AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, 1, LUT_significance_pointer, LUT_sign_pointer, shared_buffer, LUT, 
				u_auxiliar,	TData[Coeff_index-2], 	TData[Coeff_index-1],   
				m_auxiliar,				TData[Coeff_index+1],
				b_auxiliar,	TData[Coeff_index+2],	TData[Coeff_index+3]	);

		Coeff_index = (i*2)+1;

		Share_Right_Borders(TData, Coeff_index,shared_buffer, &u_auxiliar, &m_auxiliar, &b_auxiliar);

		SPP_Encoder(	TData, Coeff_index, Bit_Plane, Code_Stream, Code_Stream_pointer, Reserved_codeword, 
				AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, 2, LUT_significance_pointer, LUT_sign_pointer, shared_buffer, LUT,
				TData[Coeff_index-3],	TData[Coeff_index-2], 	u_auxiliar,   
				TData[Coeff_index-1],				m_auxiliar,
				TData[Coeff_index+1],	TData[Coeff_index+2],	b_auxiliar	);
	}

	Coeff_index = i*2;

	Share_Left_Borders(TData, Coeff_index,shared_buffer, &u_auxiliar, &m_auxiliar, &b_auxiliar);

	SPP_Encoder(		TData, Coeff_index, Bit_Plane, Code_Stream, Code_Stream_pointer, Reserved_codeword, 
				AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, 1, LUT_significance_pointer, LUT_sign_pointer, shared_buffer, LUT,
				u_auxiliar,	TData[Coeff_index-2], 	TData[Coeff_index-1],   
				m_auxiliar,				TData[Coeff_index+1],
				0,			0, 			0  	 	);

	Coeff_index = (i*2)+1;

	Share_Right_Borders(TData, Coeff_index,shared_buffer, &u_auxiliar, &m_auxiliar, &b_auxiliar);

	SPP_Encoder(		TData, Coeff_index, Bit_Plane, Code_Stream, Code_Stream_pointer, Reserved_codeword, 
				AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, 2, LUT_significance_pointer, LUT_sign_pointer, shared_buffer, LUT,
				TData[Coeff_index-3],	TData[Coeff_index-2], 	u_auxiliar,   
				TData[Coeff_index-1],				m_auxiliar,
				0,			0, 			0      		);

}


inline __device__ void SPP_Decoder_launcher(unsigned int* TData, int Mask, VOLATILE REG_DATATYPE* Code_Stream, int Code_Stream_pointer, unsigned int* Codeword, unsigned int* AC_Interval_Lower, unsigned int* AC_Interval_Size, int Shared_pointer, volatile int* Code_Stream_Shared, int LUT_significance_pointer, int LUT_sign_pointer, volatile int* shared_buffer, int const* __restrict__ LUT)
{	

	int i = 0;

	int u_auxiliar=0;
	int m_auxiliar=0;
	int b_auxiliar=0;

	int Coeff_index = i*2;

	Share_Left_Borders(TData, Coeff_index, shared_buffer, &u_auxiliar, &m_auxiliar, &b_auxiliar);

	SPP_Decoder(		TData, Coeff_index, Mask, Code_Stream, Code_Stream_pointer, Codeword, 
				AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, 1, LUT_significance_pointer, LUT_sign_pointer, shared_buffer, LUT, 
				0,		0, 			0,   
				m_auxiliar,				TData[Coeff_index+1],
				b_auxiliar,	TData[Coeff_index+2],	TData[Coeff_index+3]	);

	Coeff_index = (i*2)+1;

	Share_Right_Borders(TData, Coeff_index, shared_buffer, &u_auxiliar, &m_auxiliar, &b_auxiliar);

	SPP_Decoder(		TData, Coeff_index, Mask, Code_Stream, Code_Stream_pointer, Codeword, 						
				AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, 2, LUT_significance_pointer, LUT_sign_pointer, shared_buffer, LUT,
				0,			0, 			0,      
				TData[Coeff_index-1],				m_auxiliar,
				TData[Coeff_index+1],	TData[Coeff_index+2],	b_auxiliar	);

	for(i=1;i<(CBLOCK_LENGTH -1);i++)
	{

		Coeff_index = i*2;
	
		Share_Left_Borders(TData, Coeff_index,shared_buffer, &u_auxiliar, &m_auxiliar, &b_auxiliar);		

		SPP_Decoder(	TData, Coeff_index, Mask, Code_Stream, Code_Stream_pointer, Codeword, 
				AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, 1, LUT_significance_pointer, LUT_sign_pointer, shared_buffer, LUT, 
				u_auxiliar,	TData[Coeff_index-2], 	TData[Coeff_index-1],   
				m_auxiliar,				TData[Coeff_index+1],
				b_auxiliar,	TData[Coeff_index+2],	TData[Coeff_index+3]	);

		Coeff_index = (i*2)+1;

		Share_Right_Borders(TData, Coeff_index,shared_buffer, &u_auxiliar, &m_auxiliar, &b_auxiliar);

		SPP_Decoder(	TData, Coeff_index, Mask, Code_Stream, Code_Stream_pointer, Codeword, 
				AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, 2, LUT_significance_pointer, LUT_sign_pointer, shared_buffer, LUT,
				TData[Coeff_index-3],	TData[Coeff_index-2], 	u_auxiliar,   
				TData[Coeff_index-1],				m_auxiliar,
				TData[Coeff_index+1],	TData[Coeff_index+2],	b_auxiliar	);
	}

	Coeff_index = i*2;

	Share_Left_Borders(TData, Coeff_index,shared_buffer, &u_auxiliar, &m_auxiliar, &b_auxiliar);

	SPP_Decoder(		TData, Coeff_index, Mask, Code_Stream, Code_Stream_pointer, Codeword, 
				AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, 1, LUT_significance_pointer, LUT_sign_pointer, shared_buffer, LUT,
				u_auxiliar,	TData[Coeff_index-2], 	TData[Coeff_index-1],   
				m_auxiliar,				TData[Coeff_index+1],
				0,			0, 			0  	 	);

	Coeff_index = (i*2)+1;

	Share_Right_Borders(TData, Coeff_index,shared_buffer, &u_auxiliar, &m_auxiliar, &b_auxiliar);

	SPP_Decoder(		TData, Coeff_index, Mask, Code_Stream, Code_Stream_pointer, Codeword, 
				AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, 2, LUT_significance_pointer, LUT_sign_pointer, shared_buffer, LUT,
				TData[Coeff_index-3],	TData[Coeff_index-2], 	u_auxiliar,   
				TData[Coeff_index-1],				m_auxiliar,
				0,			0, 			0      		);

}

inline __device__ void CP_Encoder_launcher(unsigned int* TData, int Bit_Plane, VOLATILE REG_DATATYPE* Code_Stream, int Code_Stream_pointer, unsigned int* Reserved_codeword, unsigned int* AC_Interval_Lower, unsigned int* AC_Interval_Size, int Shared_pointer, volatile int* Code_Stream_Shared, int LUT_significance_pointer, int LUT_sign_pointer, volatile int* shared_buffer, int const* __restrict__ LUT)
{	

	int i = 0;

	int u_auxiliar=0;
	int m_auxiliar=0;
	int b_auxiliar=0;

	int Coeff_index = i*2;

	Share_Left_Borders(TData, Coeff_index, shared_buffer, &u_auxiliar, &m_auxiliar, &b_auxiliar);

	CP_Encoder(		TData, Coeff_index, Bit_Plane, Code_Stream, Code_Stream_pointer, Reserved_codeword, 
				AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, 1, LUT_significance_pointer, LUT_sign_pointer, shared_buffer, LUT, 
				0,		0, 			0,   
				m_auxiliar,				TData[Coeff_index+1],
				b_auxiliar,	TData[Coeff_index+2],	TData[Coeff_index+3]	);

	Coeff_index = (i*2)+1;

	Share_Right_Borders(TData, Coeff_index, shared_buffer, &u_auxiliar, &m_auxiliar, &b_auxiliar);

	CP_Encoder(		TData, Coeff_index, Bit_Plane, Code_Stream, Code_Stream_pointer, Reserved_codeword, 						
				AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, 2, LUT_significance_pointer, LUT_sign_pointer, shared_buffer, LUT,
				0,			0, 			0,      
				TData[Coeff_index-1],				m_auxiliar,
				TData[Coeff_index+1],	TData[Coeff_index+2],	b_auxiliar	);

	for(i=1;i<(CBLOCK_LENGTH -1);i++)
	{

		Coeff_index = i*2;
	
		Share_Left_Borders(TData, Coeff_index,shared_buffer, &u_auxiliar, &m_auxiliar, &b_auxiliar);		

		CP_Encoder(	TData, Coeff_index, Bit_Plane, Code_Stream, Code_Stream_pointer, Reserved_codeword, 
				AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, 1, LUT_significance_pointer, LUT_sign_pointer, shared_buffer, LUT, 
				u_auxiliar,	TData[Coeff_index-2], 	TData[Coeff_index-1],   
				m_auxiliar,				TData[Coeff_index+1],
				b_auxiliar,	TData[Coeff_index+2],	TData[Coeff_index+3]	);

		Coeff_index = (i*2)+1;

		Share_Right_Borders(TData, Coeff_index,shared_buffer, &u_auxiliar, &m_auxiliar, &b_auxiliar);

		CP_Encoder(	TData, Coeff_index, Bit_Plane, Code_Stream, Code_Stream_pointer, Reserved_codeword, 
				AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, 2, LUT_significance_pointer, LUT_sign_pointer, shared_buffer, LUT,
				TData[Coeff_index-3],	TData[Coeff_index-2], 	u_auxiliar,   
				TData[Coeff_index-1],				m_auxiliar,
				TData[Coeff_index+1],	TData[Coeff_index+2],	b_auxiliar	);
	}

	Coeff_index = i*2;

	Share_Left_Borders(TData, Coeff_index,shared_buffer, &u_auxiliar, &m_auxiliar, &b_auxiliar);

	CP_Encoder(		TData, Coeff_index, Bit_Plane, Code_Stream, Code_Stream_pointer, Reserved_codeword, 
				AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, 1, LUT_significance_pointer, LUT_sign_pointer, shared_buffer, LUT,
				u_auxiliar,	TData[Coeff_index-2], 	TData[Coeff_index-1],   
				m_auxiliar,				TData[Coeff_index+1],
				0,			0, 			0  	 	);

	Coeff_index = (i*2)+1;

	Share_Right_Borders(TData, Coeff_index,shared_buffer, &u_auxiliar, &m_auxiliar, &b_auxiliar);

	CP_Encoder(		TData, Coeff_index, Bit_Plane, Code_Stream, Code_Stream_pointer, Reserved_codeword, 
				AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, 2, LUT_significance_pointer, LUT_sign_pointer, shared_buffer, LUT,
				TData[Coeff_index-3],	TData[Coeff_index-2], 	u_auxiliar,   
				TData[Coeff_index-1],				m_auxiliar,
				0,			0, 			0      		);

}


inline __device__ void CP_Decoder_launcher(unsigned int* TData, int Mask, VOLATILE REG_DATATYPE* Code_Stream, int Code_Stream_pointer, unsigned int* Codeword, unsigned int* AC_Interval_Lower, unsigned int* AC_Interval_Size, int Shared_pointer, volatile int* Code_Stream_Shared, int LUT_significance_pointer, int LUT_sign_pointer, volatile int* shared_buffer, int const* __restrict__ LUT)
{	

	int i = 0;

	int u_auxiliar=0;
	int m_auxiliar=0;
	int b_auxiliar=0;

	int Coeff_index = i*2;

	Share_Left_Borders(TData, Coeff_index, shared_buffer, &u_auxiliar, &m_auxiliar, &b_auxiliar);

	CP_Decoder(		TData, Coeff_index, Mask, Code_Stream, Code_Stream_pointer, Codeword, 
				AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, 1, LUT_significance_pointer, LUT_sign_pointer, shared_buffer, LUT, 
				0,		0, 			0,   
				m_auxiliar,				TData[Coeff_index+1],
				b_auxiliar,	TData[Coeff_index+2],	TData[Coeff_index+3]	);

	Coeff_index = (i*2)+1;

	Share_Right_Borders(TData, Coeff_index, shared_buffer, &u_auxiliar, &m_auxiliar, &b_auxiliar);

	CP_Decoder(		TData, Coeff_index, Mask, Code_Stream, Code_Stream_pointer, Codeword, 						
				AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, 2, LUT_significance_pointer, LUT_sign_pointer, shared_buffer, LUT,
				0,			0, 			0,      
				TData[Coeff_index-1],				m_auxiliar,
				TData[Coeff_index+1],	TData[Coeff_index+2],	b_auxiliar	);

	for(i=1;i<(CBLOCK_LENGTH -1);i++)
	{

		Coeff_index = i*2;
	
		Share_Left_Borders(TData, Coeff_index,shared_buffer, &u_auxiliar, &m_auxiliar, &b_auxiliar);		

		CP_Decoder(	TData, Coeff_index, Mask, Code_Stream, Code_Stream_pointer, Codeword, 
				AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, 1, LUT_significance_pointer, LUT_sign_pointer, shared_buffer, LUT, 
				u_auxiliar,	TData[Coeff_index-2], 	TData[Coeff_index-1],   
				m_auxiliar,				TData[Coeff_index+1],
				b_auxiliar,	TData[Coeff_index+2],	TData[Coeff_index+3]	);

		Coeff_index = (i*2)+1;

		Share_Right_Borders(TData, Coeff_index,shared_buffer, &u_auxiliar, &m_auxiliar, &b_auxiliar);

		CP_Decoder(	TData, Coeff_index, Mask, Code_Stream, Code_Stream_pointer, Codeword, 
				AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, 2, LUT_significance_pointer, LUT_sign_pointer, shared_buffer, LUT,
				TData[Coeff_index-3],	TData[Coeff_index-2], 	u_auxiliar,   
				TData[Coeff_index-1],				m_auxiliar,
				TData[Coeff_index+1],	TData[Coeff_index+2],	b_auxiliar	);
	}

	Coeff_index = i*2;

	Share_Left_Borders(TData, Coeff_index,shared_buffer, &u_auxiliar, &m_auxiliar, &b_auxiliar);

	CP_Decoder(		TData, Coeff_index, Mask, Code_Stream, Code_Stream_pointer, Codeword, 
				AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, 1, LUT_significance_pointer, LUT_sign_pointer, shared_buffer, LUT,
				u_auxiliar,	TData[Coeff_index-2], 	TData[Coeff_index-1],   
				m_auxiliar,				TData[Coeff_index+1],
				0,			0, 			0  	 	);

	Coeff_index = (i*2)+1;

	Share_Right_Borders(TData, Coeff_index,shared_buffer, &u_auxiliar, &m_auxiliar, &b_auxiliar);

	CP_Decoder(		TData, Coeff_index, Mask, Code_Stream, Code_Stream_pointer, Codeword, 
				AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, 2, LUT_significance_pointer, LUT_sign_pointer, shared_buffer, LUT,
				TData[Coeff_index-3],	TData[Coeff_index-2], 	u_auxiliar,   
				TData[Coeff_index-1],				m_auxiliar,
				0,			0, 			0      		);

}



inline __device__ void MRP_Encoder_launcher(unsigned int* TData, int Bit_Plane, VOLATILE REG_DATATYPE* Code_Stream, int Code_Stream_pointer, unsigned int* Reserved_codeword, unsigned int* AC_Interval_Lower, unsigned int* AC_Interval_Size, int  Shared_pointer, volatile int* Code_Stream_Shared, int Probability){
		
	for(int i=0;i<CBLOCK_LENGTH;i++)
	{

		MRP_Encoder(	TData, i*2, Bit_Plane, Code_Stream, Code_Stream_pointer, Reserved_codeword, 
				AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, Probability);

		MRP_Encoder(	TData, (i*2)+1, Bit_Plane, Code_Stream, Code_Stream_pointer, Reserved_codeword, 
				AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, Probability);
	}
}


inline __device__ void MRP_Decoder_launcher(unsigned int* TData, int Mask, VOLATILE REG_DATATYPE* Code_Stream, int Code_Stream_pointer, unsigned int* Codeword, unsigned int* AC_Interval_Lower, unsigned int* AC_Interval_Size, int  Shared_pointer, volatile int* Code_Stream_Shared, int Probability, int Bit_Plane){
		
	for(int i=0;i<CBLOCK_LENGTH;i++)
	{

		MRP_Decoder(	TData, i*2, Mask, Code_Stream, Code_Stream_pointer, Codeword, 
				AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, Probability, Bit_Plane);

		MRP_Decoder(	TData, (i*2)+1, Mask, Code_Stream, Code_Stream_pointer, Codeword, 
				AC_Interval_Lower, AC_Interval_Size, Shared_pointer, Code_Stream_Shared, Probability, Bit_Plane);
	}
}





inline __device__ void Encode( unsigned int* TData, int MSB, VOLATILE REG_DATATYPE* Code_Stream, int Code_Stream_pointer, int Shared_pointer, volatile int* Code_Stream_Shared, int CB_decomposition_level, int CB_subband, int* shared_buffer, int const* 	__restrict__ LUT)

{	
	
//---Ini
	unsigned int Reserved_codeword = 0;
	unsigned int AC_Interval_Lower = 0;
	unsigned int AC_Interval_Size  = 0;

	int LUT_refinement_pointer = 0;
	int LUT_significance_pointer = 	0;
	int LUT_sign_pointer =  0;

	int Bit_Plane = MSB;

	Initialize_LUT_pointers(&LUT_refinement_pointer, &LUT_significance_pointer , &LUT_sign_pointer, CB_decomposition_level, CB_subband, MSB);
	

#if 	CPASS_3 == 1		

	int LUT_pointer_aux = 	((LUT_N_SUBBANDS) * (LUT_N_BITPLANES) * (N_CONTEXT_SIGNIFICANCE) * (DWT_LEVELS)) 	+ (LUT_N_BITPLANES*N_CONTEXT_SIGNIFICANCE) +
						((LUT_N_SUBBANDS) * (LUT_N_BITPLANES) * (N_CONTEXT_SIGN) * (DWT_LEVELS)) 			+ (LUT_N_BITPLANES*N_CONTEXT_SIGN);

	CP_Encoder_launcher(	TData, Bit_Plane, Code_Stream, Code_Stream_pointer, &Reserved_codeword, &AC_Interval_Lower, &AC_Interval_Size, 
				Shared_pointer, Code_Stream_Shared, LUT_significance_pointer + LUT_pointer_aux, LUT_sign_pointer + LUT_pointer_aux, shared_buffer, LUT);
	Bit_Plane--;

	Update_LUT_pointers(&LUT_refinement_pointer, &LUT_significance_pointer , &LUT_sign_pointer);
		
#endif		

	for(; Bit_Plane >= 0; Bit_Plane--)
	//for(; Bit_Plane >= MSB; Bit_Plane--)
	{
		
		SPP_Encoder_launcher(	TData, Bit_Plane, Code_Stream, Code_Stream_pointer, &Reserved_codeword, &AC_Interval_Lower, &AC_Interval_Size, 
					Shared_pointer, Code_Stream_Shared, LUT_significance_pointer, LUT_sign_pointer, shared_buffer, LUT);	

		MRP_Encoder_launcher(	TData, Bit_Plane, Code_Stream, Code_Stream_pointer, &Reserved_codeword, &AC_Interval_Lower, &AC_Interval_Size, 
				Shared_pointer, Code_Stream_Shared, LUT[LUT_refinement_pointer]);

		
#if 	CPASS_3 == 1

		CP_Encoder_launcher(	TData, Bit_Plane, Code_Stream, Code_Stream_pointer, &Reserved_codeword, &AC_Interval_Lower, &AC_Interval_Size, 
					Shared_pointer, Code_Stream_Shared, LUT_significance_pointer + LUT_pointer_aux, LUT_sign_pointer + LUT_pointer_aux, shared_buffer, LUT);
#endif		 

		Update_LUT_pointers(&LUT_refinement_pointer, &LUT_significance_pointer , &LUT_sign_pointer);	

	}
	//Flush not exhausted CW then the encoding process finish	
	Code_Stream[Reserved_codeword] = AC_Interval_Lower;
}


inline __device__ void Decode( unsigned int* TData, int MSB, VOLATILE REG_DATATYPE* Code_Stream, int Code_Stream_pointer, int Shared_pointer, volatile int* Code_Stream_Shared, int CB_decomposition_level, int CB_subband, int* shared_buffer, int const* 	__restrict__ LUT)

{	
	
//---Ini
	unsigned int Codeword = 0;
	unsigned int AC_Interval_Lower = 0;
	unsigned int AC_Interval_Size  = 0;

	int LUT_refinement_pointer = 0;
	int LUT_significance_pointer = 	0;
	int LUT_sign_pointer =  0;

	int Bit_Plane = MSB;
	int Mask = 0x3 << (Bit_Plane);

	if (Bit_Plane==0) Mask&=0x2;

	Initialize_LUT_pointers(&LUT_refinement_pointer, &LUT_significance_pointer , &LUT_sign_pointer, CB_decomposition_level, CB_subband, MSB);


#if 	CPASS_3 == 1		

	int LUT_pointer_aux = 	((LUT_N_SUBBANDS) * (LUT_N_BITPLANES) * (N_CONTEXT_SIGNIFICANCE) * (DWT_LEVELS)) 	+ (LUT_N_BITPLANES*N_CONTEXT_SIGNIFICANCE) +
						((LUT_N_SUBBANDS) * (LUT_N_BITPLANES) * (N_CONTEXT_SIGN) * (DWT_LEVELS)) 			+ (LUT_N_BITPLANES*N_CONTEXT_SIGN);

	CP_Decoder_launcher(	TData, Mask, Code_Stream, Code_Stream_pointer, &Codeword, &AC_Interval_Lower, &AC_Interval_Size, 
				Shared_pointer, Code_Stream_Shared, LUT_significance_pointer + LUT_pointer_aux, LUT_sign_pointer + LUT_pointer_aux, shared_buffer, LUT);
	Bit_Plane--;

	Mask>>=1;

	if (Bit_Plane==0) Mask=0x2;

	Update_LUT_pointers(&LUT_refinement_pointer, &LUT_significance_pointer , &LUT_sign_pointer);

		
#endif		



	for(; Bit_Plane >= 0; Bit_Plane--)
	//for(; Bit_Plane >= MSB; Bit_Plane--)	
	{
		
		SPP_Decoder_launcher(	TData, Mask, Code_Stream, Code_Stream_pointer, &Codeword, &AC_Interval_Lower, &AC_Interval_Size, 
					Shared_pointer, Code_Stream_Shared, LUT_significance_pointer, LUT_sign_pointer, shared_buffer, LUT);	

		MRP_Decoder_launcher(	TData, Mask, Code_Stream, Code_Stream_pointer, &Codeword, &AC_Interval_Lower, &AC_Interval_Size, 
				Shared_pointer, Code_Stream_Shared, LUT[LUT_refinement_pointer], Bit_Plane);

#if 	CPASS_3 == 1

		CP_Decoder_launcher(	TData, Mask, Code_Stream, Code_Stream_pointer, &Codeword, &AC_Interval_Lower, &AC_Interval_Size, 
					Shared_pointer, Code_Stream_Shared, LUT_significance_pointer + LUT_pointer_aux, LUT_sign_pointer + LUT_pointer_aux, shared_buffer, LUT);
#endif		 

		Mask>>=1;

		if (Bit_Plane==1) Mask=0x2;

		Update_LUT_pointers(&LUT_refinement_pointer, &LUT_significance_pointer , &LUT_sign_pointer);


	}
	
}


//CUDA KERNELS

//ENCODER KERNEL

__global__ void Kernel_BPC_CODER(	

									int*	Input,
									int*	Output,
									int	Input_XSize,
									int	Input_YSize,
									int 	NWarps_Block,
									int 	NWarps_X,
									int const* __restrict__	LUT

								)
{

#if 	FERMI == 1
	__shared__ REG_DATATYPE shared_buffer[THREADBLOCK_SIZE+2];
#else 
	__shared__ REG_DATATYPE shared_buffer[1];
#endif
		
	const int Shared_amount = THREADBLOCK_SIZE/WARPSIZE;

	int Shared_per_warp = Shared_amount/NWarps_Block;

	register unsigned int TData[ CBLOCK_LENGTH * NELEMENTS_THREAD_X ];
	__shared__ volatile REG_DATATYPE Code_Stream_Shared[Shared_amount];
	

	int LaneID = 			threadIdx.x & 0x1f;
	int WarpID = 			(((threadIdx.x >> 5) + (blockIdx.x * NWarps_Block)));
	int TCoordinate_X = 		(((WarpID % NWarps_X) * CBLOCK_WIDTH) + (LaneID * NELEMENTS_THREAD_X));
	int TCoordinate_Y = 		((WarpID/NWarps_X) * CBLOCK_LENGTH);
	int TCoordinate =		(Input_XSize*TCoordinate_Y) + TCoordinate_X;

	int WCoordinate = WarpID * CBLOCK_LENGTH * CBLOCK_WIDTH;

	int MSB = 0;
	//CodeBlock decomposition level: from 0 to image decomposition levels
	int CB_decomposition_level = -1;
	//CodeBlock Subband: LL = 0, HL = 0, LH = 1, HH = 2
	int CB_subband = -1;
		
	
	if(WarpID<(Input_XSize/CBLOCK_WIDTH)*(Input_YSize/CBLOCK_LENGTH)){

		Find_Subband(&CB_decomposition_level, &CB_subband, TCoordinate_X, TCoordinate_Y, Input_XSize, Input_YSize);
	
		Read_Coefficients(Input, (int*)TData, TCoordinate, Input_XSize);

		Find_MSB(TData, &MSB, shared_buffer);

		Code_Stream_Shared[(WarpID % NWarps_Block)*Shared_per_warp]=0;
	
	//We write the codestream with the following layout:
	// 	CB1 - CB2	-->	CS1 - CS1
	// 	CB1 - CB2	-->	CS2 - CS2
	// 	CB3 - CB4	-->	CS3 - CS3
	// 	CB3 - CB4	-->	CS4 - CS4

		Output[WCoordinate] = MSB;

		Encode(TData, MSB, Output, WCoordinate+1, (WarpID % NWarps_Block), Code_Stream_Shared, CB_decomposition_level, CB_subband, shared_buffer, LUT);

	}

}

//DECODER KERNEL
__global__ void Kernel_BPC_DECODER(	

									int*	Input_CodeStream,
									int*	Output_Image,
									int	Input_XSize,
									int	Input_YSize,
									int 	NWarps_Block,
									int 	NWarps_X,
									int const* __restrict__	LUT

								)
{

	
#if 	FERMI == 1
	__shared__ REG_DATATYPE shared_buffer[THREADBLOCK_SIZE+2];
#else 
	__shared__ REG_DATATYPE shared_buffer[1];
#endif
		
	const int Shared_amount = THREADBLOCK_SIZE/WARPSIZE;

	int Shared_per_warp = Shared_amount/NWarps_Block;

	register unsigned int TData[ CBLOCK_LENGTH * NELEMENTS_THREAD_X ];

	__shared__ volatile REG_DATATYPE Code_Stream_Shared[Shared_amount];
	

	int LaneID = 			threadIdx.x & 0x1f;
	int WarpID = 			(((threadIdx.x >> 5) + (blockIdx.x * NWarps_Block)));
	int TCoordinate_X = 		(((WarpID % NWarps_X) * CBLOCK_WIDTH) + (LaneID * NELEMENTS_THREAD_X));
	int TCoordinate_Y = 		((WarpID/NWarps_X) * CBLOCK_LENGTH);
	int TCoordinate =		(Input_XSize*TCoordinate_Y) + TCoordinate_X;

	int WCoordinate = WarpID * CBLOCK_LENGTH * CBLOCK_WIDTH;

	int MSB = 0;
	//CodeBlock decomposition level: from 0 to image decomposition levels
	int CB_decomposition_level = -1;
	//CodeBlock Subband: LL = 0, HL = 0, LH = 1, HH = 2
	int CB_subband = -1;
	
	
	if(WarpID<(Input_XSize/CBLOCK_WIDTH)*(Input_YSize/CBLOCK_LENGTH)){

		Find_Subband(&CB_decomposition_level, &CB_subband, TCoordinate_X, TCoordinate_Y, Input_XSize, Input_YSize);
	
		Initialize_Coefficients((int*)TData);

		Code_Stream_Shared[(WarpID % NWarps_Block)*Shared_per_warp]=0;	
	
		MSB = Input_CodeStream[WCoordinate];

		Decode(TData, MSB, Input_CodeStream, WCoordinate+1, (WarpID % NWarps_Block), Code_Stream_Shared, CB_decomposition_level, CB_subband, shared_buffer, LUT);

		Write_Coefficients(Output_Image, TData, TCoordinate, Input_XSize);

	}

}


// -----------------------------------------------------------------------
//HOST FUNCTIONS ---------------------------------------------------------
// -----------------------------------------------------------------------


static inline void Kernel_launcher(int Direction, int DSize_X, int DSize_Y, DATATYPE* DData_Initial, DATATYPE* DData_Final, int* LUT){

	int Warps_Row = 				(int)ceil(DSize_X/(float)CBLOCK_WIDTH);	
	int Warps_Column = 				(int)ceil(DSize_Y/(float)CBLOCK_LENGTH);
	
	int CUDA_number_warps = 			Warps_Row * Warps_Column;					
	int CUDA_number_blocks = 			(int)ceil((CUDA_number_warps*WARPSIZE)/(float)(THREADBLOCK_SIZE));

	switch(Direction){

		case CODE:			

			Kernel_BPC_CODER<<<CUDA_number_blocks,THREADBLOCK_SIZE, SYNTHETIC_SHARED>>>
									(	

										DData_Initial,
										DData_Final,
										DSize_X,
										DSize_Y,
										THREADBLOCK_SIZE/WARPSIZE,
										Warps_Row,
										LUT
									);	
		break;

		case DECODE:

			Kernel_BPC_DECODER<<<CUDA_number_blocks,THREADBLOCK_SIZE, SYNTHETIC_SHARED>>>
									(	

										DData_Initial,
										DData_Final,
										DSize_X,
										DSize_Y,
										THREADBLOCK_SIZE/WARPSIZE,
										Warps_Row,
										LUT

									);	
		break;

	}	

}

static inline void Device_memory_allocator(DATATYPE* HData, int HDSize_X, int HDSize_Y, DATATYPE** DData_Initial, DATATYPE** DData_Final){
	
	
	size_t 		DSize  = HDSize_X*HDSize_Y*sizeof(DATATYPE);
	
	HANDLE_ERROR(cudaMalloc ((void**) &(*DData_Initial), DSize));
	HANDLE_ERROR(cudaMalloc ((void**) &(*DData_Final), DSize));
	HANDLE_ERROR(cudaMemset (*DData_Final, -1, DSize));

	HANDLE_ERROR(cudaMemcpy(*DData_Initial, HData, DSize, cudaMemcpyHostToDevice));

}


static inline void Device_memory_deallocator(DATATYPE* HData, int HDSize_X, int HDSize_Y, DATATYPE* DData_Initial, DATATYPE* DData_Final, int* LUT){
	
	size_t 		DSize = HDSize_X*HDSize_Y*sizeof(DATATYPE);

	HANDLE_ERROR(cudaMemcpy(HData, DData_Final, DSize, cudaMemcpyDeviceToHost));

	HANDLE_ERROR(cudaFree(DData_Initial));
	HANDLE_ERROR(cudaFree(DData_Final));

	HANDLE_ERROR(cudaFree(LUT));

}


void CUDA_BPC(int Direction, DATATYPE* HData, int HDSize_X, int HDSize_Y, int* LUT)
{
			
	DATATYPE* 		DData_Initial; 
	DATATYPE* 		DData_Final;
	
	
#if PREFSHARED == 1

	HANDLE_ERROR(cudaDeviceSetCacheConfig(cudaFuncCachePreferShared));
#else

	HANDLE_ERROR(cudaDeviceSetCacheConfig(cudaFuncCachePreferL1));

#endif

	Device_memory_allocator(HData, HDSize_X, HDSize_Y, &DData_Initial, &DData_Final);
	Kernel_launcher(Direction, HDSize_X, HDSize_Y, DData_Initial, DData_Final, LUT);
	Device_memory_deallocator(HData, HDSize_X, HDSize_Y, DData_Initial, DData_Final, LUT);

}

	
	
