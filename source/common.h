#include <stdio.h>

#define HANDLE_ERROR(error) (HandleError(error, __FILE__, __LINE__ ))	

void HandleError( cudaError_t , const char * ,  int);


typedef enum BPC{
  	CODE,
	DECODE,
}BPC;


#if !defined(FERMI)
	#define FERMI	0
#endif

#if !defined(DWT_LEVELS)
	#define DWT_LEVELS				5
#endif


#if !defined(CPASS_3)
	#define CPASS_3				0
#endif
 

#define LUT_N_BITPLANES			15
#define LUT_N_SUBBANDS			3

#define N_CONTEXT_REFINEMENT		1
#define N_CONTEXT_SIGN			4
#define N_CONTEXT_SIGNIFICANCE		9



#define VOLATILE



#define WARPSIZE 				32

#define NELEMENTS_THREAD_X 			2	

#if !defined(CBLOCK_LENGTH)
	#define CBLOCK_LENGTH 			64
#endif

#if !defined(CBLOCK_WIDTH)
	#define CBLOCK_WIDTH	(NELEMENTS_THREAD_X*WARPSIZE)
#endif

#define CODEWORD_SIZE				16

#define MULT_PRECISION				7

#define SHARED_AMOUNT_RATIO			2




#if !defined(THREADBLOCK_SIZE)
	#define THREADBLOCK_SIZE		128
#endif



//-------------------------------------------------------

#if !defined(NEXPERIMNETS)
	#define NEXPERIMNETS		5
#endif

#if !defined(EXPERIMENT_INI)
	#define EXPERIMENT_INI		2048
#endif

#if !defined(EXPERIMENT_INCREMENT)
	#define EXPERIMENT_INCREMENT		2048
#endif


//-------------------------------------------------------

#if !defined(PREFSHARED)
	#define PREFSHARED				0
#endif

#if !defined(SYNTHETIC_SHARED)
	#define SYNTHETIC_SHARED		0
#endif

//-------------------------------------------------------

#if !defined(NEXPERIMNETS)
	#define NEXPERIMNETS			10
#endif

#if !defined(EXPERIMENT_INI)
	#define EXPERIMENT_INI			1024
#endif

#if !defined(EXPERIMENT_INCREMENT)
	#define EXPERIMENT_INCREMENT	1024
#endif

//-------------------------------------------------------


#if !defined(DATATYPE_16BITS_or_32BITS)
	#define DATATYPE_16BITS_or_32BITS 	1
#endif

//-------------------------------------------------------

#if !defined(READ)
	#define READ 	1
#endif

#if !defined(WRITE)
	#define WRITE 	1
#endif


//-------------------------------------------------------


#define REG_DATATYPE			int


#if DATATYPE_16BITS_or_32BITS == 1

	#define DATATYPE			int
	#define DATATYPE2 			int
	
#else

	#define DATATYPE 			int
	#define DATATYPE2 			int2

#endif


//-------------------------------------------------------
