# CUDA_BPC-PaCo

The code of this project has been developed for research pourposes. 
The "main_example.cu" contains a sample code that takes an image and a set of LUTs with probabilities, and performs a call to the BPC-PaCo encoder and decoder. 
"BPC-PaCo.cu" contains the CUDA code of this research with all the encoding and decoding functions. 
This current version is a proof of concept that takes images of size i\*2 x i\*2 with i > 1024. The implementation may not properly work for other image sizes.

The sample image used is subdivided in codeblocks and has to be unzipped before running the example.

The sample main can be compiled, executed and profiled using the nvcc compiler as follows:

	$nvcc -c -O3 -arch=sm_52 -use_fast_math -Xptxas -v -maxrregcount=32 -DCPASS_3=0 source/*.cu
	$nvcc *.o -o test
	$nvprof --print-gpu-trace ./test

Depending on your device architecture you may have to use -arch=sm_20, -arch=sm_35, etc.
	
There are some precompilation parameters defined in "common.h". Most important are:

	-DCPASS_3: =0 BPC-PaCo with 2 passes, =1 BPC-PaCo with 3 passes
	
	-DFERMI: =1 to run on Fermi or older archictures, default at 0. 
	
	-DDWT_LEVELS: Wavelet levels of the input data, default is 5.
	
	-DCBLOCK_LENGTH: Codeblock length, default is 64.
