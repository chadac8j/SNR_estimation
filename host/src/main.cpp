/******************************************************************************
*  @file    main.cpp
*  @author  Chad Cole
*  @date    6/3/2021
*  @version 1.0
*
*  @brief Host side OpenCL implementation of an SNR estimator
*
*  @section DESCRIPTION
*
*  This program sets up the host side interface for a FPGA implemetation
*  of a SNR estimator used in the DVB-S2 and DVB-S2X waveform.
*
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <iomanip>
#include <map>
#include <random>
#include <cmath>
#include "CL/opencl.h"
#include "AOCLUtils/aocl_utils.h"
#include <malloc.h>


using namespace aocl_utils;
using namespace std;

#define AOCL_ALIGNMENT  64
//#define VERBOSE


#define NORMAL_MODE     0
#define TEST_MODE       1

#define HARDWARE_PLAT   0
#define EMULATION_PLAT  1


typedef char msg_t;
typedef char msg2_t;


#define SLOT_LEN         4096 // DVB-S2 slot length is 90 symbols
//#define SLOT_LEN         1100 // DVB-S2 slot length is 90 symbols
//#define SLOT_LEN         90 // DVB-S2 slot length is 90 symbols


enum KERNELS {
K_READER,
K_SNR_EST_LUT_CORRECTION,
K_WRITER,
K_NUM_KERNELS
};

static const char* kernel_names[K_NUM_KERNELS] =
{
"data_in",
"snr_est_LUT_correction",
"data_out"
};


const char *device_kernel;

// OpenCL runtime configuration
cl_platform_id platform = NULL;
unsigned num_devices = 0;
cl_device_id device;
cl_context context = NULL;
cl_command_queue queue[K_NUM_KERNELS];
cl_program program = NULL;
cl_kernel kernel[K_NUM_KERNELS];
cl_mem input_noisy_message_I_buf;
cl_mem input_noisy_message_Q_buf;
cl_mem output_buf;

unsigned int slotLen = SLOT_LEN;
unsigned int numFrames = 0;
unsigned int SNR_in = 0;
unsigned int SNR_expected = 0;
unsigned int SOF_ind = 0;
int input_file_size = 0;
int num_output_frames = 1;
//int num_output_frames = 4;

void *noisyDataIn_I_array_ptr = alignedMalloc((slotLen)*sizeof(char));
void *noisyDataIn_Q_array_ptr = alignedMalloc((slotLen)*sizeof(char));
void *dout_snr_est_ptr = alignedMalloc(num_output_frames*sizeof(short));


char *noisyDataIn_I = (char *)noisyDataIn_I_array_ptr;
char *noisyDataIn_Q = (char *)noisyDataIn_Q_array_ptr;
// I pass the numerator and denominator out as unsigned longs
short *dout_snr_est = (short *)dout_snr_est_ptr;


// Function prototypes
int read_test_vector_file_char(const char *filename, char *din_array);
int read_test_vector_file_short(const char *filename, short *din_array);
bool init_opencl();
void run();
void cleanup();
int verify_output();


//test vector data files

const char *input_noisy_sym_file_I;
const char *input_noisy_sym_file_Q;
const char *output_data_file;


char ptype = EMULATION_PLAT;


//****************************************
// Option Parser structures and routines
//****************************************
struct Arg : public option::Arg
{
	static void printError(const char* msg1, const option::Option& opt, const char* msg2)
	{
		fprintf(stderr, "%s", msg1);
		fwrite(opt.name, opt.namelen, 1, stderr);
		fprintf(stderr, "%s", msg2);
	}

	static option::ArgStatus Unknown(const option::Option& option, bool msg)
	{
		if (msg) printError("Unknown option '", option, "'\n");
		return option::ARG_ILLEGAL;
	}

	static option::ArgStatus Required(const option::Option& option, bool msg)
	{
		if (option.arg != 0)
			return option::ARG_OK;

		if (msg) printError("Option '", option, "' requires an argument\n");
		return option::ARG_ILLEGAL;
	}

	static option::ArgStatus NonEmpty(const option::Option& option, bool msg)
	{
		if (option.arg != 0 && option.arg[0] != 0)
			return option::ARG_OK;

		if (msg) printError("Option '", option, "' requires a non-empty argument\n");
		return option::ARG_ILLEGAL;
	}

	static option::ArgStatus Numeric(const option::Option& option, bool msg)
	{
		char* endptr = 0;
		if (option.arg != 0 && strtol(option.arg, &endptr, 10)){};
		if (endptr != option.arg && *endptr == 0)
			return option::ARG_OK;

		if (msg) printError("Option '", option, "' requires a numeric argument\n");
		return option::ARG_ILLEGAL;
	}
};

enum  optionIndex { UNKNOWN, HELP, NFRAME, EMODE, HMODE, N_FRAMES, SNR };
const option::Descriptor usage[] = {
	{ UNKNOWN, 0, "", "", Arg::Unknown, "USAGE: example_arg [options]\n\n"
	"Options:" },
	{ HELP, 0, "", "help", Arg::None, "  \t--help  \tPrint usage and exit." },
	{ N_FRAMES, 0, "f", "number of frames", Arg::Required, "  -f <arg>, \t--required=<arg>  \tNumber of DVB-S2 frames to decode." },
	{ EMODE, 0, "e", "run in emulation mode", Arg::None, "  -e\t\tRun in emulation mode" },
	{ HMODE, 0, "h", "run on hardware", Arg::None, "  -h\t\tRun on hardware" },
	{ SNR, 0, "s", "SNR test input val ", Arg::Required, "  -s <arg>, \t--required=<arg>  \tSNR_in\n \t\t0 = 3 dB\n, \t\t1 = 6 dB,\n \t\t2 = 9 dB,\n"
	"\t\t3 = 12 dB\n, \t\t4 = NA,\n \t\t5 = NA." },
	{ 0, 0, 0, 0, 0, 0 } };


/*************************************************************************

@brief The main function serves as the launch point for initializing
an OpenCL application.  The function finds a device, sets up
the context and read/write buffers.  It then launches the
kernel and verifes the results against Matlab generated
test vectors.

@param argc argument count
@param argv argument varaibles
@return int if a negative value is returned the function failed

**************************************************************************/
int main(int argc, char **argv) {

	
	argc -= (argc>0); argv += (argc>0); // skip program name argv[0] if present
	option::Stats stats(usage, argc, argv);

#ifdef __GNUC__
	// GCC supports C99 VLAs for C++ with proper constructor calls.
	option::Option options[stats.options_max], buffer[stats.buffer_max];
#else
	// use calloc() to allocate 0-initialized memory. It's not the same
	// as properly constructed elements, but good enough. Obviously in an
	// ordinary C++ program you'd use new[], but this file demonstrates that
	// TLMC++OP can be used without any dependency on the C++ standard library.
	option::Option* options = (option::Option*)calloc(stats.options_max, sizeof(option::Option));
	option::Option* buffer = (option::Option*)calloc(stats.buffer_max, sizeof(option::Option));
#endif

	option::Parser parse(usage, argc, argv, options, buffer);

	if (parse.error())
		return 1;

	if (options[HELP] || argc == 0)
	{
		int columns = getenv("COLUMNS") ? atoi(getenv("COLUMNS")) : 80;
		option::printUsage(fwrite, stdout, usage, columns);
		return 0;
	}

	for (int i = 0; i < parse.optionsCount(); ++i)
	{
		option::Option& opt = buffer[i];
		switch (opt.index())
		{
		case HELP:
		case N_FRAMES:
			numFrames = atoi(opt.arg);
			break;
		case SNR:
			SNR_in = atoi(opt.arg);
			break;
		case EMODE:
			ptype = EMULATION_PLAT;
			break;
		case HMODE:
			ptype = HARDWARE_PLAT;
			break;
		case UNKNOWN:
			// not possible because Arg::Unknown returns ARG_ILLEGAL
			// which aborts the parse with an error
			break;
		}
	}

	for (int i = 0; i < parse.nonOptionsCount(); ++i)
		fprintf(stdout, "Non-option argument #%d is %s\n", i, parse.nonOption(i));

// These are I/Q test input files at various SNR's and # of samples	
	if (SNR_in == 0) {
//	input_noisy_sym_file_I = "noisy_sym_IN_I_3dB.txt";
//	input_noisy_sym_file_Q = "noisy_sym_IN_Q_3dB.txt";
	  input_noisy_sym_file_I = "noisy_sym_IN_I_3dB_1100Samp.txt";
	  input_noisy_sym_file_Q = "noisy_sym_IN_Q_3dB_1100Samp.txt"; 
	  SNR_expected = 3; }
	if (SNR_in == 1) {
//	input_noisy_sym_file_I = "noisy_sym_IN_I_6dB.txt";
//	input_noisy_sym_file_Q = "noisy_sym_IN_Q_6dB.txt"; 
	  input_noisy_sym_file_I = "noisy_sym_IN_I_6dB_1100Samp.txt";
      input_noisy_sym_file_Q = "noisy_sym_IN_Q_6dB_1100Samp.txt"; 
      SNR_expected = 6; }
	if (SNR_in == 2) {
//	input_noisy_sym_file_I = "noisy_sym_IN_I_9dB.txt";
//	input_noisy_sym_file_Q = "noisy_sym_IN_Q_9dB.txt";
	input_noisy_sym_file_I = "noisy_sym_IN_I_9dB_1100Samp.txt";
	input_noisy_sym_file_Q = "noisy_sym_IN_Q_9dB_1100Samp.txt"; 
	SNR_expected = 9; }
	if (SNR_in == 3) {
	input_noisy_sym_file_I = "noisy_sym_IN_I_highSNR_freqOffset_1100Samp_pilots.txt";
	input_noisy_sym_file_Q = "noisy_sym_IN_Q_highSNR_freqOffset_1100Samp_pilots.txt"; 
	SNR_expected = 9; }
	if (SNR_in == 4) {
	input_noisy_sym_file_I = "noisy_sym_IN_I_minus3dB_1100Samp_21Mag.txt";
	input_noisy_sym_file_Q = "noisy_sym_IN_Q_minus3dB_1100Samp_21Mag.txt"; 
	SNR_expected = -3; }
	if (SNR_in == 5) {
	input_noisy_sym_file_I = "noisy_sym_IN_I_minus6dB_1100Samp_21Mag.txt";
	input_noisy_sym_file_Q = "noisy_sym_IN_Q_minus6dB_1100Samp_21Mag.txt"; 
	SNR_expected = -6; }
	if (SNR_in == 6) {
	input_noisy_sym_file_I = "noisy_sym_IN_I_minus10dB_1100Samp_21Mag.txt";
	input_noisy_sym_file_Q = "noisy_sym_IN_Q_minus10dB_1100Samp_21Mag.txt"; 
	SNR_expected = -10; }
	if (SNR_in == 7) {
	input_noisy_sym_file_I = "noisy_sym_IN_I_highSNR_freqOffset_4096Samp_pilots.txt";
	input_noisy_sym_file_Q = "noisy_sym_IN_Q_highSNR_freqOffset_4096Samp_pilots.txt"; 
	SNR_expected = 30; }
	
	output_data_file = "snr_est_OUT.txt";

	if (ptype == EMULATION_PLAT)
		device_kernel = "SNR_estimator_LUT_correction_top";
//		device_kernel = "snr_estimator_em";
	else
		device_kernel = "snr_estimator";

#ifdef VERBOSE
	//************************************
	// Display Frame Type and Code Rate
	//************************************
	printf("\n\n****************************************************************************\n");
	if (ptype == HARDWARE_PLAT) printf("Platform=Hardware\n");
	else printf("Platform=Emulator\n");
	printf("Number of symbols per slot=%d\n", slotLen);
	printf("****************************************************************************\n\n");
#endif
		

	//**********************************************
	//get input and output test vectors from files
	//**********************************************
/*	if (read_test_vector_file_short(input_data_file, msgIn) < 0)
	{
		printf("Error opening input data vector file\n");
		return -1;
	} */
	// read I and Q in same file or separate?
	if (read_test_vector_file_char(input_noisy_sym_file_I, noisyDataIn_I) < 0)
	{
		printf("Error opening input noisy data I vector file\n");
		return -1;
	}
	if (read_test_vector_file_char(input_noisy_sym_file_Q, noisyDataIn_Q) < 0)
	{
		printf("Error opening input noisy data Q vector file\n");
		return -1;
	}

	num_output_frames = input_file_size/1024;
	
	//_***********************
	// Initialize OpenCL.
	//_***********************
	if (!init_opencl()) {
		return -1;
	}
	
	//_************
	// Run Kernel
	//_************
	run();

	//_******************
	// Verify Results
	//_******************
	if (verify_output() == 0)
	{
		printf("Estimated SNR is within +/-1 of real value.... PASSED!\n");
	}
	else
	{
		printf("Estimated SNR not within +/-1 of real value.... FAILED!\n");
	}
	
	//_*********************************
	// Free the resources allocated
	//_*********************************
	//cleanup();

	return 0;
}






/************************************************************************

@brief The read_test_vector_file function reads a text file of input
parameters intended to be used as inputs to an OpenCL module

@param filename is a const char pointer containg the name of the
file to be parsed
@param din_array is a float array that holds the fields parsed data
@return int if a negative value is returned the function failed

*************************************************************************/


int read_test_vector_file_char(const char *filename, char *din_array)
{
	FILE* file = fopen(filename, "rt");
	if (file == NULL) {
		printf("File %s could not be opened\n", filename);
		return -1;
	}
	char line[256];
	input_file_size = 0;  // reset for each input such as I and Q
	while (fgets(line, sizeof(line), file))
	{
		int r = strtol(line, nullptr, 10);
		char tmp = (char)r;
		*din_array++ = (char)r;
		input_file_size += 1;
	}

	fclose(file);
	return 0;
}


int read_test_vector_file_short(const char *filename, short *din_array)
{
	FILE* file = fopen(filename, "rt");
	char line[256];
	const char *s = ", ";
	char *token;
	char* pEnd;
	float tmp;
	char ctmp;
	int cnt = 0;

	while (fgets(line, sizeof(line), file))
	{
		int r = strtol(line, nullptr, 10);
		short tmp = (short)r;
		*din_array++ = (short)r;

	}

	fclose(file);
	return 0;
}



/**************************************************************

@brief The verify_output function verifies the kernels outputs
against a text file containg outptus generated in Matlab

@return int if less than 0 an error has occurred

**************************************************************/
int verify_output()
{

//    printf("in verify_output:  The estimated SNR numerator is =%lu \n", dout_snr_est[0]);
//	printf("in verify_output:  The estimated SNR denominator is =%lu \n", dout_snr_est[1]);
	int bool_val = 1;
	printf("in verify_output: The input_file_size is =%d\n", input_file_size);
	for (int i = 0; i < num_output_frames; i++) {
		printf("in verify_output: The estimated SNR is =%f\n", (double)dout_snr_est[i]/10);
//	if (abs(((double)dout_snr_est[0]/(double)dout_snr_est[1]) - SNR_expected) < 1)
		if (abs(((double)dout_snr_est[i]/10) - SNR_expected) > 1)
			bool_val = 0; //failed
	}
	return bool_val;
}

/*************************************************************************

@brief The init_opencl function intializes the OpenCL objects.
Essentially it looks for an Altera OpenCl device, creates a
context for that device.  The AOC compiled kernel is pointed to
using the global variable "device_kernel".  Provided a valid .aocx
compiled kernel is found, the FPGA is programmed and command and data
queues are generated.

@return bool True if successful, otherwise false

**************************************************************************/
bool init_opencl() {
	cl_int status;

	//printf("Initializing OpenCL\n");

	if (!setCwdToExeDir()) {
		return false;
	}


	//****************************
	// Get the OpenCL platform.
	//****************************
/*	platform = findPlatform("Intel(R) FPGA SDK for OpenCL(TM)");
	
	if (platform == NULL) {
		printf("ERROR: Unable to find Intel(R) FPGA OpenCL platform.\n");
		return false;
	} */
    if (ptype == EMULATION_PLAT) {
             // new 'fast' emulator
        platform = findPlatform("Intel(R) FPGA Emulation Platform for OpenCL(TM)");
        if (platform == NULL) {
            printf("ERROR: Unable to find Intel(R) FPGA Emulation Platform for OpenCL(TM).\n");
         // For legacy emulator
            platform = findPlatform("Intel(R) FPGA SDK for OpenCL(TM)");
            if (platform == NULL) {
                printf("ERROR: Unable to find Intel(R) FPGA Legacy Emulation Platform for OpenCL(TM).\n");
                return false;
            }
        }
    }else{    //use for hardware and simulation
        platform = findPlatform("Intel(R) FPGA SDK for OpenCL(TM)");
        if (platform == NULL) {
            printf("ERROR: Unable to find Intel(R) FPGA SDK for OpenCL(TM).\n");
            return false;
        }
    }



	//*************************************
	// Query the available OpenCL devices.
	//*************************************
	scoped_array<cl_device_id> devices;
	devices.reset(getDevices(platform, CL_DEVICE_TYPE_ALL, &num_devices));
	device = devices[0];
#ifdef VERBOSE
	printf("Platform: %s\n", getPlatformName(platform).c_str());
	printf("Using %d device(s)\n", num_devices);
	for (unsigned i = 0; i < num_devices; ++i) {
		printf("  %s\n", getDeviceName(device).c_str());
	}
#endif

	//*********************
	// Create the context.
	//*********************
	context = clCreateContext(NULL, num_devices, &device, &oclContextCallback, NULL, &status);
	checkError(status, "Failed to create context");

	//************************************
	// Create the program for all device. 
	//************************************
	std::string binary_file = getBoardBinaryFile(device_kernel, device);
#ifdef VERBOSE
	printf("Using AOCX: %s\n", binary_file.c_str());
#endif
	program = createProgramFromBinary(context, binary_file.c_str(), &device, num_devices);

	//*******************************************
	// Build the program that was just created.
	//*******************************************
	status = clBuildProgram(program, 0, NULL, "", NULL, NULL);
	checkError(status, "Failed to build program");

	//**********************
	// Create Command queue.
	//**********************
	for (int i = 0; i < K_NUM_KERNELS; ++i)
	{
		queue[i] = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &status);
		checkError(status, "Failed to create command queue %d", i);
	}

	//*****************
	// Create Kernel.
	//*****************
	for (int i = 0; i < K_NUM_KERNELS; ++i)
	{
		kernel[i] = clCreateKernel(program, kernel_names[i], &status);
		checkError(status, "Failed to create kernel %s", kernel_names[i]);
	}

	//**********************
	// Create Input buffers.
	//**********************
	//message buffer - holds encoded message 
	// CAC - Need one for each of two complex input messages?
	input_noisy_message_I_buf = clCreateBuffer(context, CL_MEM_READ_ONLY,
		(SLOT_LEN) * sizeof(char), NULL, &status);
	input_noisy_message_Q_buf = clCreateBuffer(context, CL_MEM_READ_ONLY,
		(SLOT_LEN) * sizeof(char), NULL, &status);
	
	//**********************  
	// Create Output buffer.
	//**********************
	output_buf = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
		num_output_frames*sizeof(short), NULL, &status);
	checkError(status, "Failed to create buffer for output");


	return true;
}


/*************************************************************************

@brief The run function sets kernel arguments and then launches
the kernels.  Data is then read from the kernels write buffer.

@return void

**************************************************************************/
void run() {
	cl_int status;
	cl_event kernel_event;
	cl_event finish_event;
	cl_event events[8];

	const double start_time = getCurrentTimestamp();

	// Copy data from host to device
	status = clEnqueueWriteBuffer(queue[K_READER], input_noisy_message_I_buf, CL_TRUE,
		0, (SLOT_LEN)* sizeof(char), noisyDataIn_I, 0, NULL, &events[0]);
	checkError(status, "Failed to transfer input noisy I data");
	// Try putting Q data in same buffer, but moved over 90 bytes
	status = clEnqueueWriteBuffer(queue[K_READER], input_noisy_message_Q_buf, CL_TRUE,
		0, (SLOT_LEN)* sizeof(char), noisyDataIn_Q, 0, NULL, &events[0]);
	checkError(status, "Failed to transfer input  noisy Q data");
	

	//***********************************
	// Set kernel arguments.
	//***********************************
	//SNR Estimation Reader Kernel
	status = clSetKernelArg(kernel[K_READER], 0, sizeof(cl_mem), &input_noisy_message_I_buf);
	checkError(status, "Failed to set input reader arg 0");
	status = clSetKernelArg(kernel[K_READER], 1, sizeof(cl_mem), &input_noisy_message_Q_buf);
	checkError(status, "Failed to set input reader arg 1");
	status = clSetKernelArg(kernel[K_READER], 2, sizeof(int), &input_file_size);
	checkError(status, "Failed to set K_READER arg 2");

	//SNR Estimation
	status = clSetKernelArg(kernel[K_SNR_EST_LUT_CORRECTION], 0, sizeof(unsigned int), &slotLen);
	checkError(status, "Failed to set K_SNR_EST_LUT_CORRECTION arg 0");	

	//SNR Estimation Writer Kernel
	status = clSetKernelArg(kernel[K_WRITER], 0, sizeof(cl_mem), &output_buf);  //store final SNR estimate here
	checkError(status, "Failed to set K_WRITER arg 0");
	status = clSetKernelArg(kernel[K_WRITER], 1, sizeof(int), &num_output_frames);
	checkError(status, "Failed to set K_WRITER arg 1");
	
	
	//***********************************
	// Enqueue kernel.
	//*********************************** 
	status = clEnqueueTask(queue[K_READER], kernel[K_READER], 0, NULL, NULL);
	checkError(status, "Failed to launch K_READER");

	//Filter
	status = clEnqueueTask(queue[K_SNR_EST_LUT_CORRECTION], kernel[K_SNR_EST_LUT_CORRECTION], 0, NULL, NULL);
	checkError(status, "Failed to launch K_SNR_EST_LUT_CORRECTION");
	
	// Write
	status = clEnqueueTask(queue[K_WRITER], kernel[K_WRITER], 0, NULL, &kernel_event);
	checkError(status, "Failed to launch K_WRITER");

	printf("Before clFinish\n");

	//***************************************************
	// Wait for command queue to complete pending events
	//***************************************************
	status = clFinish(queue[K_WRITER]);
	checkError(status, "Failed to finish (%d: %s)", K_WRITER, kernel_names[K_WRITER]);
	
	printf("Before clEnqueueReadBuffer\n");
	//********************************************
	// Read the result. This the final operation.
	//********************************************
	status = clEnqueueReadBuffer(queue[K_WRITER], output_buf, CL_TRUE,
		0, num_output_frames*sizeof(short), dout_snr_est, 0, NULL, NULL);


	// Get kernel times using the OpenCL event profiling API.
	cl_ulong time_ns = getStartEndTime(kernel_event);
	#ifdef VERBOSE
		printf("Kernel time: %0.3f ms\n", double(time_ns) * 1e-6);
	#else
		printf("%f\t", (numFrames*slotLen)/(double(time_ns)* 1e-9));
	#endif
}


/*********************************************************

@brief Free the resources allocated during initialization

@return void
*********************************************************/
void cleanup() {

	for (int i = 0; i<K_NUM_KERNELS; ++i) {
		if (kernel[i])
			clReleaseKernel(kernel[i]);
	}
	for (int i = 0; i<K_NUM_KERNELS; ++i) {
		if (queue[i])
			clReleaseCommandQueue(queue[i]);
	}
	
	if (input_noisy_message_I_buf && input_noisy_message_I_buf) {
		clReleaseMemObject(input_noisy_message_I_buf);
	}
	
	if (input_noisy_message_Q_buf && input_noisy_message_Q_buf) {
		clReleaseMemObject(input_noisy_message_Q_buf);
	}
	
	if (output_buf && output_buf) {
		clReleaseMemObject(output_buf);
	}
	if (program) {
		clReleaseProgram(program);
	}
	if (context) {
		clReleaseContext(context);
	}
}
