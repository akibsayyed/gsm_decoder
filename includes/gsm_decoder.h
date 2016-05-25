

#ifndef GSM_DECODER_H_
#define GSM_DECODER_H_
#include <stdlib.h>
#include<stdio.h>
#include<unistd.h>
#include<string.h>
#include <sys/types.h>
#include <signal.h>
#include <gnuradio/gr_complex.h>
#include<vector>
#include<liquid/liquid.h>
#include <boost/circular_buffer.hpp>
#include <math.h>
#include <gnuradio/math.h>
#include <algorithm>
#include <numeric>
#include <stdint.h>
#include <sch.h>

/*
 * Local functions header files
 */
#include <extra_functions.h>
#include <gsm_constants.h>
#include <gsm_decoder.h>
#include <receiver_config.h>
#include <receiver.h>
states d_state=fcch_search;

#define DEFAULT_SAMPLE_RATE		1000000
#define DEFAULT_BUF_LENGTH		(16 * 16384)
#define MINIMAL_BUF_LENGTH		512
#define MAXIMAL_BUF_LENGTH		(256 * 16384)


int abc;
int samples_len=0;
int d_counter_bursts;
void parse_options(int argc,char **argv);
void usage();
#endif /* TEMP_H_ */
