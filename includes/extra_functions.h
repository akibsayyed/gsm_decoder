/*
 * extra_functions.cpp
 *
 *  Created on: 23-May-2016
 *      Author: akib
 */

#ifndef EXTRA_FUNCTIONS_CPP_
#define EXTRA_FUNCTIONS_CPP_
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef _WIN32
#include <unistd.h>
#else
#include <windows.h>
#include <fcntl.h>
#include <io.h>
#define _USE_MATH_DEFINES
#endif

#include <math.h>
#include <receiver_config.h>
#include <gnuradio/gr_complex.h>


#define DEFAULT_BUF_LENGTH		(16 * 16384)
#define SYNC_SEARCH_RANGE 30




double atofs(char *s);
inline float compute_phase_diff(gr_complex val1, gr_complex val2);
int consume_each(gr_complex *in,int to_consume);
bool reach_sch_burst(gr_complex *in,const int nitems);
int get_sch_chan_imp_resp(const gr_complex *input, gr_complex * chan_imp_resp);
gr_complex correlate_sequence(const gr_complex * sequence, int length, const gr_complex * input);
void detect_burst(const gr_complex * input, gr_complex * chan_imp_resp, int burst_start, unsigned char * output_binary);
double compute_freq_offset(const gr_complex * input, unsigned first_sample, unsigned last_sample);
inline void autocorrelation(const gr_complex * input, gr_complex * out, int nitems);
inline void mafi(const gr_complex * input, int nitems, gr_complex * filter, int filter_length, gr_complex * output);
void viterbi_detector(const gr_complex * input, unsigned int samples_num, gr_complex * rhh, unsigned int start_state, const unsigned int * stop_states, unsigned int stops_num, float * output);
int get_norm_chan_imp_resp(const gr_complex *input, gr_complex * chan_imp_resp, float *corr_max, int bcc);
void send_burst(burst_counter burst_nr, const unsigned char * burst_binary);





#endif /* EXTRA_FUNCTIONS_CPP_ */
