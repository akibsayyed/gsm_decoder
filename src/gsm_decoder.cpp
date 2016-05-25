
#include <stdio.h>
#include <stdlib.h>
#include <gsm_decoder.h>

uint32_t samp_rate = DEFAULT_SAMPLE_RATE;
uint32_t out_block_size = DEFAULT_BUF_LENGTH;
char *filename="";
FILE *file=NULL;


#include<extra_functions.h>
#include <gsm_constants.h>
#include <unistd.h>
#include <math.h>
#include <gnuradio/math.h>
#include <algorithm>
#include <numeric>
#include <gnuradio/gr_complex.h>
#include <cmath>
#include <receiver_config.h>
double atofs(char *s)
/* standard suffixes */
{
	char last;
	int len;
	double suff = 1.0;
	len = strlen(s);
	last = s[len-1];
	s[len-1] = '\0';
	switch (last) {
		case 'g':
		case 'G':
			suff *= 1e3;
		case 'm':
		case 'M':
			suff *= 1e3;
		case 'k':
		case 'K':
			suff *= 1e3;
			suff *= atof(s);
			s[len-1] = last;
			return suff;
	}
	s[len-1] = last;
	return atof(s);
}
inline float compute_phase_diff(gr_complex val1, gr_complex val2)
{
    gr_complex conjprod = val1 * conj(val2);
   // return fast_atan2f(imag(conjprod), real(conjprod));
    return gr::fast_atan2f(imag(conjprod), real(conjprod));

}

int consume_each(gr_complex *in,int to_consume){
	int i;
	int counter=to_consume;
	for (i=0;i<DEFAULT_BUF_LENGTH/2;i++){
		in[i]=in[counter];
		counter++;
	}
	return to_consume;

}

bool reach_sch_burst(gr_complex *in,const int nitems)
{
    //it just consumes samples to get near to a SCH burst
    int to_consume = 0;
    bool result = false;
    unsigned sample_nr_near_sch_start = d_fcch_start_pos + (FRAME_BITS - SAFETY_MARGIN) * d_OSR;

    //consume samples until d_counter will be equal to sample_nr_near_sch_start
    if (d_counter < sample_nr_near_sch_start)
    {
        if (d_counter + nitems >= sample_nr_near_sch_start)
        {
            to_consume = sample_nr_near_sch_start - d_counter;
        }
        else
        {
            to_consume = nitems;
        }
        result = false;
    }
    else
    {
        to_consume = 0;
        result = true;
    }
fprintf(stdout,"\nTo Consume %d",to_consume);
    d_counter += to_consume;
    samples_len=nitems-to_consume;
    consume_each(in,to_consume);
    return result;
}





int get_sch_chan_imp_resp(const gr_complex *input, gr_complex * chan_imp_resp)
{
    std::vector<gr_complex> correlation_buffer;
    std::vector<float> power_buffer;
    std::vector<float> window_energy_buffer;

    int strongest_window_nr;
    int burst_start = 0;
    int chan_imp_resp_center = 0;
    float max_correlation = 0;
    float energy = 0;

    for (int ii = SYNC_POS * d_OSR; ii < (SYNC_POS + SYNC_SEARCH_RANGE) *d_OSR; ii++)
    {
        gr_complex correlation = correlate_sequence(&d_sch_training_seq[5], N_SYNC_BITS - 10, &input[ii]);
        correlation_buffer.push_back(correlation);
        power_buffer.push_back(std::pow(abs(correlation), 2));
    }
    //compute window energies
    std::vector<float>::iterator iter = power_buffer.begin();
    bool loop_end = false;
    while (iter != power_buffer.end())
    {
        std::vector<float>::iterator iter_ii = iter;
        energy = 0;

        for (int ii = 0; ii < (d_chan_imp_length) *d_OSR; ii++, iter_ii++)
        {
            if (iter_ii == power_buffer.end())
            {
                loop_end = true;
                break;
            }
            energy += (*iter_ii);
        }
        if (loop_end)
        {
            break;
        }
        iter++;
        window_energy_buffer.push_back(energy);
    }

    strongest_window_nr = max_element(window_energy_buffer.begin(), window_energy_buffer.end()) - window_energy_buffer.begin();
    //   d_channel_imp_resp.clear();

    max_correlation = 0;
    for (int ii = 0; ii < (d_chan_imp_length) *d_OSR; ii++)
    {
        gr_complex correlation = correlation_buffer[strongest_window_nr + ii];
        if (abs(correlation) > max_correlation)
        {
            chan_imp_resp_center = ii;
            max_correlation = abs(correlation);
        }
        //     d_channel_imp_resp.push_back(correlation);
        chan_imp_resp[ii] = correlation;
    }

    burst_start = strongest_window_nr + chan_imp_resp_center - 48 * d_OSR - 2 * d_OSR + 2 + SYNC_POS * d_OSR;
    int test=0;
    return burst_start;
}



gr_complex correlate_sequence(const gr_complex * sequence, int length, const gr_complex * input)
{
    gr_complex result(0.0, 0.0);
    int sample_number = 0;

    for (int ii = 0; ii < length; ii++)
    {
        sample_number = (ii * d_OSR) ;
        result += sequence[ii] * conj(input[sample_number]);
    }

    result = result / gr_complex(length, 0);
    return result;
}


void detect_burst(const gr_complex * input, gr_complex * chan_imp_resp, int burst_start, unsigned char * output_binary)
{
    float output[BURST_SIZE];
    std::vector<gr_complex> rhh_temp(CHAN_IMP_RESP_LENGTH*d_OSR);
    gr_complex rhh[CHAN_IMP_RESP_LENGTH];
    gr_complex filtered_burst[BURST_SIZE];
    int start_state = 3;
    unsigned int stop_states[2] = {4, 12};

    autocorrelation(chan_imp_resp, &rhh_temp[0], d_chan_imp_length*d_OSR);
    for (int ii = 0; ii < (d_chan_imp_length); ii++)
    {
        rhh[ii] = conj(rhh_temp[ii*d_OSR]);
    }

    mafi(&input[burst_start], BURST_SIZE, chan_imp_resp, d_chan_imp_length*d_OSR, filtered_burst);

    viterbi_detector(filtered_burst, BURST_SIZE, rhh, start_state, stop_states, 2, output);

    //fprintf(stdout,"\n");
    for (int i = 0; i < BURST_SIZE ; i++)
    {
    	//fprintf(stdout,"%x",(output[i] > 0));
        output_binary[i] = (output[i] > 0);
    }
}


double compute_freq_offset(const gr_complex * input, unsigned first_sample, unsigned last_sample)
{
    double phase_sum = 0;
    unsigned ii;

    for (ii = first_sample; ii < last_sample; ii++)
    {
        double phase_diff = compute_phase_diff(input[ii], input[ii-1]) - (M_PI / 2) / d_OSR;
        phase_sum += phase_diff;
    }

    double phase_offset = phase_sum / (last_sample - first_sample);
    double freq_offset = phase_offset * 1625000.0 / (12.0 * M_PI);
    return freq_offset;
}


//computes autocorrelation for positive arguments
inline void autocorrelation(const gr_complex * input, gr_complex * out, int nitems)
{
    int i, k;
    for (k = nitems - 1; k >= 0; k--)
    {
        out[k] = gr_complex(0, 0);
        for (i = k; i < nitems; i++)
        {
            out[k] += input[i] * conj(input[i-k]);
        }
    }
}


inline void mafi(const gr_complex * input, int nitems, gr_complex * filter, int filter_length, gr_complex * output)
{
    int ii = 0, n, a;

    for (n = 0; n < nitems; n++)
    {
        a = n * d_OSR;
        output[n] = 0;
        ii = 0;

        while (ii < filter_length)
        {
            if ((a + ii) >= nitems*d_OSR){
                break;
            }
            output[n] += input[a+ii] * filter[ii];
            ii++;
        }
    }
}





/* -*- c++ -*- */
/*
 * @file
 * @author Piotr Krysik <ptrkrysik@gmail.com>
 * @section LICENSE
 *
 * Gr-gsm is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * Gr-gsm is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with gr-gsm; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

/*
 * viterbi_detector:
 *           This part does the detection of received sequnece.
 *           Employed algorithm is viterbi Maximum Likehood Sequence Estimation.
 *           At this moment it gives hard decisions on the output, but
 *           it was designed with soft decisions in mind.
 *
 * SYNTAX:   void viterbi_detector(
 *                                  const gr_complex * input,
 *                                  unsigned int samples_num,
 *                                  gr_complex * rhh,
 *                                  unsigned int start_state,
 *                                  const unsigned int * stop_states,
 *                                  unsigned int stops_num,
 *                                  float * output)
 *
 * INPUT:    input:       Complex received signal afted matched filtering.
 *           samples_num: Number of samples in the input table.
 *           rhh:         The autocorrelation of the estimated channel
 *                        impulse response.
 *           start_state: Number of the start point. In GSM each burst
 *                        starts with sequence of three bits (0,0,0) which
 *                        indicates start point of the algorithm.
 *           stop_states: Table with numbers of possible stop states.
 *           stops_num:   Number of possible stop states
 *
 *
 * OUTPUT:   output:      Differentially decoded hard output of the algorithm:
 *                        -1 for logical "0" and 1 for logical "1"
 *
 * SUB_FUNC: none
 *
 * TEST(S):  Tested with real world normal burst.
 */



#define PATHS_NUM (1 << (CHAN_IMP_RESP_LENGTH-1))

void viterbi_detector(const gr_complex * input, unsigned int samples_num, gr_complex * rhh, unsigned int start_state, const unsigned int * stop_states, unsigned int stops_num, float * output)
{
   float increment[8];
   float path_metrics1[16];
   float path_metrics2[16];
   float paths_difference;
   float * new_path_metrics;
   float * old_path_metrics;
   float * tmp;
   float trans_table[BURST_SIZE][16];
   float pm_candidate1, pm_candidate2;
   bool real_imag;
   float input_symbol_real, input_symbol_imag;
   unsigned int i, sample_nr;

/*
* Setup first path metrics, so only state pointed by start_state is possible.
* Start_state metric is equal to zero, the rest is written with some very low value,
* which makes them practically impossible to occur.
*/
   for(i=0; i<PATHS_NUM; i++){
      path_metrics1[i]=(-10e30);
   }
   path_metrics1[start_state]=0;

/*
* Compute Increment - a table of values which does not change for subsequent input samples.
* Increment is table of reference levels for computation of branch metrics:
*    branch metric = (+/-)received_sample (+/-) reference_level
*/
   increment[0] = -rhh[1].imag() -rhh[2].real() -rhh[3].imag() +rhh[4].real();
   increment[1] = rhh[1].imag() -rhh[2].real() -rhh[3].imag() +rhh[4].real();
   increment[2] = -rhh[1].imag() +rhh[2].real() -rhh[3].imag() +rhh[4].real();
   increment[3] = rhh[1].imag() +rhh[2].real() -rhh[3].imag() +rhh[4].real();
   increment[4] = -rhh[1].imag() -rhh[2].real() +rhh[3].imag() +rhh[4].real();
   increment[5] = rhh[1].imag() -rhh[2].real() +rhh[3].imag() +rhh[4].real();
   increment[6] = -rhh[1].imag() +rhh[2].real() +rhh[3].imag() +rhh[4].real();
   increment[7] = rhh[1].imag() +rhh[2].real() +rhh[3].imag() +rhh[4].real();


/*
* Computation of path metrics and decisions (Add-Compare-Select).
* It's composed of two parts: one for odd input samples (imaginary numbers)
* and one for even samples (real numbers).
* Each part is composed of independent (parallelisable) statements like
* this one:
*      pm_candidate1 = old_path_metrics[0] -input_symbol_imag +increment[2];
*      pm_candidate2 = old_path_metrics[8] -input_symbol_imag -increment[5];
*      paths_difference=pm_candidate2-pm_candidate1;
*      new_path_metrics[1]=(paths_difference<0) ? pm_candidate1 : pm_candidate2;
*      trans_table[sample_nr][1] = paths_difference;
* This is very good point for optimisations (SIMD or OpenMP) as it's most time
* consuming part of this function.
*/
   sample_nr=0;
   old_path_metrics=path_metrics1;
   new_path_metrics=path_metrics2;
   while(sample_nr<samples_num){
      //Processing imag states
      real_imag=1;
      input_symbol_imag = input[sample_nr].imag();

      pm_candidate1 = old_path_metrics[0] +input_symbol_imag -increment[2];
      pm_candidate2 = old_path_metrics[8] +input_symbol_imag +increment[5];
      paths_difference=pm_candidate2-pm_candidate1;
      new_path_metrics[0]=(paths_difference<0) ? pm_candidate1 : pm_candidate2;
      trans_table[sample_nr][0] = paths_difference;

      pm_candidate1 = old_path_metrics[0] -input_symbol_imag +increment[2];
      pm_candidate2 = old_path_metrics[8] -input_symbol_imag -increment[5];
      paths_difference=pm_candidate2-pm_candidate1;
      new_path_metrics[1]=(paths_difference<0) ? pm_candidate1 : pm_candidate2;
      trans_table[sample_nr][1] = paths_difference;

      pm_candidate1 = old_path_metrics[1] +input_symbol_imag -increment[3];
      pm_candidate2 = old_path_metrics[9] +input_symbol_imag +increment[4];
      paths_difference=pm_candidate2-pm_candidate1;
      new_path_metrics[2]=(paths_difference<0) ? pm_candidate1 : pm_candidate2;
      trans_table[sample_nr][2] = paths_difference;

      pm_candidate1 = old_path_metrics[1] -input_symbol_imag +increment[3];
      pm_candidate2 = old_path_metrics[9] -input_symbol_imag -increment[4];
      paths_difference=pm_candidate2-pm_candidate1;
      new_path_metrics[3]=(paths_difference<0) ? pm_candidate1 : pm_candidate2;
      trans_table[sample_nr][3] = paths_difference;

      pm_candidate1 = old_path_metrics[2] +input_symbol_imag -increment[0];
      pm_candidate2 = old_path_metrics[10] +input_symbol_imag +increment[7];
      paths_difference=pm_candidate2-pm_candidate1;
      new_path_metrics[4]=(paths_difference<0) ? pm_candidate1 : pm_candidate2;
      trans_table[sample_nr][4] = paths_difference;

      pm_candidate1 = old_path_metrics[2] -input_symbol_imag +increment[0];
      pm_candidate2 = old_path_metrics[10] -input_symbol_imag -increment[7];
      paths_difference=pm_candidate2-pm_candidate1;
      new_path_metrics[5]=(paths_difference<0) ? pm_candidate1 : pm_candidate2;
      trans_table[sample_nr][5] = paths_difference;

      pm_candidate1 = old_path_metrics[3] +input_symbol_imag -increment[1];
      pm_candidate2 = old_path_metrics[11] +input_symbol_imag +increment[6];
      paths_difference=pm_candidate2-pm_candidate1;
      new_path_metrics[6]=(paths_difference<0) ? pm_candidate1 : pm_candidate2;
      trans_table[sample_nr][6] = paths_difference;

      pm_candidate1 = old_path_metrics[3] -input_symbol_imag +increment[1];
      pm_candidate2 = old_path_metrics[11] -input_symbol_imag -increment[6];
      paths_difference=pm_candidate2-pm_candidate1;
      new_path_metrics[7]=(paths_difference<0) ? pm_candidate1 : pm_candidate2;
      trans_table[sample_nr][7] = paths_difference;

      pm_candidate1 = old_path_metrics[4] +input_symbol_imag -increment[6];
      pm_candidate2 = old_path_metrics[12] +input_symbol_imag +increment[1];
      paths_difference=pm_candidate2-pm_candidate1;
      new_path_metrics[8]=(paths_difference<0) ? pm_candidate1 : pm_candidate2;
      trans_table[sample_nr][8] = paths_difference;

      pm_candidate1 = old_path_metrics[4] -input_symbol_imag +increment[6];
      pm_candidate2 = old_path_metrics[12] -input_symbol_imag -increment[1];
      paths_difference=pm_candidate2-pm_candidate1;
      new_path_metrics[9]=(paths_difference<0) ? pm_candidate1 : pm_candidate2;
      trans_table[sample_nr][9] = paths_difference;

      pm_candidate1 = old_path_metrics[5] +input_symbol_imag -increment[7];
      pm_candidate2 = old_path_metrics[13] +input_symbol_imag +increment[0];
      paths_difference=pm_candidate2-pm_candidate1;
      new_path_metrics[10]=(paths_difference<0) ? pm_candidate1 : pm_candidate2;
      trans_table[sample_nr][10] = paths_difference;

      pm_candidate1 = old_path_metrics[5] -input_symbol_imag +increment[7];
      pm_candidate2 = old_path_metrics[13] -input_symbol_imag -increment[0];
      paths_difference=pm_candidate2-pm_candidate1;
      new_path_metrics[11]=(paths_difference<0) ? pm_candidate1 : pm_candidate2;
      trans_table[sample_nr][11] = paths_difference;

      pm_candidate1 = old_path_metrics[6] +input_symbol_imag -increment[4];
      pm_candidate2 = old_path_metrics[14] +input_symbol_imag +increment[3];
      paths_difference=pm_candidate2-pm_candidate1;
      new_path_metrics[12]=(paths_difference<0) ? pm_candidate1 : pm_candidate2;
      trans_table[sample_nr][12] = paths_difference;

      pm_candidate1 = old_path_metrics[6] -input_symbol_imag +increment[4];
      pm_candidate2 = old_path_metrics[14] -input_symbol_imag -increment[3];
      paths_difference=pm_candidate2-pm_candidate1;
      new_path_metrics[13]=(paths_difference<0) ? pm_candidate1 : pm_candidate2;
      trans_table[sample_nr][13] = paths_difference;

      pm_candidate1 = old_path_metrics[7] +input_symbol_imag -increment[5];
      pm_candidate2 = old_path_metrics[15] +input_symbol_imag +increment[2];
      paths_difference=pm_candidate2-pm_candidate1;
      new_path_metrics[14]=(paths_difference<0) ? pm_candidate1 : pm_candidate2;
      trans_table[sample_nr][14] = paths_difference;

      pm_candidate1 = old_path_metrics[7] -input_symbol_imag +increment[5];
      pm_candidate2 = old_path_metrics[15] -input_symbol_imag -increment[2];
      paths_difference=pm_candidate2-pm_candidate1;
      new_path_metrics[15]=(paths_difference<0) ? pm_candidate1 : pm_candidate2;
      trans_table[sample_nr][15] = paths_difference;
      tmp=old_path_metrics;
      old_path_metrics=new_path_metrics;
      new_path_metrics=tmp;

      sample_nr++;
      if(sample_nr==samples_num)
         break;

      //Processing real states
      real_imag=0;
      input_symbol_real = input[sample_nr].real();

      pm_candidate1 = old_path_metrics[0] -input_symbol_real -increment[7];
      pm_candidate2 = old_path_metrics[8] -input_symbol_real +increment[0];
      paths_difference=pm_candidate2-pm_candidate1;
      new_path_metrics[0]=(paths_difference<0) ? pm_candidate1 : pm_candidate2;
      trans_table[sample_nr][0] = paths_difference;

      pm_candidate1 = old_path_metrics[0] +input_symbol_real +increment[7];
      pm_candidate2 = old_path_metrics[8] +input_symbol_real -increment[0];
      paths_difference=pm_candidate2-pm_candidate1;
      new_path_metrics[1]=(paths_difference<0) ? pm_candidate1 : pm_candidate2;
      trans_table[sample_nr][1] = paths_difference;

      pm_candidate1 = old_path_metrics[1] -input_symbol_real -increment[6];
      pm_candidate2 = old_path_metrics[9] -input_symbol_real +increment[1];
      paths_difference=pm_candidate2-pm_candidate1;
      new_path_metrics[2]=(paths_difference<0) ? pm_candidate1 : pm_candidate2;
      trans_table[sample_nr][2] = paths_difference;

      pm_candidate1 = old_path_metrics[1] +input_symbol_real +increment[6];
      pm_candidate2 = old_path_metrics[9] +input_symbol_real -increment[1];
      paths_difference=pm_candidate2-pm_candidate1;
      new_path_metrics[3]=(paths_difference<0) ? pm_candidate1 : pm_candidate2;
      trans_table[sample_nr][3] = paths_difference;

      pm_candidate1 = old_path_metrics[2] -input_symbol_real -increment[5];
      pm_candidate2 = old_path_metrics[10] -input_symbol_real +increment[2];
      paths_difference=pm_candidate2-pm_candidate1;
      new_path_metrics[4]=(paths_difference<0) ? pm_candidate1 : pm_candidate2;
      trans_table[sample_nr][4] = paths_difference;

      pm_candidate1 = old_path_metrics[2] +input_symbol_real +increment[5];
      pm_candidate2 = old_path_metrics[10] +input_symbol_real -increment[2];
      paths_difference=pm_candidate2-pm_candidate1;
      new_path_metrics[5]=(paths_difference<0) ? pm_candidate1 : pm_candidate2;
      trans_table[sample_nr][5] = paths_difference;

      pm_candidate1 = old_path_metrics[3] -input_symbol_real -increment[4];
      pm_candidate2 = old_path_metrics[11] -input_symbol_real +increment[3];
      paths_difference=pm_candidate2-pm_candidate1;
      new_path_metrics[6]=(paths_difference<0) ? pm_candidate1 : pm_candidate2;
      trans_table[sample_nr][6] = paths_difference;

      pm_candidate1 = old_path_metrics[3] +input_symbol_real +increment[4];
      pm_candidate2 = old_path_metrics[11] +input_symbol_real -increment[3];
      paths_difference=pm_candidate2-pm_candidate1;
      new_path_metrics[7]=(paths_difference<0) ? pm_candidate1 : pm_candidate2;
      trans_table[sample_nr][7] = paths_difference;

      pm_candidate1 = old_path_metrics[4] -input_symbol_real -increment[3];
      pm_candidate2 = old_path_metrics[12] -input_symbol_real +increment[4];
      paths_difference=pm_candidate2-pm_candidate1;
      new_path_metrics[8]=(paths_difference<0) ? pm_candidate1 : pm_candidate2;
      trans_table[sample_nr][8] = paths_difference;

      pm_candidate1 = old_path_metrics[4] +input_symbol_real +increment[3];
      pm_candidate2 = old_path_metrics[12] +input_symbol_real -increment[4];
      paths_difference=pm_candidate2-pm_candidate1;
      new_path_metrics[9]=(paths_difference<0) ? pm_candidate1 : pm_candidate2;
      trans_table[sample_nr][9] = paths_difference;

      pm_candidate1 = old_path_metrics[5] -input_symbol_real -increment[2];
      pm_candidate2 = old_path_metrics[13] -input_symbol_real +increment[5];
      paths_difference=pm_candidate2-pm_candidate1;
      new_path_metrics[10]=(paths_difference<0) ? pm_candidate1 : pm_candidate2;
      trans_table[sample_nr][10] = paths_difference;

      pm_candidate1 = old_path_metrics[5] +input_symbol_real +increment[2];
      pm_candidate2 = old_path_metrics[13] +input_symbol_real -increment[5];
      paths_difference=pm_candidate2-pm_candidate1;
      new_path_metrics[11]=(paths_difference<0) ? pm_candidate1 : pm_candidate2;
      trans_table[sample_nr][11] = paths_difference;

      pm_candidate1 = old_path_metrics[6] -input_symbol_real -increment[1];
      pm_candidate2 = old_path_metrics[14] -input_symbol_real +increment[6];
      paths_difference=pm_candidate2-pm_candidate1;
      new_path_metrics[12]=(paths_difference<0) ? pm_candidate1 : pm_candidate2;
      trans_table[sample_nr][12] = paths_difference;

      pm_candidate1 = old_path_metrics[6] +input_symbol_real +increment[1];
      pm_candidate2 = old_path_metrics[14] +input_symbol_real -increment[6];
      paths_difference=pm_candidate2-pm_candidate1;
      new_path_metrics[13]=(paths_difference<0) ? pm_candidate1 : pm_candidate2;
      trans_table[sample_nr][13] = paths_difference;

      pm_candidate1 = old_path_metrics[7] -input_symbol_real -increment[0];
      pm_candidate2 = old_path_metrics[15] -input_symbol_real +increment[7];
      paths_difference=pm_candidate2-pm_candidate1;
      new_path_metrics[14]=(paths_difference<0) ? pm_candidate1 : pm_candidate2;
      trans_table[sample_nr][14] = paths_difference;

      pm_candidate1 = old_path_metrics[7] +input_symbol_real +increment[0];
      pm_candidate2 = old_path_metrics[15] +input_symbol_real -increment[7];
      paths_difference=pm_candidate2-pm_candidate1;
      new_path_metrics[15]=(paths_difference<0) ? pm_candidate1 : pm_candidate2;
      trans_table[sample_nr][15] = paths_difference;

      tmp=old_path_metrics;
      old_path_metrics=new_path_metrics;
      new_path_metrics=tmp;

      sample_nr++;
   }

/*
* Find the best from the stop states by comparing their path metrics.
* Not every stop state is always possible, so we are searching in
* a subset of them.
*/
   unsigned int best_stop_state;
   float stop_state_metric, max_stop_state_metric;
   best_stop_state = stop_states[0];
   max_stop_state_metric = old_path_metrics[best_stop_state];
   for(i=1; i< stops_num; i++){
      stop_state_metric = old_path_metrics[stop_states[i]];
      if(stop_state_metric > max_stop_state_metric){
         max_stop_state_metric = stop_state_metric;
         best_stop_state = stop_states[i];
      }
   }

/*
* This table was generated with hope that it gives a litle speedup during
* traceback stage.
* Received bit is related to the number of state in the trellis.
* I've numbered states so their parity (number of ones) is related
* to a received bit.
*/
   static const unsigned int parity_table[PATHS_NUM] = { 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0,  };

/*
* Table of previous states in the trellis diagram.
* For GMSK modulation every state has two previous states.
* Example:
*   previous_state_nr1 = prev_table[current_state_nr][0]
*   previous_state_nr2 = prev_table[current_state_nr][1]
*/
   static const unsigned int prev_table[PATHS_NUM][2] = { {0,8}, {0,8}, {1,9}, {1,9}, {2,10}, {2,10}, {3,11}, {3,11}, {4,12}, {4,12}, {5,13}, {5,13}, {6,14}, {6,14}, {7,15}, {7,15},  };

/*
* Traceback and differential decoding of received sequence.
* Decisions stored in trans_table are used to restore best path in the trellis.
*/
   sample_nr=samples_num;
   unsigned int state_nr=best_stop_state;
   unsigned int decision;
   bool out_bit=0;

   while(sample_nr>0){
      sample_nr--;
      decision = (trans_table[sample_nr][state_nr]>0);

      if(decision != out_bit)
         output[sample_nr]=-trans_table[sample_nr][state_nr];
      else
         output[sample_nr]=trans_table[sample_nr][state_nr];

      out_bit = out_bit ^ real_imag ^ parity_table[state_nr];
      state_nr = prev_table[state_nr][decision];
      real_imag = !real_imag;
   }
}


//especially computations of strongest_window_nr
int get_norm_chan_imp_resp(const gr_complex *input, gr_complex * chan_imp_resp, float *corr_max, int bcc)
{
    std::vector<gr_complex> correlation_buffer;
    std::vector<float> power_buffer;
    std::vector<float> window_energy_buffer;

    int strongest_window_nr;
    int burst_start = 0;
    int chan_imp_resp_center = 0;
    float max_correlation = 0;
    float energy = 0;

    int search_center = (int)((TRAIN_POS + GUARD_PERIOD) * d_OSR);
    int search_start_pos = search_center + 1 - 5*d_OSR;
    //   int search_start_pos = search_center -  d_chan_imp_length * d_OSR;
    int search_stop_pos = search_center + d_chan_imp_length * d_OSR + 5 * d_OSR;

    for(int ii = search_start_pos; ii < search_stop_pos; ii++)
    {
        gr_complex correlation = correlate_sequence(&d_norm_training_seq[bcc][TRAIN_BEGINNING], N_TRAIN_BITS - 10, &input[ii]);
        correlation_buffer.push_back(correlation);
        power_buffer.push_back(std::pow(abs(correlation), 2));
    }
//    plot(power_buffer);
    //compute window energies
    std::vector<float>::iterator iter = power_buffer.begin();
    bool loop_end = false;
    while (iter != power_buffer.end())
    {
        std::vector<float>::iterator iter_ii = iter;
        energy = 0;

        for (int ii = 0; ii < (d_chan_imp_length - 2)*d_OSR; ii++, iter_ii++)
        {
            if (iter_ii == power_buffer.end())
            {
                loop_end = true;
                break;
            }
            energy += (*iter_ii);
        }
        if (loop_end)
        {
            break;
        }
        iter++;

        window_energy_buffer.push_back(energy);
    }

    strongest_window_nr = max_element(window_energy_buffer.begin(), window_energy_buffer.end()-((d_chan_imp_length)*d_OSR)) - window_energy_buffer.begin();
    //strongest_window_nr = strongest_window_nr-d_OSR;
    if(strongest_window_nr<0){
       strongest_window_nr = 0;
    }

    max_correlation = 0;
    for (int ii = 0; ii < (d_chan_imp_length)*d_OSR; ii++)
    {
        gr_complex correlation = correlation_buffer[strongest_window_nr + ii];
        if (abs(correlation) > max_correlation)
        {
            chan_imp_resp_center = ii;
            max_correlation = abs(correlation);
        }
        //     d_channel_imp_resp.push_back(correlation);
        chan_imp_resp[ii] = correlation;
    }

    *corr_max = max_correlation;

    //DCOUT("strongest_window_nr_new: " << strongest_window_nr);
    burst_start = search_start_pos + strongest_window_nr - TRAIN_POS * d_OSR; //compute first sample posiiton which corresponds to the first sample of the impulse response

    //DCOUT("burst_start: " << burst_start);
    return burst_start;
}


void send_burst(burst_counter burst_nr, const unsigned char * burst_binary){
	fprintf(stdout,"\nframe-nr= %x",d_burst_nr.get_frame_nr());
	for(int i=0;i<BURST_SIZE;i++){
		fprintf(stdout,"%x",burst_binary[i]);
	}
}




/*
 *
 *
 *
 */
void usage(){
	printf("\nEnter proper format ./gsm_decoder -s sample_rate -f filename");
}
void parse_options(int argc,char **argv){


	int r, opt;
	while ((opt = getopt(argc, argv, "d:f:g:s:b:n:p:S:R")) != -1) {
			switch (opt) {

			case 's':
				samp_rate = (uint32_t)atofs(optarg);
				printf("\nChosen sample rate=%d",samp_rate);
				break;
			case 'f':
				filename = optarg;
				printf("\nChosen filename rate=%s",filename);
				break;
			default:
				usage();
				break;
			}
		}

}
void configure_receiver()
{
    d_channel_conf.set_multiframe_type(TIMESLOT0, multiframe_51);
    d_channel_conf.set_burst_types(TIMESLOT0, TEST51, sizeof(TEST51) / sizeof(unsigned), dummy_or_normal);

    d_channel_conf.set_burst_types(TIMESLOT0, TEST_CCH_FRAMES, sizeof(TEST_CCH_FRAMES) / sizeof(unsigned), dummy_or_normal);
    d_channel_conf.set_burst_types(TIMESLOT0, FCCH_FRAMES, sizeof(FCCH_FRAMES) / sizeof(unsigned), fcch_burst);
    d_channel_conf.set_burst_types(TIMESLOT0, SCH_FRAMES, sizeof(SCH_FRAMES) / sizeof(unsigned), sch_burst);

    d_channel_conf.set_multiframe_type(TIMESLOT1, multiframe_51);
    d_channel_conf.set_burst_types(TIMESLOT1, TEST51, sizeof(TEST51) / sizeof(unsigned), dummy_or_normal);
    d_channel_conf.set_multiframe_type(TIMESLOT2, multiframe_51);
    d_channel_conf.set_burst_types(TIMESLOT2, TEST51, sizeof(TEST51) / sizeof(unsigned), dummy_or_normal);
    d_channel_conf.set_multiframe_type(TIMESLOT3, multiframe_51);
    d_channel_conf.set_burst_types(TIMESLOT3, TEST51, sizeof(TEST51) / sizeof(unsigned), dummy_or_normal);
    d_channel_conf.set_multiframe_type(TIMESLOT4, multiframe_51);
    d_channel_conf.set_burst_types(TIMESLOT4, TEST51, sizeof(TEST51) / sizeof(unsigned), dummy_or_normal);
    d_channel_conf.set_multiframe_type(TIMESLOT5, multiframe_51);
    d_channel_conf.set_burst_types(TIMESLOT5, TEST51, sizeof(TEST51) / sizeof(unsigned), dummy_or_normal);
    d_channel_conf.set_multiframe_type(TIMESLOT6, multiframe_51);
    d_channel_conf.set_burst_types(TIMESLOT6, TEST51, sizeof(TEST51) / sizeof(unsigned), dummy_or_normal);
    d_channel_conf.set_multiframe_type(TIMESLOT7, multiframe_51);
    d_channel_conf.set_burst_types(TIMESLOT7, TEST51, sizeof(TEST51) / sizeof(unsigned), dummy_or_normal);
}

void gmsk_mapper(const unsigned char * input, int nitems, gr_complex * gmsk_output, gr_complex start_point){

	  gr_complex j = gr_complex(0.0, 1.0);

	    int current_symbol;
	    int encoded_symbol;
	    int previous_symbol = 2 * input[0] - 1;
	    gmsk_output[0] = start_point;

	    for (int i = 1; i < nitems; i++)
	    {
	        //change bits representation to NRZ
	        current_symbol = 2 * input[i] - 1;
	        //differentially encode
	        encoded_symbol = current_symbol * previous_symbol;
	        //and do gmsk mapping
	        gmsk_output[i] = j * gr_complex(encoded_symbol, 0.0) * gmsk_output[i-1];
	        previous_symbol = current_symbol;
	    }
}
void init_gmsk(){
	 int i;
	                                                                      //don't send samples to the receiver until there are at least samples for one
	   // set_output_multiple(floor((TS_BITS + 2 * GUARD_PERIOD) * d_OSR)); // burst and two gurad periods (one gurard period is an arbitrary overlap)
     //don't send samples to the receiver until there are at least samples for one
	//set_output_multiple(floor((TS_BITS + 2 * GUARD_PERIOD) * d_OSR)); // burst and two gurad periods (one gurard period is an arbitrary overlap)
	gmsk_mapper(SYNC_BITS, N_SYNC_BITS, d_sch_training_seq, gr_complex(0.0, -1.0));
	for (i = 0; i < TRAIN_SEQ_NUM; i++)
	{
		gr_complex startpoint = (train_seq[i][0]==0) ? gr_complex(1.0, 0.0) : gr_complex(-1.0, 0.0); //if first bit of the seqeunce ==0  first symbol ==1
											   //if first bit of the seqeunce ==1  first symbol ==-1
		gmsk_mapper(train_seq[i], N_TRAIN_BITS, d_norm_training_seq[i], startpoint);
	}

	    configure_receiver();  //configure the receiver - tell it where to find which burst type



}



bool find_fcch_burst(gr_complex *input, const int nitems, double & computed_freq_offset)
{
    boost::circular_buffer<float> phase_diff_buffer(FCCH_HITS_NEEDED * d_OSR); //circular buffer used to scan throug signal to find
    //best match for FCCH burst
    float phase_diff = 0;
    gr_complex conjprod;
    int start_pos = -1;
    int hit_count = 0;
    int miss_count = 0;
    float min_phase_diff;
    float max_phase_diff;
    double best_sum = 0;
    float lowest_max_min_diff = 99999;

    int to_consume = 0;
    int sample_number = 0;
    bool end = false;
    bool result = false;
    boost::circular_buffer<float>::iterator buffer_iter;

    /**@name Possible states of FCCH search algorithm*/
    //@{
    enum states
    {
        init,               ///< initialize variables
        search,             ///< search for positive samples
        found_something,    ///< search for FCCH and the best position of it
        fcch_found,         ///< when FCCH was found
        search_fail         ///< when there is no FCCH in the input vector
    } fcch_search_state;
    //@}

    fcch_search_state = init;

    while (!end)
    {
        switch (fcch_search_state)
        {

        case init: //initialize variables
            hit_count = 0;
            miss_count = 0;
            start_pos = -1;
            lowest_max_min_diff = 99999;
            phase_diff_buffer.clear();
            fcch_search_state = search;

            break;

        case search:                                                // search for positive samples
            sample_number++;

            if (sample_number > nitems - FCCH_HITS_NEEDED * d_OSR)   //if it isn't possible to find FCCH because
            {
                                                                       //there's too few samples left to look into,
                to_consume = sample_number;                            //don't do anything with those samples which are left
                                                                       //and consume only those which were checked
                fcch_search_state = search_fail;
            }
            else
            {
                phase_diff = compute_phase_diff(input[sample_number], input[sample_number-1]);

                if (phase_diff > 0)                                   //if a positive phase difference was found
                {
                    to_consume = sample_number;
                    fcch_search_state = found_something;                //switch to state in which searches for FCCH
                }
                else
                {
                    fcch_search_state = search;
                }
            }

            break;

        case found_something:  // search for FCCH and the best position of it
        {
            if (phase_diff > 0)
            {
                hit_count++;       //positive phase differencies increases hits_count
            }
            else
            {
                miss_count++;      //negative increases miss_count
            }

            if ((miss_count >= FCCH_MAX_MISSES * d_OSR) && (hit_count <= FCCH_HITS_NEEDED * d_OSR))
            {
                //if miss_count exceeds limit before hit_count
                fcch_search_state = init;       //go to init
                continue;
            }
            else if (((miss_count >= FCCH_MAX_MISSES * d_OSR) && (hit_count > FCCH_HITS_NEEDED * d_OSR)) || (hit_count > 2 * FCCH_HITS_NEEDED * d_OSR))
            {
                //if hit_count and miss_count exceeds limit then FCCH was found
                fcch_search_state = fcch_found;
                continue;
            }
            else if ((miss_count < FCCH_MAX_MISSES * d_OSR) && (hit_count > FCCH_HITS_NEEDED * d_OSR))
            {
                //find difference between minimal and maximal element in the buffer
                //for FCCH this value should be low
                //this part is searching for a region where this value is lowest
                min_phase_diff = * (min_element(phase_diff_buffer.begin(), phase_diff_buffer.end()));
                max_phase_diff = * (max_element(phase_diff_buffer.begin(), phase_diff_buffer.end()));

                if (lowest_max_min_diff > max_phase_diff - min_phase_diff)
                {
                    lowest_max_min_diff = max_phase_diff - min_phase_diff;
                    start_pos = sample_number - FCCH_HITS_NEEDED * d_OSR - FCCH_MAX_MISSES * d_OSR; //store start pos
                    best_sum = 0;

                    for (buffer_iter = phase_diff_buffer.begin();
                            buffer_iter != (phase_diff_buffer.end());
                            buffer_iter++)
                    {
                        best_sum += *buffer_iter - (M_PI / 2) / d_OSR;   //store best value of phase offset sum
                    }
                }
            }

            sample_number++;

            if (sample_number >= nitems)      //if there's no single sample left to check
            {
                fcch_search_state = search_fail;//FCCH search failed
                continue;
            }

            phase_diff = compute_phase_diff(input[sample_number], input[sample_number-1]);
            phase_diff_buffer.push_back(phase_diff);
            fcch_search_state = found_something;
        }
        break;

        case fcch_found:
        {
            to_consume = start_pos + FCCH_HITS_NEEDED * d_OSR + 1; //consume one FCCH burst

            d_fcch_start_pos = d_counter + start_pos;

            //compute frequency offset
            double phase_offset = best_sum / FCCH_HITS_NEEDED;
            double freq_offset = phase_offset * 1625000.0/6 / (2 * M_PI); //1625000.0/6 - GMSK symbol rate in GSM
            computed_freq_offset = freq_offset;

            end = true;
            result = true;
            break;
        }

        case search_fail:
            end = true;
            result = false;
            break;
        }
    }

    d_counter += to_consume;
    //consume_each(to_consume);
consume_each(input,to_consume);
samples_len=nitems-to_consume;
    //fprintf(stdout,"\nFCCH To Consume %d",to_consume);
    return result;
}

void gsm_decode(gr_complex* input,int noutput_items)
{
	init_gmsk();
	int z=0;
//    std::vector<const gr_complex *> iii = (std::vector<const gr_complex *>) input_items; // jak zrobić to rzutowanie poprawnie
    //gr_complex * input = (gr_complex *) input_items[0];
    //std::vector<tag_t> freq_offset_tags;
    //uint64_t start = nitems_read(0);
    //uint64_t stop = start + noutput_items;

   // float current_time = static_cast<float>(start)/(GSM_SYMBOL_RATE*d_OSR);
   // if((current_time - d_last_time) > 0.1)
   // {
   //     pmt::pmt_t msg = pmt::make_tuple(pmt::mp("current_time"),pmt::from_double(current_time));
   //     message_port_pub(pmt::mp("measurements"), msg);
   //     d_last_time = current_time;
   // }

    //pmt::pmt_t key = pmt::string_to_symbol("setting_freq_offset");
    //get_tags_in_range(freq_offset_tags, 0, start, stop, key);
    //bool freq_offset_tag_in_fcch = false;
    //uint64_t tag_offset=-1; //-1 - just some clearly invalid value

    /*if(!freq_offset_tags.empty()){
        tag_t freq_offset_tag = freq_offset_tags[0];
        tag_offset = freq_offset_tag.offset - start;

        burst_type b_type = d_channel_conf.get_burst_type(d_burst_nr);
        if(d_state == synchronized && b_type == fcch_burst){
            uint64_t last_sample_nr = ceil((GUARD_PERIOD + 2.0 * TAIL_BITS + 156.25) * d_OSR) + 1;
            if(tag_offset < last_sample_nr){
                freq_offset_tag_in_fcch = true;
            }
            d_freq_offset_setting = pmt::to_double(freq_offset_tag.value);
        } else {
            d_freq_offset_setting = pmt::to_double(freq_offset_tag.value);
        }
    }*/
//for (z=0;z<noutput_items;z++){
	switch (d_state)
	    {
	        //bootstrapping
	    case fcch_search:
	    {
	        double freq_offset_tmp;
	        if (find_fcch_burst(input, noutput_items,freq_offset_tmp))
	        {
	            //pmt::pmt_t msg = pmt::make_tuple(pmt::mp("freq_offset"),pmt::from_double(freq_offset_tmp-d_freq_offset_setting),pmt::mp("fcch_search"));
	           // message_port_pub(pmt::mp("measurements"), msg);
	        	fprintf(stderr,"found fcch offset is %f\n",freq_offset_tmp);
	            d_state = sch_search;
	        }
	        else
	        {
	            d_state = fcch_search;
	        }
	        break;
	    }

	    case sch_search:
	    {
	    	//fprintf(stderr,"found sch\n");
	        std::vector<gr_complex> channel_imp_resp(CHAN_IMP_RESP_LENGTH*d_OSR);
	        int t1, t2, t3;
	        int burst_start = 0;
	        unsigned char output_binary[BURST_SIZE];

	        if (reach_sch_burst(input,noutput_items))                                //wait for a SCH burst
	        {
	        	//fprintf(stderr,"inside sch search\n");
	            burst_start = get_sch_chan_imp_resp(input, &channel_imp_resp[0]); //get channel impulse response from it
			printf("d_c0_burst_start %d", burst_start);
				//exit(1);
	            detect_burst(input, &channel_imp_resp[0], burst_start, output_binary); //detect bits using MLSE detection
	            if (decode_sch(&output_binary[3], &t1, &t2, &t3, &d_ncc, &d_bcc) == 0)   //decode SCH burst
	            {
	                d_burst_nr.set(t1, t2, t3, 0);                                  //set counter of bursts value
	                d_burst_nr++;

	                consume_each(input,burst_start + BURST_SIZE * d_OSR + 4*d_OSR);   //consume samples up to next guard period
	                d_state = synchronized;
	                fprintf(stderr,"synchronized\n");
	            }
	            else
	            {
	                d_state = fcch_search;                       //if there is error in the sch burst go back to fcch search phase
	            }
	        }
	        else
	        {
	            d_state = sch_search;
	        }
	        break;
	    }
	    //in this state receiver is synchronized and it processes bursts according to burst type for given burst number
	    case synchronized:
	    {
	    	///fprintf(dump_processed,"=========================");
	    	//dump_complex(input,70000,0);
	    	//ctr++;
	    	//if (ctr>2)
	    	//	exit(1);
	        std::vector<gr_complex> channel_imp_resp(CHAN_IMP_RESP_LENGTH*d_OSR);
	        int offset = 0;
	        int to_consume = 0;
	        unsigned char output_binary[BURST_SIZE];
	        unsigned int inputs_to_process=d_cell_allocation.size();

	        burst_type b_type;

	        /*
	         * for uplink processing
	         */
	        if(0)
	        {
	            inputs_to_process = 2*inputs_to_process;
	        }
	        for(int input_nr=0; input_nr<1; input_nr++)
	        {
	            double signal_pwr = 0;
	            //input = (gr_complex *)input_items[input_nr];

	            for(int ii=GUARD_PERIOD;ii<TS_BITS;ii++)
	            {
	                signal_pwr += abs(input[ii])*abs(input[ii]);
	            }
	            signal_pwr = signal_pwr/(TS_BITS);
	            d_signal_dbm = round(10*log10(signal_pwr/50));
	            if(input_nr==0){
	                d_c0_signal_dbm = d_signal_dbm;
	            }

	            if(input_nr==0) //for c0 channel burst type is controlled by channel configuration
	            {
	                b_type = d_channel_conf.get_burst_type(d_burst_nr); //get burst type for given burst number
	            }
	            else
	            {
	                b_type = normal_or_noise; //for the rest it can be only normal burst or noise (at least at this moment of development)
	            }

	            switch (b_type)
	            {
	            case fcch_burst:                                                                      //if it's FCCH  burst
	            {
	            	fprintf(stderr,"\nfcch_burst!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
	            	d_counter_bursts++;
	                const unsigned first_sample = ceil((GUARD_PERIOD + 2 * TAIL_BITS) * d_OSR) + 1;
	                const unsigned last_sample = first_sample + USEFUL_BITS * d_OSR - TAIL_BITS * d_OSR;
	                double freq_offset_tmp = compute_freq_offset(input, first_sample, last_sample);       //extract frequency offset from it

	                //send_burst(d_burst_nr, fc_fb, GSMTAP_BURST_FCCH, input_nr);

	                //pmt::pmt_t msg = pmt::make_tuple(pmt::mp("freq_offset"),pmt::from_double(freq_offset_tmp-d_freq_offset_setting),pmt::mp("synchronized"));
	                //message_port_pub(pmt::mp("measurements"), msg);
	                break;
	            }
	            case sch_burst:                                                                      //if it's SCH burst
	            {	d_counter_bursts++;
	            fprintf(stderr,"\nsch_burst!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
	                int t1, t2, t3, d_ncc, d_bcc;
	                d_c0_burst_start = get_sch_chan_imp_resp(input, &channel_imp_resp[0]);                //get channel impulse response
			//printf("d_c0_burst_start %d",d_c0_burst_start);
			//exit(1);
	                detect_burst(input, &channel_imp_resp[0], d_c0_burst_start, output_binary);           //MLSE detection of bits
	                //send_burst(d_burst_nr, output_binary, GSMTAP_BURST_SCH, input_nr);
	                if (decode_sch(&output_binary[3], &t1, &t2, &t3, &d_ncc, &d_bcc) == 0)           //and decode SCH data
	                {
	                    // d_burst_nr.set(t1, t2, t3, 0);                                              //but only to check if burst_start value is correct
	                    d_failed_sch = 0;
	                    offset =  d_c0_burst_start - floor((GUARD_PERIOD) * d_OSR);                         //compute offset from burst_start - burst should start after a guard period
	                    to_consume += offset;                                                          //adjust with offset number of samples to be consumed
	                }
	                else
	                {
	                    d_failed_sch++;
	                    if (d_failed_sch >= MAX_SCH_ERRORS)
	                    {
	                        d_state = fcch_search;
	                        //pmt::pmt_t msg = pmt::make_tuple(pmt::mp("freq_offset"),pmt::from_double(0.0),pmt::mp("sync_loss"));
	                        //message_port_pub(pmt::mp("measurements"), msg);
	                        fprintf(stderr,"Re-Synchronization!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
	                    }
	                }
	                break;
	            }
	            case normal_burst:
	            {
	            	d_counter_bursts++;
	            	fprintf(stderr,"\nNormal burst!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
	                float normal_corr_max;                                                    //if it's normal burst
	                d_c0_burst_start = get_norm_chan_imp_resp(input, &channel_imp_resp[0], &normal_corr_max, d_bcc); //get channel impulse response for given training sequence number - d_bcc
	                detect_burst(input, &channel_imp_resp[0], d_c0_burst_start, output_binary);            //MLSE detection of bits
	                //send_burst(d_burst_nr, output_binary, GSMTAP_BURST_NORMAL, input_nr);
	                break;
	            }
	            case dummy_or_normal:
	            {
	            	d_counter_bursts++;
	            	fprintf(stderr,"\nDummy or Normal burst!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
	                unsigned int normal_burst_start, dummy_burst_start;
	                float dummy_corr_max, normal_corr_max;

	                dummy_burst_start = get_norm_chan_imp_resp(input, &channel_imp_resp[0], &dummy_corr_max, TS_DUMMY);
	                normal_burst_start = get_norm_chan_imp_resp(input, &channel_imp_resp[0], &normal_corr_max, d_bcc);

	                if (normal_corr_max > dummy_corr_max)
	                {
	                    d_c0_burst_start = normal_burst_start;
	                    detect_burst(input, &channel_imp_resp[0], normal_burst_start, output_binary);
	                    send_burst(d_burst_nr, output_binary);
	                }
	                else
	                {
	                    d_c0_burst_start = dummy_burst_start;
	                    send_burst(d_burst_nr, dummy_burst);
	                }
	                break;
	            }
	            case rach_burst:
	                break;
	            case dummy:
	            	d_counter_bursts++;
	                //send_burst(d_burst_nr, dummy_burst, GSMTAP_BURST_DUMMY, input_nr);
	                break;
	            case normal_or_noise:
	            {
	            	d_counter_bursts++;
	            	fprintf(stderr,"\n Noise or Normal burst!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
	                unsigned int burst_start;
	                float normal_corr_max_tmp;
	                float normal_corr_max=-1e6;
	                int max_tn;
	                //std::vector<gr_complex> v(input, input + noutput_items);
	                if(d_signal_dbm>=d_c0_signal_dbm-13)
	                {
	                    if(d_tseq_nums.size()==0)              //there is no information about training sequence
	                    {                                      //however the receiver can detect it
	                        get_norm_chan_imp_resp(input, &channel_imp_resp[0], &normal_corr_max, 0);
	                        float ts_max=normal_corr_max;     //with use of a very simple algorithm based on finding
	                        int ts_max_num=0;                 //maximum correlation
	                        for(int ss=1; ss<=7; ss++)
	                        {
	                            get_norm_chan_imp_resp(input, &channel_imp_resp[0], &normal_corr_max, ss);
	                            if(ts_max<normal_corr_max)
	                            {
	                                ts_max = normal_corr_max;
	                                ts_max_num = ss;
	                            }
	                        }
	                        d_tseq_nums.push_back(ts_max_num);
	                    }
	                    int tseq_num;
	                    if(input_nr<=d_tseq_nums.size()){
	                        tseq_num = d_tseq_nums[input_nr-1];
	                    } else {
	                        tseq_num = d_tseq_nums.back();
	                    }
	                    burst_start = get_norm_chan_imp_resp(input, &channel_imp_resp[0], &normal_corr_max, tseq_num);
	//                  if(abs(d_c0_burst_start-burst_start)<=2){ //unused check/filter based on timing
	                    if((normal_corr_max/sqrt(signal_pwr))>=0.9){
	                        detect_burst(input, &channel_imp_resp[0], burst_start, output_binary);
	                        //send_burst(d_burst_nr, output_binary, GSMTAP_BURST_NORMAL, input_nr);
	                    }
	                }
	                break;
	            }
	            case empty:   //if it's empty burst
	                break;      //do nothing
	            }

	          //  if(input_nr==input_items.size()-1)
	          // {
	                d_burst_nr++;   //go to next burst
	                to_consume += TS_BITS * d_OSR + d_burst_nr.get_offset();  //consume samples of the burst up to next guard period
	                fprintf(stderr,"\nto consume in sync-%d---%d",to_consume,d_burst_nr.get_offset());
	                consume_each(input,to_consume);
	          //  }
	            //and add offset which is introduced by
	            //0.25 fractional part of a guard period
	        }
	    }
	    break;
	    }
//}

    //return 0;
}
void resampler(gr_complex *in,int len,double sample_rate,int osr)
{

	float samp_out=1625000.0/6.0*osr;



float r= samp_out/sample_rate;
float As=100.0f; // resampling filter stop-band attenuation [dB]
//unsigned int n= NUM_SAMPLES; // number of input samples

msresamp_crcf q = msresamp_crcf_create(r,As);
//msresamp_crcf_print(q);
float delay = msresamp_crcf_get_delay(q);

// number of input samples (zero-padded)
unsigned int nx = len + (int)ceilf(delay) + 10;

// output buffer with extra padding for good measure
unsigned int ny_alloc = (unsigned int) (2*(float)nx * r);  // allocation for output

printf("nyalloc=%d",ny_alloc);

gr_complex *tempout;
tempout = (gr_complex *)malloc(ny_alloc * sizeof(gr_complex));

//float complex r[nx]; // received signal
//float complex y[ny_alloc]; // resampled signal

//while (1) {
    unsigned int ny;





    //read_complex_iq_block(r, NUM_SAMPLES); // made up function
    msresamp_crcf_execute(q, in, len, tempout, &ny);
    //fwrite(out,sizeof(gr_complex),ny,file2);
    // here "y" is signal with sample rate of OUT_SAMPLERATE and "ny" shows how many samples are in "y"
//}
    int i=0;
    int first=0;
    //dump_complex(out,ny,1);


    if (d_state==synchronized){
 sync:   	 for (;i<=ny;){
    			  //d_counter=0;
    			    //dump_complex(out,70000,0);
    		    	gsm_decode(tempout,70000);
    		    	//fprintf(stderr,"test1");
    		    	//consume_each(out,i)
    		    	//if(ctr>4)
    		    	//	exit(1);
    		    	i=i+625;
    		    	//ctr++;
    		    }
    }
    else{
    	 for (i=	70000;i<=ny;){
    			  //d_counter=0;
    			    //dump_complex(out,70000,0);
    		// fprintf(stderr,"\ntest3");
    		    	gsm_decode(tempout,70000);
    		    	//fprintf(stderr,"\ntest2");
    		    	//consume_each(out,i)
    		    	//if(ctr>4)
    		    	///	exit(1);
    		    	//printf("\ni=%d",i);
    		    	if (d_state==synchronized){
    		    		goto sync;
    		    	}
    		    	i=i+625;
    		    	//ctr++;

    		    }
    }



msresamp_crcf_destroy(q);

}


void  filter( gr_complex *in,int len){

	// options

	//fwrite(in,sizeof(gr_complex),len,file3);

	gr_complex *out;
	out = (gr_complex *)malloc(DEFAULT_BUF_LENGTH * sizeof(gr_complex));
	double samplerate=1e6;
	double cutofffreq=125e3;

	float ft=0.1f; //filter transition
	float As = 120.0f; // stop-band attenuation
	float mu=0.0f;//fractional timing offset
	// estimate required filter length and generate filter
	unsigned int h_len = estimate_req_filter_len(ft,As);
	float h[h_len];
	//liquid_firdes_kaiser(h_len,fc,As,mu,h);


	float fc = cutofffreq/samplerate; // cutoff frequency

	//printf("\nLen=%d\n",h_len);



	// design filter from prototype and scale to bandwidth
	firfilt_crcf q = firfilt_crcf_create_kaiser(h_len, fc, As, 0.0f);
	//firfilt_crcf_set_scale(q, 2.0f*cutofffreq);
	int i=0;


	for (i=0;i<len;i++){
		firfilt_crcf_push(q, in[i]);
		firfilt_crcf_execute(q, &out[i]);
	}

	resampler(out,len,samp_rate,4);

    // here lowpass filtered signal "rf" can be given to demodulator as shown above
}



int main(int argc, char **argv)
{
	gr_complex *buffer;
	int n_read;

	parse_options(argc,argv);

	buffer = (gr_complex*)malloc(out_block_size * sizeof(gr_complex));


	file = fopen(filename, "rb");

	if (!file){
		fprintf(stderr,"unable to open file");
		return EXIT_SUCCESS;
	}



	do{
		n_read=fread(buffer,sizeof(gr_complex),out_block_size*sizeof(gr_complex),file);
		if(n_read!=0){
			gr_complex out[n_read];
			memcpy(out, buffer, n_read);
			fprintf(stderr,"\nSending Samples1");

			filter( out,n_read);
			fprintf(stderr,"\nSending Samples2");

		}


	}while ( n_read >= 1) ;// expecting 1

	return EXIT_SUCCESS;
}
