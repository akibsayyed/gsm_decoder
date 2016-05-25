/*
 * receiver.h
 *
 *  Created on: 23-May-2016
 *      Author: akib
 */

/*
 * receiver.h
 *
 *  Created on: 20-Jul-2015
 *      Author: akib
 */

#ifndef RECEIVER_H_
#define RECEIVER_H_
#include <gsm_constants.h>
#include <receiver_config.h>
#include <stdio.h>
#include <list>
      unsigned int d_c0_burst_start;
        float d_c0_signal_dbm;
        /**@name Configuration of the receiver */
        //@{
         int d_OSR=4; ///< oversampling ratio
         int d_chan_imp_length=CHAN_IMP_RESP_LENGTH; ///< channel impulse length
        float d_signal_dbm=-120;
        std::vector<int> d_tseq_nums; ///< stores training sequence numbers for channels different than C0
        std::vector<int> d_cell_allocation; ///< stores cell allocation - absolute rf channel numbers (ARFCNs) assigned to the given cell. The variable should at least contain C0 channel number.
        //@}

        gr_complex d_sch_training_seq[N_SYNC_BITS]; ///<encoded training sequence of a SCH burst
        gr_complex d_norm_training_seq[TRAIN_SEQ_NUM][N_TRAIN_BITS]; ///<encoded training sequences of a normal and dummy burst

        float d_last_time;

        /** Counts samples consumed by the receiver
         *
         * It is used in beetween find_fcch_burst and reach_sch_burst calls.
         * My intention was to synchronize this counter with some internal sample
         * counter of the USRP. Simple access to such USRP's counter isn't possible
         * so this variable isn't used in the "synchronized" state of the receiver yet.
         */
        unsigned d_counter=0;

        /**@name Variables used to store result of the find_fcch_burst fuction */
        //@{
        unsigned d_fcch_start_pos=0; ///< position of the first sample of the fcch burst
        float d_freq_offset_setting=0; ///< frequency offset set in frequency shifter located upstream
        //@}
        std::list<double> d_freq_offset_vals;

        /**@name Identifiers of the BTS extracted from the SCH burst */
        //@{
        int d_ncc; ///< network color code
        int d_bcc; ///< base station color code
        //@}

        /**@name Internal state of the gsm receiver */
        //@{
        enum states {
          fcch_search, sch_search, // synchronization search part
          synchronized // receiver is synchronized in this state
        } ;

        //@}

        /**@name Variables which make internal state in the "synchronized" state */
        //@{
        burst_counter d_burst_nr(4); ///< frame number and timeslot number
        channel_configuration d_channel_conf; ///< mapping of burst_counter to burst_type
        //@}

        unsigned d_failed_sch=0; ///< number of subsequent erroneous SCH bursts






#endif /* RECEIVER_H_ */
