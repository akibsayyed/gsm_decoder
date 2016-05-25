

#ifndef __SCH_H__
#define __SCH_H__ 1
#include <gsm_constants.h>
//#include <grgsm/api.h>
#ifdef gnuradio_gsm_EXPORTS
#  define GSM_API __GR_ATTR_EXPORT
#else
#  define GSM_API __GR_ATTR_IMPORT
#endif

#ifdef __cplusplus
extern "C"
{
#endif

   int decode_sch(const unsigned char *buf, int * t1_o, int * t2_o, int * t3_o, int * ncc, int * bcc);

#ifdef __cplusplus
}
#endif

#endif
