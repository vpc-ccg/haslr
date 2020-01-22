/*****************************************************
 * Author: Ehsan Haghshenas (ehaghshe AT sfu DOT ca) *
 *****************************************************/

#ifndef __COMPRESSED_SEQUENCE__
#define __COMPRESSED_SEQUENCE__

#include <string>
#include <stdint.h>

using namespace std;

void set_compressed_dna(char *src, int32_t len, uint8_t *dst);
string get_uncompressed_dna(uint8_t *comp_seq, int32_t len, int32_t comp_len);

#endif // __COMPRESSED_SEQUENCE__
