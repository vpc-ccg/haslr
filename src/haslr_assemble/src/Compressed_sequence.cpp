/*****************************************************
 * Author: Ehsan Haghshenas (ehaghshe AT sfu DOT ca) *
 *****************************************************/

#include "Compressed_sequence.hpp"
#include "Common.hpp"
#include <fstream>
#include <vector>

uint8_t _dna_tableVal[128] = {
    4, 4, 4, 4,    4, 4, 4, 4,    4, 4, 4, 4,    4, 4, 4, 4, 
    4, 4, 4, 4,    4, 4, 4, 4,    4, 4, 4, 4,    4, 4, 4, 4, 
    4, 4, 4, 4,    4, 4, 4, 4,    4, 4, 4, 4,    4, 4, 4, 4, 
    4, 4, 4, 4,    4, 4, 4, 4,    4, 4, 4, 4,    4, 4, 4, 4, 
    4, 0, 4, 1,    4, 4, 4, 2,    4, 4, 4, 4,    4, 4, 4, 4, 
    4, 4, 4, 4,    3, 4, 4, 4,    4, 4, 4, 4,    4, 4, 4, 4, 
    4, 0, 4, 1,    4, 4, 4, 2,    4, 4, 4, 4,    4, 4, 4, 4, 
    4, 4, 4, 4,    3, 4, 4, 4,    4, 4, 4, 4,    4, 4, 4, 4
};

string get_uncompressed_dna(uint8_t *comp_seq, int32_t len, int32_t comp_len)
{
    // int32_t comp_len = (len / 4) + (len % 4 ? 1 : 0);
    string seq;
    int32_t i = comp_len - 1;
    int32_t j;
    if(len % 4)
    {
        for(j = 0; j < len % 4; j++)
        {
            seq += "ACGT"[((comp_seq[i] >> j*2) & 3)];
        }
        i--;
    }
    while(i >= 0)
    {
        for(j = 0; j < 4; j++)
        {
            seq += "ACGT"[((comp_seq[i] >> j*2) & 3)];
        }
        i--;
    }
    return seq;
}

void set_compressed_dna(char *src, int32_t len, uint8_t *dst)
{
    char *ptr_seq = src + len - 1;
    int32_t pos = 0;
    int32_t i = 0, j;
    uint8_t tmp;
    while(pos < len)
    {
        tmp = 0;
        for(j = 0; j < 4 && pos < len; j++)
        {
            tmp = ((tmp << 2) | (_dna_tableVal[(uint8_t)(*ptr_seq--)] & 3));
            pos++;
        }
        dst[i++] = tmp;
    }
}
