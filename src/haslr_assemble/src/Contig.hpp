/*****************************************************
 * Author: Ehsan Haghshenas (ehaghshe AT sfu DOT ca) *
 *****************************************************/

#ifndef __CONTIG__
#define __CONTIG__

#include <string>
#include <vector>
#include "Compressed_sequence.hpp"

using namespace std;

typedef struct
{
	uint32_t len;
	uint32_t comp_len;
	uint32_t kmer_count;
	double mean_kmer;
    uint8_t *comp_seq;
} Contig_t;

typedef struct
{
	Contig_t *contigs;
	uint64_t  contigs_size;
	uint64_t  contigs_cap;
	uint64_t  contigs_increment;
	uint8_t  *block;
	uint64_t  block_size;
	uint64_t  block_cap;
	uint64_t  block_increment;
} Contig_List_t;

void initialize_contig(Contig_List_t &contig_list);
void finalize_contig(Contig_List_t &contig_list);
void load_contig_compressed(string path, Contig_List_t &contig_list);
void write_contig_index(string path, Contig_List_t &contig_list);
void read_contig_index(string path, Contig_List_t &contig_list);
void calc_uniq_freq(Contig_List_t &contig_list);
void print_loaded_contigs(Contig_List_t &contig_list);

#endif // __CONTIG__
