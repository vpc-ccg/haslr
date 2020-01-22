/*****************************************************
 * Author: Ehsan Haghshenas (ehaghshe AT sfu DOT ca) *
 *****************************************************/

#include "Contig.hpp"
#include "Common.hpp"
// #include <fstream>
#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

void initialize_contig(Contig_List_t &contig_list)
{
    // 
    contig_list.contigs_increment = 1000000;
    contig_list.contigs_cap = contig_list.contigs_increment;
    contig_list.contigs = (Contig_t*) malloc(contig_list.contigs_cap * sizeof(Contig_t));
    contig_list.contigs_size = 0;
    // 
    contig_list.block_increment = 1000000;
    contig_list.block_cap = contig_list.block_increment;
    contig_list.block = (uint8_t*) malloc(contig_list.block_cap * sizeof(uint8_t));
    contig_list.block_size = 0;
}

void finalize_contig(Contig_List_t &contig_list)
{
    free(contig_list.contigs);
    free(contig_list.block);
}

void update_contigs(Contig_List_t &contig_list)
{
    uint64_t i;
    uint64_t offset = 0;
    for(i = 0; i < contig_list.contigs_size; i++)
    {
        contig_list.contigs[i].comp_seq = contig_list.block + offset;
        offset += contig_list.contigs[i].comp_len;
    }
}

void load_contig_compressed(string path, Contig_List_t &contig_list)
{
    gzFile fp = gzopen(path.c_str(), "r");
    if(fp == NULL)
    {
        fprintf(stderr, "[ERROR] (Contig::load_contig_compressed) could not open file: %s\n", path.c_str());
        exit(EXIT_FAILURE);
    }
    // 
    fprintf(stderr, "       processing file: %s...", path.c_str());
    double cputime = get_cpu_time();
    double realtime = get_real_time();
    vector<string> comm_cols;
    kseq_t *seq = kseq_init(fp);
    while (kseq_read(seq) >= 0)
    {
        uint32_t comp_len = (seq->seq.l / 4) + (seq->seq.l % 4 ? 1 : 0);
        uint32_t kc = 0;
        double km = 0;
        char *p1;
        p1 = strstr(seq->comment.s, "KC:i:");
        kc = strtoul(p1+5, NULL, 10);
        p1 = strstr(seq->comment.s, "km:f:");
        km = strtod(p1+5, NULL);
        // km = atof(p1+5);
        // km = str2type<double>(p1+5);
        // fprintf(stdout, "%lu\t%lf\n", contig_list.contigs_size, km);

        // str_split(seq->comment.s, ' ', comm_cols);
        // for(uint32_t i = 0; i < comm_cols.size(); i++)
        // {
        //     if(comm_cols[i].substr(0, 5) == "km:f:")
        //     {
        //         km = str2type<double>(comm_cols[i].substr(5));
        //         // fprintf(stdout, "%lu\t%lf\n", contig_list.contigs_size, km);
        //         break;
        //     }
        // }
        // 
        if(contig_list.contigs_size >= contig_list.contigs_cap)
        {
            contig_list.contigs_cap += contig_list.contigs_increment;
            contig_list.contigs = (Contig_t*) realloc(contig_list.contigs, contig_list.contigs_cap * sizeof(Contig_t));
        }
        // 
        if(contig_list.block_size + comp_len > contig_list.block_cap)
        {
            contig_list.block_cap += contig_list.block_increment + comp_len;
            contig_list.block = (uint8_t*) realloc(contig_list.block, contig_list.block_cap * sizeof(uint8_t));
        }
        // 
        contig_list.contigs[contig_list.contigs_size].len = seq->seq.l;
        contig_list.contigs[contig_list.contigs_size].comp_len = comp_len;
        contig_list.contigs[contig_list.contigs_size].kmer_count = kc;
        contig_list.contigs[contig_list.contigs_size].mean_kmer = km;
        // _ctg_contigs[_ctg_contigs_size].comp_seq = _ctg_block + _ctg_block_size;
        set_compressed_dna(seq->seq.s, seq->seq.l, contig_list.block + contig_list.block_size);
        contig_list.block_size += comp_len;
        contig_list.contigs_size++;
    }
    kseq_destroy(seq);
    gzclose(fp);
    // 
    update_contigs(contig_list);
    fprintf(stderr, " Done in %.2lf CPU seconds (%.2lf real seconds)\n", get_cpu_time() - cputime, get_real_time() - realtime);
    // 
    // // uint64_t offset = 0;
    // for(uint32_t i = 0; i < _ctg_contigs_size; i++)
    // {
    //     fprintf(stdout, ">%u\n", i);
    //     string s2 = get_uncompressed_dna(_ctg_contigs[i].comp_seq, _ctg_contigs[i].len, _ctg_contigs[i].comp_len);
    //     fprintf(stdout, "%s\n", s2.c_str());
    //     // offset += _ctg_contigs[i].comp_len;
    // }
}

void write_contig_index(string path, Contig_List_t &contig_list)
{
    FILE *fp = fopen(path.c_str(), "wb");
    if(fp == NULL)
    {
        fprintf(stderr, "[ERROR] (Contig::write_contig_index) could not open file: %s\n", path.c_str());
        exit(EXIT_FAILURE);
    }
    fwrite(&contig_list.contigs_size, sizeof(uint64_t), 1, fp);
    fwrite(contig_list.contigs, sizeof(Contig_t), contig_list.contigs_size, fp);
    fwrite(&contig_list.block_size, sizeof(uint64_t), 1, fp);
    fwrite(contig_list.block, sizeof(uint8_t), contig_list.block_size, fp);
    fclose(fp);
}

void read_contig_index(string path, Contig_List_t &contig_list)
{
    FILE *fp = fopen(path.c_str(), "rb");
    if(fp == NULL)
    {
        fprintf(stderr, "[ERROR] (Contig::read_contig_index) could not open file: %s\n", path.c_str());
        exit(EXIT_FAILURE);
    }
    fread(&contig_list.contigs_size, sizeof(uint64_t), 1, fp);
    contig_list.contigs = (Contig_t*) malloc(contig_list.contigs_size * sizeof(Contig_t));
    fread(contig_list.contigs, sizeof(Contig_t), contig_list.contigs_size, fp);
    fread(&contig_list.block_size, sizeof(uint64_t), 1, fp);
    contig_list.block = (uint8_t*) malloc(contig_list.block_size * sizeof(uint8_t));
    fread(contig_list.block, sizeof(uint8_t), contig_list.block_size, fp);
    fclose(fp);
    // 
    update_contigs(contig_list);
    // // uint64_t offset = 0;
    // for(uint32_t i = 0; i < _ctg_contigs_size; i++)
    // {
    //     fprintf(stdout, ">%u\n", i);
    //     string s2 = get_uncompressed_dna(_ctg_contigs[i].comp_seq, _ctg_contigs[i].len, _ctg_contigs[i].comp_len);
    //     fprintf(stdout, "%s\n", s2.c_str());
    //     // offset += _ctg_contigs[i].comp_len;
    // }
}

// estimates the kmer frequency of unique regions
void calc_uniq_freq(Contig_List_t &contig_list)
{
    uint32_t i;
    vector<pair<uint32_t, double> > contig_freq(contig_list.contigs_size);
    for(i = 0; i < contig_list.contigs_size; i++)
        contig_freq[i] = {contig_list.contigs[i].len, contig_list.contigs[i].mean_kmer};
    sort(contig_freq.begin(), contig_freq.end(), greater<pair<uint32_t, double> >());
    double freq = 0;
    // TODO: what if there are only very few contigs?
    for(i = 0; i < 20 && i < contig_freq.size(); i++)
        freq += contig_freq[i].second;
    gopt.uniq_freq = freq / i;
}

void print_loaded_contigs(Contig_List_t &contig_list)
{
    for(uint32_t i = 0; i < contig_list.contigs_size; i++)
    {
        fprintf(stdout, ">%u\n", i);
        string s2 = get_uncompressed_dna(contig_list.contigs[i].comp_seq, contig_list.contigs[i].len, contig_list.contigs[i].comp_len);
        fprintf(stdout, "%s\n", s2.c_str());
    }
}
