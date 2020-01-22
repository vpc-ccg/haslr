/*****************************************************
 * Author: Ehsan Haghshenas (ehaghshe AT sfu DOT ca) *
 *****************************************************/

#include "Longread.hpp"
#include "Common.hpp"
#include <fstream>
#include <algorithm>
#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

// Longread_List_t *_read_lr_list;

void initialize_longread(Longread_List_t &lr_list)
{
    // 
    lr_list.reads_increment = 1000000;
    lr_list.reads_cap = lr_list.reads_increment;
    lr_list.reads = (Longread_t*) malloc(lr_list.reads_cap * sizeof(Longread_t));
    lr_list.reads_size = 0;
    // 
    lr_list.seqs_increment = 1000000;
    lr_list.seqs_cap = lr_list.seqs_increment;
    lr_list.seqs = (uint8_t*) malloc(lr_list.seqs_cap * sizeof(uint8_t));
    lr_list.seqs_size = 0;
    // 
    lr_list.alignments_increment = 1000000;
    lr_list.alignments_cap = lr_list.alignments_increment;
    lr_list.alignments = (Align_Seq_t*) malloc(lr_list.alignments_cap * sizeof(Align_Seq_t));
    lr_list.alignments_size = 0;
    // 
    lr_list.cigars_increment = 1000000;
    lr_list.cigars_cap = lr_list.alignments_increment;
    lr_list.cigars = (char*) malloc(lr_list.cigars_cap * sizeof(char));
    lr_list.cigars_size = 0;
}

void finalize_longread(Longread_List_t &lr_list)
{
    free(lr_list.reads);
    free(lr_list.seqs);
    free(lr_list.alignments);
    free(lr_list.cigars);
}

bool compare_Align_Seg(const Align_Seq_t &a, const Align_Seq_t &b)
{
    return (a.q_end < b.q_end) || (a.q_end == b.q_end && a.q_start < b.q_start);
}

bool compare_Align_Seg2(const Align_Seq2_t &a, const Align_Seq2_t &b)
{
    return (a.q_end < b.q_end) || (a.q_end == b.q_end && a.q_start < b.q_start);
}

void update_longreads(Longread_List_t &lr_list)
{
    uint64_t i, j;
    uint64_t off_seq = 0; // offset for seqs
    uint64_t off_aln = 0; // offset for alignments
    uint64_t off_cg = 0; // offset for cigars
    for(i = 0; i < lr_list.reads_size; i++)
    {
        // 
        lr_list.reads[i].comp_seq = lr_list.seqs + off_seq;
        off_seq += lr_list.reads[i].comp_len;
        // 
        if(lr_list.reads[i].contig_aln_num > 0)
        {
            lr_list.reads[i].contig_aln = lr_list.alignments + off_aln;
            off_aln += lr_list.reads[i].contig_aln_num;
            for(j = 0; j < lr_list.reads[i].contig_aln_num; j++)
            {
                lr_list.reads[i].contig_aln[j].cigar = lr_list.cigars + off_cg;
                off_cg += lr_list.reads[i].contig_aln[j].cigar_len + 1;
            }
        }
        else
        {
            lr_list.reads[i].contig_aln = NULL;
        }
    }
}

void load_longread_comments(string path, vector<string> &comm_list)
{
    gzFile fp = gzopen(path.c_str(), "r");
    if(fp == NULL)
    {
        fprintf(stderr, "[ERROR] (Longread::load_longread_comments) could not open file: %s\n", path.c_str());
        exit(EXIT_FAILURE);
    }
    // fprintf(stderr, "       processing file: %s...", path.c_str());
    // double cputime = get_cpu_time();
    // double realtime = get_real_time();
    kseq_t *seq = kseq_init(fp);
    int32_t cnt = 0;
    while (kseq_read(seq) >= 0)
    {
        comm_list[cnt++] = seq->comment.s;
    }
    kseq_destroy(seq);
    gzclose(fp);
    // 
    // fprintf(stderr, " Done in %.2lf CPU seconds (%.2lf real seconds)\n", get_cpu_time() - cputime, get_real_time() - realtime);
}

void load_longread_compressed(string path, Longread_List_t &lr_list)
{
    gzFile fp = gzopen(path.c_str(), "r");
    if(fp == NULL)
    {
        fprintf(stderr, "[ERROR] (DNA-Seq::load_seq_compressed) could not open file: %s\n", path.c_str());
        exit(EXIT_FAILURE);
    }
    fprintf(stderr, "       processing file: %s...", path.c_str());
    double cputime = get_cpu_time();
    double realtime = get_real_time();
    kseq_t *seq = kseq_init(fp);
    while (kseq_read(seq) >= 0)
    {
        // lr_list.push_back(Longread(seq->seq.s, seq->seq.l));
        uint32_t comp_len = (seq->seq.l / 4) + (seq->seq.l % 4 ? 1 : 0);
        // realloc reads if needed
        if(lr_list.reads_size >= lr_list.reads_cap)
        {
            lr_list.reads_cap += lr_list.reads_increment;
            lr_list.reads = (Longread_t*) realloc(lr_list.reads, lr_list.reads_cap * sizeof(Longread_t));
        }
        // realloc seqs if needed
        if(lr_list.seqs_size + comp_len > lr_list.seqs_cap)
        {
            lr_list.seqs_cap += lr_list.seqs_increment + comp_len;
            lr_list.seqs = (uint8_t*) realloc(lr_list.seqs, lr_list.seqs_cap * sizeof(uint8_t));
        }
        // 
        lr_list.reads[lr_list.reads_size].len = seq->seq.l;
        lr_list.reads[lr_list.reads_size].comp_len = comp_len;
        lr_list.reads[lr_list.reads_size].contig_aln_num = 0;
        // _ctg_reads[_ctg_reads_size].comp_seq = _ctg_block + _ctg_block_size;
        // lr_list.reads[lr_list.reads_size].comp_seq = NULL;
        // lr_list.reads[lr_list.reads_size].contig_aln = NULL;
        set_compressed_dna(seq->seq.s, seq->seq.l, lr_list.seqs + lr_list.seqs_size);
        lr_list.seqs_size += comp_len;
        lr_list.reads_size++;
    }
    kseq_destroy(seq);
    gzclose(fp);
    // 
    // update_longreads(lr_list);
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

void load_longread_compressed_fofn(string path, Longread_List_t &lr_list)
{
    ifstream fin(path.c_str());
    if(fin.is_open() == false)
    {
        fprintf(stderr, "[ERROR] (Longread::load_longread_compressed_fofn) could not open file: %s\n", path.c_str());
        exit(EXIT_FAILURE);
    }
    string line;
    while(getline(fin, line))
    {
        load_longread_compressed(line, lr_list);
    }
    fin.close();
    // 
    update_longreads(lr_list);
}

void process_lr_alignment_group(vector<Align_Seq2_t> &alns, Contig_List_t &contig_list, Longread_List_t &lr_list)
{
    if(alns.size() <= 1) return;

    // check for palindromic sequence first
    unordered_map<uint32_t, uint32_t> uniq_list;
    for(uint32_t i = 0; i < alns.size(); i++)
    {
        uint32_t tid = alns[i].t_id;
        if(contig_list.contigs[tid].mean_kmer < gopt.uniq_freq * (1 + gopt.max_uniq_dev))
        {
            if(uniq_list.count(tid) > 0)
            {
                alns.resize(i);
            }
            else
            {
                uniq_list[tid] = i;
            }
        }
    }

    for(uint32_t i = 0; i < alns.size(); i++)
    {
        // filter 5: alignments of contigs to the middle of the read that cover less than 90% of the whole contig
        if(i > 0 && i < alns.size() - 1 && (alns[i].t_end - alns[i].t_start) / (double)alns[i].t_len < 0.8) continue;
        // add alignment
        uint32_t cg_len = (uint32_t)alns[i].cigar.size();
        // realloc alignments if needed
        if(lr_list.alignments_size >= lr_list.alignments_cap)
        {
            lr_list.alignments_cap += lr_list.alignments_increment;
            lr_list.alignments = (Align_Seq_t*) realloc(lr_list.alignments, lr_list.alignments_cap * sizeof(Align_Seq_t));
        }
        lr_list.alignments[lr_list.alignments_size] = {alns[i].q_id, alns[i].q_start, alns[i].q_end, alns[i].t_id, alns[i].t_start, alns[i].t_end, 
            alns[i].n_match, alns[i].n_block, alns[i].is_rev, alns[i].mapq, 0, cg_len, NULL};
        lr_list.alignments_size++;
        // add cigar
        // realloc cigars if needed
        if(lr_list.cigars_size + cg_len + 1 > lr_list.cigars_cap)
        {
            lr_list.cigars_cap += lr_list.cigars_increment;
            lr_list.cigars = (char*) realloc(lr_list.cigars, lr_list.cigars_cap * sizeof(char));
        }
        strcpy(lr_list.cigars + lr_list.cigars_size, alns[i].cigar.c_str());
        lr_list.cigars_size += cg_len + 1;
        // 
        uint32_t rid = alns[i].q_id;
        lr_list.reads[rid].contig_aln_num++;
    }
}

void load_alignment(string path, Contig_List_t &contig_list, Longread_List_t &lr_list)
{
    ifstream fin(path.c_str());
    if(fin.is_open() == false)
    {
        fprintf(stderr, "[ERROR] (Longread::load_alignment) could not open file: %s\n", path.c_str());
        exit(EXIT_FAILURE);
    }
    // 
    fprintf(stderr, "       processing file: %s... ", path.c_str());
    double cputime = get_cpu_time();
    double realtime = get_real_time();
    string line;
    vector<string> fields;
    vector<Align_Seq2_t> lr_aln_list; // all alignment of a single long read
    string last_lr = "";
    while(getline(fin, line))
    {
        str_split(line, '\t', fields);
        if(last_lr != fields[0] && lr_aln_list.size() > 0) // process previous group
        {
            // sort by the q_end; required for weighted maximum interval scheduling
            sort(lr_aln_list.begin(), lr_aln_list.end(), compare_Align_Seg2);
            process_lr_alignment_group(lr_aln_list, contig_list, lr_list);
            lr_aln_list.clear();
        }

        // filter 1: alignment to small contigs are ignored
        if(str2type<uint32_t>(fields[10]) < gopt.min_aln_block) continue;
        
        // filter 2: alignment with low identity are ignored
        if(str2type<double>(fields[9]) / str2type<double>(fields[10]) < gopt.min_aln_sim) continue;

        // filter 3: alignments whose MAPQ is not high are ignored
        if(str2type<uint32_t>(fields[11]) < gopt.min_aln_mapq) continue;

        // filter 4: keep only alignments to unique contigs
        uint32_t tid = str2type<uint32_t>(fields[5]);
        if(contig_list.contigs[tid].mean_kmer > gopt.uniq_freq * (3 + gopt.max_uniq_dev)) continue;

        // find the cigar
        string cg = "";
        for(size_t i = 12; i < fields.size(); i++)
        {
            if(fields[i].substr(0, 5) == "cg:Z:")
            {
                cg = fields[i].substr(5);
                break;
            }
        }
        last_lr = fields[0];
        // add the new alignment
        lr_aln_list.push_back({str2type<uint32_t>(fields[0]), str2type<uint32_t>(fields[1]), str2type<uint32_t>(fields[2]), str2type<uint32_t>(fields[3]), 
            str2type<uint32_t>(fields[5]), str2type<uint32_t>(fields[6]), str2type<uint32_t>(fields[7]), str2type<uint32_t>(fields[8]), 
            str2type<uint32_t>(fields[9]), str2type<uint32_t>(fields[10]), (uint8_t)(fields[4][0] == '-' ? 1 : 0),
            (uint8_t)str2type<uint32_t>(fields[11]), cg});
        // NOTE: if I use str2type<uint8_t>(fields[11]) 60 is converted into 54! why?!
    }
    if(lr_aln_list.size() > 0) // process the last group
    {
        // sort by the q_end; required for weighted maximum interval scheduling
        sort(lr_aln_list.begin(), lr_aln_list.end(), compare_Align_Seg2);
        process_lr_alignment_group(lr_aln_list, contig_list, lr_list);
        lr_aln_list.clear();
    }
    fin.close();
    // update_longreads(lr_list);
    fprintf(stderr, "Done in %.2lf CPU seconds (%.2lf real seconds)\n", get_cpu_time() - cputime, get_real_time() - realtime);
}

void load_alignment_fofn(string path, Contig_List_t &contig_list, Longread_List_t &lr_list)
{
    ifstream fin(path.c_str());
    if(fin.is_open() == false)
    {
        fprintf(stderr, "[ERROR] (Longread::load_alignment_fofn) could not open file: %s\n", path.c_str());
        exit(EXIT_FAILURE);
    }
    string line;
    while(getline(fin, line))
    {
        load_alignment(line, contig_list, lr_list);
    }
    fin.close();
    // 
    update_longreads(lr_list);
}

void write_longread_index(string path, Longread_List_t &lr_list)
{
    FILE *fp = fopen(path.c_str(), "wb");
    if(fp == NULL)
    {
        fprintf(stderr, "[ERROR] (Longread::write_longread_index) could not open file: %s\n", path.c_str());
        exit(EXIT_FAILURE);
    }
    fwrite(&lr_list.reads_size, sizeof(uint64_t), 1, fp);
    fwrite(lr_list.reads, sizeof(Longread_t), lr_list.reads_size, fp);
    fwrite(&lr_list.seqs_size, sizeof(uint64_t), 1, fp);
    fwrite(lr_list.seqs, sizeof(uint8_t), lr_list.seqs_size, fp);
    fwrite(&lr_list.alignments_size, sizeof(uint64_t), 1, fp);
    fwrite(lr_list.alignments, sizeof(Align_Seq_t), lr_list.alignments_size, fp);
    fwrite(&lr_list.cigars_size, sizeof(uint64_t), 1, fp);
    fwrite(lr_list.cigars, sizeof(char), lr_list.cigars_size, fp);
    fclose(fp);
}

void read_longread_index(string path, Longread_List_t &lr_list)
{
    FILE *fp = fopen(path.c_str(), "rb");
    if(fp == NULL)
    {
        fprintf(stderr, "[ERROR] (Longread::read_longread_index) could not open file: %s\n", path.c_str());
        exit(EXIT_FAILURE);
    }
    fread(&lr_list.reads_size, sizeof(uint64_t), 1, fp);
    lr_list.reads = (Longread_t*) malloc(lr_list.reads_size * sizeof(Longread_t));
    fread(lr_list.reads, sizeof(Longread_t), lr_list.reads_size, fp);
    fread(&lr_list.seqs_size, sizeof(uint64_t), 1, fp);
    lr_list.seqs = (uint8_t*) malloc(lr_list.seqs_size * sizeof(uint8_t));
    fread(lr_list.seqs, sizeof(uint8_t), lr_list.seqs_size, fp);
    fread(&lr_list.alignments_size, sizeof(uint64_t), 1, fp);
    lr_list.alignments = (Align_Seq_t*) malloc(lr_list.alignments_size * sizeof(Align_Seq_t));
    fread(lr_list.alignments, sizeof(Align_Seq_t), lr_list.alignments_size, fp);
    fread(&lr_list.cigars_size, sizeof(uint64_t), 1, fp);
    lr_list.cigars = (char*) malloc(lr_list.cigars_size * sizeof(char));
    fread(lr_list.cigars, sizeof(char), lr_list.cigars_size, fp);
    fclose(fp);
    // 
    update_longreads(lr_list);
    // // uint64_t offset = 0;
    // for(uint32_t i = 0; i < _ctg_contigs_size; i++)
    // {
    //     fprintf(stdout, ">%u\n", i);
    //     string s2 = get_uncompressed_dna(_ctg_contigs[i].comp_seq, _ctg_contigs[i].len, _ctg_contigs[i].comp_len);
    //     fprintf(stdout, "%s\n", s2.c_str());
    //     // offset += _ctg_contigs[i].comp_len;
    // }
}

// lr_step and c_step: +1 / -1
void find_contig_pos(string &cigar_str, uint32_t &lr_curr, uint32_t &c_curr, int lr_step, int c_step, uint32_t lr_pos)
{
    // fprintf(stdout, "[find_contig_pos] cigarLen:%lu, lr_curr:%u, c_curr:%u, lr_step:%d, c_step:%d, lr_pos:%u\n", cigar_str.size(), lr_curr, c_curr, lr_step, c_step, lr_pos);
    uint32_t i;
    for(i = 0; i < cigar_str.size(); i++)
    {
        // fprintf(stdout, "[debug] %c %u %u\n", cigar_str[i], c_curr, lr_curr);
        if(lr_curr == lr_pos)
            break;
        if(cigar_str[i] == 'M')
        {
            c_curr += c_step;
            lr_curr += lr_step;
        }
        else if(cigar_str[i] == 'I')
        {
            lr_curr += lr_step;
        }
        else // if(cigar_str[i] == 'D')
        {
            c_curr += c_step;
        }
    }
    // fix the cigar if not ending at a match
    while(cigar_str[i] != 'M')
    {
        if(cigar_str[i-1] == 'M')
        {
            c_curr -= c_step;
            lr_curr -= lr_step;
        }
        else if(cigar_str[i-1] == 'I')
        {
            lr_curr -= lr_step;
        }
        else if(cigar_str[i-1] == 'D')
        {
            c_curr -= c_step;
        }
        i--;
    }
    // fprintf(stdout, "[debug] %u\n", i + 1);
    if(i + 1 < cigar_str.size())
        cigar_str.erase(i+1);
    // return c_curr;
}

uint32_t count_matches_expanded_cigar(string &exp_cigar)
{
    uint32_t n_match = 0;
    for(uint32_t i = 0; i < exp_cigar.size(); i++)
        n_match += (exp_cigar[i] == 'M');
    return n_match;
}

void fix_overlapping_alignments(Align_Seq_t *aln, int num)
{
    for(int i = 0; i < num - 1; i++)
    {
        if(aln[i].q_end > aln[i+1].q_start) // overlap between i and i+1
        {
            // TODO: give a warning if the overlap is above 50% of the shorter segment
            long long ovLen = (long long)aln[i].q_end - (long long)aln[i+1].q_start;
            // fprintf(stdout, "overlap\trid:%u\tcid_1:%u\tcid_2:%u\tlen:%lld\taln1_len:%u\taln2_len:%u\tmiddle_rpos:%lld\n", aln[i].q_id, aln[i].t_id, aln[i+1].t_id, ovLen, aln[i].n_block, aln[i+1].n_block, (long long)aln[i+1].q_start + (ovLen - ovLen/2));
            // fix first alignment
            string cigar_exp; // = expand_cigar(aln[i].cigar);
            // string cigar_exp_rev; // = reverseString(cigar_exp);
            uint32_t res_qpos;
            uint32_t res_tpos;
            string tmp_str;
            if(aln[i].is_rev == 0)
            {
                cigar_exp = expand_cigar(aln[i].cigar);
                res_qpos = aln[i].q_start;
                res_tpos = aln[i].t_start;
                find_contig_pos(cigar_exp, res_qpos, res_tpos, +1, +1, aln[i].q_end - ovLen/2 - 1);
                aln[i].q_end = res_qpos + 1;
                aln[i].t_end = res_tpos + 1;
                aln[i].n_block = cigar_exp.size();
                aln[i].n_match = count_matches_expanded_cigar(cigar_exp);
                tmp_str = collapse_cigar(cigar_exp);
                strcpy(aln[i].cigar, tmp_str.c_str());
            }
            else
            {
                cigar_exp = expand_cigar(aln[i].cigar);
                reverse(cigar_exp.begin(), cigar_exp.end());
                // cigar_exp_rev = reverseString(cigar_exp);
                res_qpos = aln[i].q_start;
                res_tpos = aln[i].t_end - 1;
                find_contig_pos(cigar_exp, res_qpos, res_tpos, +1, -1, aln[i].q_end - ovLen/2 - 1);
                aln[i].q_end = res_qpos + 1;
                aln[i].t_start = res_tpos;
                aln[i].n_block = cigar_exp.size();
                aln[i].n_match = count_matches_expanded_cigar(cigar_exp);
                reverse(cigar_exp.begin(), cigar_exp.end());
                tmp_str = collapse_cigar(cigar_exp);
                strcpy(aln[i].cigar, tmp_str.c_str());
            }
            // fprintf(stdout, "    1 query:%lld contig:%u\n", (long long)aln[i].q_end - ovLen/2 - 1, res_pos);
            // fix second alignment
            // cigar_exp = expand_cigar(aln[i+1].cigar);
            // cigar_exp_rev = reverseString(cigar_exp);
            if(aln[i+1].is_rev == 0)
            {
                cigar_exp = expand_cigar(aln[i+1].cigar);
                reverse(cigar_exp.begin(), cigar_exp.end());
                res_qpos = aln[i+1].q_end - 1;
                res_tpos = aln[i+1].t_end - 1;
                find_contig_pos(cigar_exp, res_qpos, res_tpos, -1, -1, (long long)aln[i+1].q_start + (ovLen - ovLen/2));
                aln[i+1].q_start = res_qpos;
                aln[i+1].t_start = res_tpos;
                aln[i+1].n_block = cigar_exp.size();
                aln[i+1].n_match = count_matches_expanded_cigar(cigar_exp);
                reverse(cigar_exp.begin(), cigar_exp.end());
                tmp_str = collapse_cigar(cigar_exp);
                strcpy(aln[i+1].cigar, tmp_str.c_str());
            }
            else
            {
                cigar_exp = expand_cigar(aln[i+1].cigar);
                res_qpos = aln[i+1].q_end - 1;
                res_tpos = aln[i+1].t_start;
                find_contig_pos(cigar_exp, res_qpos, res_tpos, -1, +1, (long long)aln[i+1].q_start + (ovLen - ovLen/2));
                aln[i+1].q_start = res_qpos;
                aln[i+1].t_end = res_tpos + 1;
                aln[i+1].n_block = cigar_exp.size();
                aln[i+1].n_match = count_matches_expanded_cigar(cigar_exp);
                tmp_str = collapse_cigar(cigar_exp);
                strcpy(aln[i+1].cigar, tmp_str.c_str());
            }
            // fprintf(stdout, "    ovlen:%lld\n", ovLen);
            // // print the update alignments
            // fprintf(stdout, "        %u-%u %c %u:%u-%u cg:%s\n", aln[i].q_start, aln[i].q_end - 1, (aln[i].is_rev ? '-' : '+'), aln[i].t_id, aln[i].t_start, aln[i].t_end - 1, aln[i].cigar);
            // fprintf(stdout, "        %u-%u %c %u:%u-%u cg:%s\n", aln[i+1].q_start, aln[i+1].q_end - 1, (aln[i+1].is_rev ? '-' : '+'), aln[i+1].t_id, aln[i+1].t_start, aln[i+1].t_end - 1, aln[i+1].cigar);
        }
    }
}

int latest_compatible(Align_Seq_t *aln_uniq[], int curr) 
{ 
    for (int j = curr - 1; j >= 0; j--) 
    { 
        if (aln_uniq[j]->q_end <= aln_uniq[curr]->q_start) 
            return j; 
    } 
    return -1; 
}

void find_best_scheduling(Align_Seq_t *aln, int num, vector<Align_Seq_t*> &compact_lr, Contig_List_t &contig_list, uint32_t min_aln_block, uint32_t copy_count)
{
    // fix_overlapping_alignments(aln, num);

    // dp[i] maximum sum of weights for first i intervals
    uint32_t dp[10000];
    Align_Seq_t* aln_uniq[10000];
    int num_uniq = 0;
    for(int i = 0; i < num; i++)
    {
        // filter 1: alignment to small contigs are ignored
        if(aln[i].n_block < min_aln_block) continue;

        // filter 4: keep only alignments to unique contigs
        uint32_t tid = aln[i].t_id;
        if(contig_list.contigs[tid].mean_kmer > gopt.uniq_freq * (copy_count + gopt.max_uniq_dev)) continue;
        // fprintf(stdout, "good_tid %u %lf\n", tid, contig_list.contigs[tid].mean_kmer);

        // // filter 1: alignment to small contigs are ignored
        // if(aln[i].n_block < gopt.minAlnBlock_uniq) continue;
        
        // // filter 2: alignment with low identity are ignored
        // if((double)aln[i].n_match / aln[i].n_block < gopt.minAlnSim_uniq) continue;

        // // filter 3: alignments whose MAPQ is not high are ignored
        // if(aln[i].mapq < gopt.minAlnMapq_uniq) continue;

        // // filter 4: keep only alignments to unique contigs
        // uint32_t tid = aln[i].t_id;
        // if(contig_list.contigs[tid].mean_kmer > gopt.uniqFreq * (1 + gopt.maxUniqDev)) continue;

        aln_uniq[num_uniq] = aln + i;
        num_uniq++;
    }

    if(num_uniq == 0) return;
    // fprintf(stdout, ">%u %d\n", rid, num_uniq);
    // for(int i = 0; i < num_uniq; i++)
    //     fprintf(stdout, "%u\n", aln_uniq[i]->t_id);

    vector<vector<int> > track;
    track.resize(num_uniq);
    dp[0] = aln_uniq[0]->n_match;
    track[0].push_back(0);

    // dynamic programming
    for (int i = 1; i < num_uniq; i++)
    {
        int j = latest_compatible(aln_uniq, i);
        // fprintf(stdout, "[debug] i:%d j:%d\n", i, j);
        if(j != -1)
        {
            if(aln_uniq[i]->n_match + dp[j] > dp[i-1])
            {
                dp[i] = aln_uniq[i]->n_match + dp[j];
                track[i] = track[j];
                track[i].push_back(i);
            }
            else
            {
                dp[i] = dp[i-1];
                track[i] = track[i-1];
            }
        }
        else
        {
            if(aln_uniq[i]->n_match > dp[i-1])
            {
                dp[i] = aln_uniq[i]->n_match;
                track[i].push_back(i);
            }
            else
            {
                dp[i] = dp[i-1];
                track[i] = track[i-1];
            }
        }
    }

    for(uint32_t i = 0; i < track[num_uniq - 1].size(); i++)
    {
        // if(i > 0 && i < track[num_uniq - 1].size() - 1 && (double)(aln_uniq[track[num_uniq - 1][i]]->t_end - aln_uniq[track[num_uniq - 1][i]]->t_start) / contig_list.contigs[aln_uniq[track[num_uniq - 1][i]]->t_id].len < 0.9)
        //     continue;
        compact_lr.push_back(aln_uniq[track[num_uniq - 1][i]]);
        // fprintf(stdout, "added %u %lf\n", compact_lr.back()->t_id, contig_list.contigs[compact_lr.back()->t_id].mean_kmer);
    }
}

void build_compact_longreads(Longread_List_t &lr_list, vector<vector<Align_Seq_t*> > &compact_lr_list_uniq, Contig_List_t &contig_list, uint32_t min_aln_block, uint32_t copy_count)
{
    compact_lr_list_uniq.resize(lr_list.reads_size);
    uint32_t i;
    for(i = 0; i < lr_list.reads_size; i++)
    {
        if(lr_list.reads[i].contig_aln_num > 0)
        {
            // sort(lr_list.reads[i].contig_aln, lr_list.reads[i].contig_aln + lr_list.reads[i].contig_aln_num, compare_Align_Seg);
            find_best_scheduling(lr_list.reads[i].contig_aln, lr_list.reads[i].contig_aln_num, compact_lr_list_uniq[i], contig_list, min_aln_block, copy_count);
        }
    }
}

void fix_alignments(Longread_List_t &lr_list)
{
    for(uint32_t i = 0; i < lr_list.reads_size; i++)
    {
        if(lr_list.reads[i].contig_aln_num > 0)
        {
            fix_overlapping_alignments(lr_list.reads[i].contig_aln, lr_list.reads[i].contig_aln_num);
        }
    }
}

// void* fix_alignments_ST(void *args_ptr)
// {
//     tuple<int, Longread_List_t*> args = *(tuple<int, Longread_List_t*> *)args_ptr;
//     int th_id = std::get<0>(args);
//     Longread_List_t *lr_list = std::get<1>(args);
//     // 
//     // fprintf(stdout, "%d\t%zu\t%zu\n", th_id, lr_list, lr_list->reads_size);
//     uint32_t div = lr_list->reads_size / gopt.num_threads;
//     div += (lr_list->reads_size % gopt.num_threads ? 1 : 0);
//     // for(uint32_t i = th_id; i < lr_list->reads_size; i += gopt.num_threads)
//     for(uint32_t i = th_id * div; i < (th_id + 1) * div; i++)
//     {
//         if(lr_list->reads[i].contig_aln_num > 0)
//         {
//             fix_overlapping_alignments(lr_list->reads[i].contig_aln, lr_list->reads[i].contig_aln_num);
//         }
//     }
//     return NULL;
// }

// void fix_alignments_MT(Longread_List_t &lr_list)
// {
//     // init arguments
//     vector<tuple<int, Longread_List_t*>> args_MT(gopt.num_threads);
//     for(uint32_t i = 0; i < gopt.num_threads; i++)
//     {
//         std::get<0>(args_MT[i]) = i;
//         std::get<1>(args_MT[i]) = &lr_list;
//     }
//     // spawn threads and join
//     pthread_t *sort_threads = (pthread_t*)malloc(gopt.num_threads * sizeof(pthread_t));
//     for(uint32_t i = 0; i < gopt.num_threads; i++)
//         pthread_create(sort_threads + i, NULL, fix_alignments_ST, &args_MT[i]);
//     for (uint32_t i = 0; i < gopt.num_threads; i++)
//         pthread_join(sort_threads[i], NULL);
//     free(sort_threads);
// }

void print_compact_longreads(vector<vector<Align_Seq_t*> > &compact_lr_list, string path)
{
    FILE *fp = fopen(path.c_str(), "w");
    if(fp == NULL)
    {
        fprintf(stderr, "[ERROR] cannot open file: %s\n", path.c_str());
        exit(EXIT_FAILURE);
    }
    // print compact long reads uniq
    for(uint32_t i=0; i < compact_lr_list.size(); i++)
    {
        fprintf(fp, ">%u\t", i);
        for(uint32_t j = 0; j < compact_lr_list[i].size(); j++)
            fprintf(fp, "%u-%u:%u:%c:%u-%u\t", compact_lr_list[i][j]->q_start, compact_lr_list[i][j]->q_end, compact_lr_list[i][j]->t_id, (compact_lr_list[i][j]->is_rev ? '-' : '+'), compact_lr_list[i][j]->t_start, compact_lr_list[i][j]->t_end);
        fprintf(fp, "\n");
    }
    // 
    fclose(fp);
}

void print_loaded_lrs(Longread_List_t &lr_list)
{
    for(uint32_t i = 0; i < lr_list.reads_size; i++)
    {
        fprintf(stdout, ">%u\n", i);
        string s2 = get_uncompressed_dna(lr_list.reads[i].comp_seq, lr_list.reads[i].len, lr_list.reads[i].comp_len);
        fprintf(stdout, "%s\n", s2.c_str());
    }
}

void print_loaded_alignments(Longread_List_t &lr_list, string logpath)
{
    FILE *fp = file_open_write(logpath.c_str());
    for(uint32_t i = 0; i < lr_list.reads_size; i++)
    {
        for(uint32_t j = 0; j < lr_list.reads[i].contig_aln_num; j++)
        {
            fprintf(fp, "%u\t%u\t%u\t%c\t%u\t%u\t%u\t%u\t%u\t%u\tcg:Z:%s\n", lr_list.reads[i].contig_aln[j].q_id, lr_list.reads[i].contig_aln[j].q_start, lr_list.reads[i].contig_aln[j].q_end,
                (lr_list.reads[i].contig_aln[j].is_rev ? '-' : '+'), lr_list.reads[i].contig_aln[j].t_id, lr_list.reads[i].contig_aln[j].t_start, 
                lr_list.reads[i].contig_aln[j].t_end, lr_list.reads[i].contig_aln[j].n_match, lr_list.reads[i].contig_aln[j].n_block, lr_list.reads[i].contig_aln[j].mapq, lr_list.reads[i].contig_aln[j].cigar);
        }
    }
    fclose(fp);
}