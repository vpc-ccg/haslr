// /*****************************************************
//  * Author: Ehsan Haghshenas (ehaghshe AT sfu DOT ca) *
//  *****************************************************/

#include "Assemble.hpp"
#include "spoa.hpp"

#define SPOA_ALGORITHM 1
#define SPOA_SCORE_MATCH 5
#define SPOA_SCORE_MIS -4
#define SPOA_SCORE_GAP -8

vector<BBG_Node_t> *_asm_bbgraph;
Contig_List_t *_asm_contig_list;
Longread_List_t *_asm_lr_list;
vector<vector<Align_Seq_t*>> *_asm_compact_lr_list;
// 
pthread_mutex_t _asm_coord_lock;
uint32_t _asm_vertex_id;
uint32_t _asm_vertex_max_id;
map<uint32_t, BBG_Edge_t>::iterator _asm_it;
uint32_t _asm_flag_visited;

void asm_best_supported_interval_contig1(vector<pair<uint32_t, uint32_t>> &beg_list, vector<pair<uint32_t, uint32_t>> &end_list, pair<uint32_t, uint32_t> &best_int, set<uint32_t> &best_lrs)
{
    sort(beg_list.begin(), beg_list.end());
    sort(end_list.begin(), end_list.end());
    // 
    int curr_supp = 0;
    int best_supp = 0;
    int i = 0;
    int j = 0;
    int len = beg_list.size();
    uint32_t beg_best;
    uint32_t end_best;
    bool interval_started = false;
    set<uint32_t> curr_lrs;
    // 
    while(i < len && j < len)
    {
        if(beg_list[i].first < end_list[j].first)
        {
            curr_supp++;
            curr_lrs.insert(beg_list[i].second);
            if(curr_supp >= best_supp)
            {
                best_supp = curr_supp;
                beg_best = beg_list[i].first;
                best_lrs = curr_lrs;
                interval_started = true;
            }
            i++;
        }
        else
        {
            if(interval_started == true)
            {
                end_best = end_list[j].first;
                interval_started = false;
            }
            curr_supp--;
            curr_lrs.erase(end_list[j].second);
            j++;
        }
    }
    if(interval_started == true)
    {
        end_best = end_list[j].first;
    }
    // fprintf(stdout, "[best interval] i:%d j:%d\n", i, j);
    // for(set<uint32_t>::iterator it = best_lrs.begin(); it != best_lrs.end(); it++)
    //     fprintf(stdout, "    lr: %u\n", *it);
    best_int = {beg_best, end_best};
}

void asm_best_supported_interval_contig2(vector<pair<uint32_t, uint32_t> > &beg_list, vector<pair<uint32_t, uint32_t> > &end_list, pair<uint32_t, uint32_t> &best_int, set<uint32_t> &best_lrs)
{
    sort(beg_list.begin(), beg_list.end());
    sort(end_list.begin(), end_list.end());
    // 
    int curr_supp = 0;
    int best_supp = 0;
    int i = 0;
    int j = 0;
    int len = beg_list.size();
    uint32_t beg_best;
    uint32_t end_best;
    bool interval_started = false;
    set<uint32_t> curr_lrs;
    // 
    while(i < len && j < len)
    {
        if(beg_list[i].first < end_list[j].first)
        {
            curr_supp++;
            curr_lrs.insert(beg_list[i].second);
            if(curr_supp > best_supp)
            {
                best_supp = curr_supp;
                beg_best = beg_list[i].first;
                best_lrs = curr_lrs;
                interval_started = true;
            }
            i++;
        }
        else
        {
            if(interval_started == true)
            {
                end_best = end_list[j].first;
                interval_started = false;
            }
            curr_supp--;
            curr_lrs.erase(end_list[j].second);
            j++;
        }
    }
    if(interval_started == true)
    {
        end_best = end_list[j].first;
    }
    // fprintf(stdout, "[best interval] i:%d j:%d\n", i, j);
    // for(set<uint32_t>::iterator it = best_lrs.begin(); it != best_lrs.end(); it++)
    //     fprintf(stdout, "    lr: %u\n", *it);
    best_int = {beg_best, end_best};
}

// lr_step and c_step: +1 / -1
long long asm_find_lr_pos(string cigar_str, uint32_t lr_curr, uint32_t c_curr, int lr_step, int c_step, uint32_t contig_pos)
{
    // fprintf(stdout, "[find_lr_pos] cigarLen:%zu, lr_curr:%u, c_curr:%u, lr_step:%d, c_step:%d, contig_pos:%u\n", cigar_str.size(), lr_curr, c_curr, lr_step, c_step, contig_pos);
    if((c_step > 0 && c_curr > contig_pos) || (c_step < 0 && c_curr < contig_pos))
        return -1;
    for(uint32_t i = 0; i < cigar_str.size(); i++)
    {
        // fprintf(stdout, "[debug] %c %u %u\n", cigar_str[i], c_curr, lr_curr);
        if(c_curr == contig_pos)
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
    // fprintf(stdout, "[find_lr_pos] lr_curr:%u c_curr:%u\n", lr_curr, c_curr);
    return lr_curr;
}

void asm_calc_single_edge_coordinates(uint32_t node1, uint32_t rev1, uint32_t node2, uint32_t rev2, FILE *fp_log)
{
    uint32_t i;
    // 
    uint32_t to_node1 = (node2 << 1) | rev2;
    if((*_asm_bbgraph)[node1].edges[rev1].count(to_node1) == 0)
    {
        fprintf(stderr, "[ERROR] (Assemble::asm_calc_single_edge_coordinates) outgoing edge not found! %u:%c ===> %u:%c\n", node1, "+-"[rev1], node2, "+-"[rev2]);
        return;
    }
    BBG_Edge_t &edge1 = (*_asm_bbgraph)[node1].edges[rev1][to_node1];
    // 
    uint32_t to_node2 = (node1 << 1) | (1 - rev1);
    if((*_asm_bbgraph)[node2].edges[1-rev2].count(to_node2) == 0)
    {
        fprintf(stderr, "[ERROR] (Assemble::asm_calc_single_edge_coordinates) outgoing edge not found! %u:%c ===> %u:%c\n", node2, "+-"[1-rev2], node1, "+-"[1-rev1]);
        return;
    }
    BBG_Edge_t &edge2 = (*_asm_bbgraph)[node2].edges[1-rev2][to_node2];
    fprintf(fp_log, "edge      %u:%c -> %u:%c\n", node1, "+-"[rev1], node2, "+-"[rev2]);
    fprintf(fp_log, "edge_twin %u:%c -> %u:%c\n", node2, "+-"[1-rev2], node1, "+-"[1-rev1]);
    // fflush(fp_log);
    //
    vector<Edge_Supp_t> &edge_supp = edge1.edge_supp;
    fprintf(fp_log, "\tedge_supp size:%zu\n", edge_supp.size());
    // fflush(fp_log);
    for(i = 0; i < edge_supp.size(); i++)
    {
        uint32_t rid = edge_supp[i].lr_id;
        uint32_t cmpId1 = edge_supp[i].cmp_head_id;
        uint32_t cmpId2 = edge_supp[i].cmp_tail_id;
        // fprintf(fp_log, "\t\tsupp_lr lr_id:%u lr_size:%zu cmpId1:%u cmpId2:%u\n", rid, (*_cl_compact_lr_list_uniq)[rid].size(), cmpId1, cmpId2);
        // fflush(fp_log);
        fprintf(fp_log, "\tsupp_detail head\t%u\t%u\t%c\ttail\t%u\t%u\t%c\n", 
            (*_asm_compact_lr_list)[rid][cmpId1]->t_start, (*_asm_compact_lr_list)[rid][cmpId1]->t_end, "+-"[(*_asm_compact_lr_list)[rid][cmpId1]->is_rev], 
            (*_asm_compact_lr_list)[rid][cmpId2]->t_start, (*_asm_compact_lr_list)[rid][cmpId2]->t_end, "+-"[(*_asm_compact_lr_list)[rid][cmpId2]->is_rev]);
        // fflush(fp_log);
    }
    // find the best supported interval for contig1
    vector<pair<uint32_t, uint32_t>> beg_list1;
    vector<pair<uint32_t, uint32_t>> end_list1;
    pair<uint32_t, uint32_t> best_int1;
    set<uint32_t> best_lrs1;
    for(i = 0; i < edge_supp.size(); i++)
    {
        uint32_t rid = edge_supp[i].lr_id;
        uint32_t cmpId1 = edge_supp[i].cmp_head_id;
        beg_list1.push_back({(*_asm_compact_lr_list)[rid][cmpId1]->t_start, i});
        end_list1.push_back({(*_asm_compact_lr_list)[rid][cmpId1]->t_end, i});
    }
    asm_best_supported_interval_contig1(beg_list1, end_list1, best_int1, best_lrs1);
    fprintf(fp_log, "    @@@ best interval contig1 %u %u\n", best_int1.first, best_int1.second);
    // fflush(fp_log);
    // find the best supported interval for contig2
    vector<pair<uint32_t, uint32_t> > beg_list2;
    vector<pair<uint32_t, uint32_t> > end_list2;
    pair<uint32_t, uint32_t> best_int2;
    set<uint32_t> best_lrs2;
    for(i = 0; i < edge_supp.size(); i++)
    {
        uint32_t rid = edge_supp[i].lr_id;
        uint32_t cmpId2 = edge_supp[i].cmp_tail_id;
        beg_list2.push_back({(*_asm_compact_lr_list)[rid][cmpId2]->t_start, i});
        end_list2.push_back({(*_asm_compact_lr_list)[rid][cmpId2]->t_end, i});
    }
    asm_best_supported_interval_contig2(beg_list2, end_list2, best_int2, best_lrs2);
    fprintf(fp_log, "    @@@ best_interval contig2 %u %u\n", best_int2.first, best_int2.second);
    // fflush(fp_log);
    // 
    uint32_t contig1_pos; // end position of shared region of contig1 (inclusive)
    uint32_t contig2_pos; // start position of shared region of contig2 (inclusive)
    if(rev1 == 0)
        contig1_pos = best_int1.second - 1;
    else
        contig1_pos = best_int1.first;
    if(rev2 == 0)
        contig2_pos = best_int2.first;
    else
        contig2_pos = best_int2.second - 1;
    // get the list of long reads that support the best supported interval
    vector<uint32_t> best_lrs;
    std::set_intersection(best_lrs1.begin(), best_lrs1.end(), best_lrs2.begin(), best_lrs2.end(), std::inserter(best_lrs, best_lrs.begin()));

    fprintf(fp_log, "coordinates contig1_pos: %u\tcontig2_pos: %u\n", contig1_pos, contig2_pos);
    fprintf(fp_log, "supproting_lr: %lu\n", best_lrs.size());
    // fflush(fp_log);

    if(best_lrs.size() == 0)
    {
        edge1.cns_supp.clear();
        edge2.cns_supp.clear();
        edge1.head_end = edge2.tail_beg = (rev1 == 0 ? _asm_contig_list->contigs[node1].len - 1 : 0);
        edge1.tail_beg = edge2.head_end = (rev2 == 0 ? 0 : _asm_contig_list->contigs[node2].len - 1);
        return;
    }

    for(i = 0; i < best_lrs.size(); i++)
    {
        uint32_t rid = edge_supp[best_lrs[i]].lr_id;
        uint32_t rlen = _asm_lr_list->reads[rid].len;
        uint32_t cmpId1 = edge_supp[best_lrs[i]].cmp_head_id;
        uint32_t cmpId2 = edge_supp[best_lrs[i]].cmp_tail_id;
        Align_Seq_t *aln_cmp1 = (*_asm_compact_lr_list)[rid][cmpId1];
        Align_Seq_t *aln_cmp2 = (*_asm_compact_lr_list)[rid][cmpId2];
        // 
        uint32_t rstrand = (rev1 == aln_cmp1->is_rev ? 0 : 1);
        fprintf(fp_log, "    +++ lr:%u len:%u strand:%c\n", rid, _asm_lr_list->reads[rid].len, "+-"[rstrand]);
        // fflush(fp_log);
        long long lr_start = -1; // exlusive
        long long lr_end = -1; // exlusive
        string cigar_exp;
        string cigar_exp_rev;
        if(rstrand == 0)
        {
            cigar_exp = expand_cigar(aln_cmp1->cigar);
            cigar_exp_rev = reverseString(cigar_exp);
            if(rev1 == 0)
            {
                fprintf(fp_log, "        case 1\n");
                lr_start = asm_find_lr_pos(cigar_exp, aln_cmp1->q_start, aln_cmp1->t_start, +1, +1, contig1_pos);
            }
            else
            {
                fprintf(fp_log, "        case 2\n");
                lr_start = asm_find_lr_pos(cigar_exp_rev, aln_cmp1->q_start, aln_cmp1->t_end - 1, +1, -1, contig1_pos);
            }
            // 
            cigar_exp = expand_cigar(aln_cmp2->cigar);
            cigar_exp_rev = reverseString(cigar_exp);
            if(rev2 == 0)
            {
                fprintf(fp_log, "        case 3\n");
                lr_end = asm_find_lr_pos(cigar_exp_rev, aln_cmp2->q_end - 1, aln_cmp2->t_end - 1, -1, -1, contig2_pos);
            }
            else
            {
                fprintf(fp_log, "        case 4\n");
                lr_end = asm_find_lr_pos(cigar_exp, aln_cmp2->q_end - 1, aln_cmp2->t_start, -1, +1, contig2_pos);
            }
        }
        else // rstrand == 1
        {
            cigar_exp = expand_cigar(aln_cmp1->cigar);
            cigar_exp_rev = reverseString(cigar_exp);
            if(rev1 == 0)
            {
                fprintf(fp_log, "        case 5\n");
                lr_start = asm_find_lr_pos(cigar_exp, _asm_lr_list->reads[rid].len - aln_cmp1->q_end, aln_cmp1->t_start, +1, +1, contig1_pos);
            }
            else
            {
                fprintf(fp_log, "        case 6\n");
                lr_start = asm_find_lr_pos(cigar_exp_rev, _asm_lr_list->reads[rid].len - aln_cmp1->q_end, aln_cmp1->t_end - 1 , +1, -1, contig1_pos);
            }
            // 
            cigar_exp = expand_cigar(aln_cmp2->cigar);
            cigar_exp_rev = reverseString(cigar_exp);
            if(rev2 == 0)
            {
                fprintf(fp_log, "        case 7\n");
                lr_end = asm_find_lr_pos(cigar_exp_rev, _asm_lr_list->reads[rid].len - aln_cmp2->q_start - 1, aln_cmp2->t_end - 1, -1, -1, contig2_pos);
            }
            else
            {
                fprintf(fp_log, "        case 8\n");
                lr_end = asm_find_lr_pos(cigar_exp, _asm_lr_list->reads[rid].len - aln_cmp2->q_start - 1, aln_cmp2->t_start, -1, +1, contig2_pos);
            }
        }
        // 
        if(lr_start != -1 && lr_end != -1)
        {
            fprintf(fp_log, "        [coordinate] subseq_len:%lld lr_start:%lld lr_end:%lld\n", lr_end - lr_start - 1, lr_start + 1, lr_end - 1);
            // fflush(fp_log);
            edge1.cns_supp.push_back({rid, rstrand, uint32_t(lr_start + 1), uint32_t(lr_end - 1)});
            edge2.cns_supp.push_back({rid, 1 - rstrand, uint32_t(rlen - (lr_end - 1) - 1), uint32_t(rlen - (lr_start + 1) - 1)});
        }
        else
        {
            fprintf(fp_log, "        [coordinate] could not extract subseq\n");
            // fflush(fp_log);
        }
    } // END for(i = 0; i < best_lrs.size(); i++)
    if(edge1.cns_supp.size() > 0)
    {
        // int64_t sum_len = 0;
        // for(i = 0; i < cns_info.lr_seqs.size(); i++)
        // {
        //     // fprintf(fp_log, "[debug2] rid:%u rstrand:%u spos:%u epos:%u\n", cns_info.lr_seqs[i].lr_id, cns_info.lr_seqs[i].lr_strand, cns_info.lr_seqs[i].spos, cns_info.lr_seqs[i].epos);
        //     sum_len += (cns_info.lr_seqs[i].epos - cns_info.lr_seqs[i].spos + 1);
        //     // fprintf(fp_log, "[debug3] sum_len:%lld\n", sum_len);
        // }
        // cns_info.len_estimate = (int32_t)ceil((double)sum_len/cns_info.lr_seqs.size());
        // fprintf(fp_log, "        [debug] len_estimate: %d\n", cns_info.len_estimate);
        // fflush(fp_log);
        edge1.head_end = edge2.tail_beg = contig1_pos;
        edge1.tail_beg = edge2.head_end = contig2_pos;
    }
    else
    {
        edge1.cns_supp.clear();
        edge2.cns_supp.clear();
        // cns_info.len_estimate = -1;
        edge1.head_end = edge2.tail_beg = (rev1 == 0 ? _asm_contig_list->contigs[node1].len - 1 : 0);
        edge1.tail_beg = edge2.head_end = (rev2 == 0 ? 0 : _asm_contig_list->contigs[node2].len - 1);
    }
    fprintf(fp_log, "\n");
}

bool asm_init_vertex_index(int flag)
{
    _asm_flag_visited = flag;
    _asm_vertex_max_id = _asm_bbgraph->size() * 2;
    _asm_vertex_id = 0;
    for( ; _asm_vertex_id < _asm_vertex_max_id; _asm_vertex_id++)
    {
        uint32_t node_id = _asm_vertex_id / 2;
        uint32_t node_rev = _asm_vertex_id % 2;
        for(map<uint32_t, BBG_Edge_t>::iterator it = (*_asm_bbgraph)[node_id].edges[node_rev].begin(); it != (*_asm_bbgraph)[node_id].edges[node_rev].end(); it++)
        {
            if(it->second.flag != _asm_flag_visited)
            {
                _asm_it = it;
                return true;
            }
        }
    }
    return false;
}

bool asm_get_next_edge(pair<uint32_t, map<uint32_t, BBG_Edge_t>::iterator> &p1)
{
    bool ret_val;
    pthread_mutex_lock(&_asm_coord_lock);
    if(_asm_vertex_id < _asm_vertex_max_id)
    {
        // prepare return
        ret_val = true;
        p1.first = _asm_vertex_id;
        p1.second = _asm_it;
        // set flag
        uint32_t node1 = _asm_vertex_id / 2;
        uint32_t rev1 = _asm_vertex_id % 2;
        uint32_t node2 = _asm_it->first >> 1;
        uint32_t rev2 = _asm_it->first & 1;
        uint32_t to_node1 = (node2 << 1) | rev2;
        uint32_t to_node2 = (node1 << 1) | (1 - rev1);
        (*_asm_bbgraph)[node1].edges[rev1][to_node1].flag = _asm_flag_visited;
        (*_asm_bbgraph)[node2].edges[1-rev2][to_node2].flag = _asm_flag_visited;
        // 
        _asm_it++;
        while(_asm_it != (*_asm_bbgraph)[node1].edges[rev1].end() && _asm_it->second.flag == _asm_flag_visited) _asm_it++;
        if(_asm_it == (*_asm_bbgraph)[node1].edges[rev1].end()) // no other unvisited edge found on this vertex (i.e. node+strand)
        {
            for(_asm_vertex_id = _asm_vertex_id + 1; _asm_vertex_id < _asm_vertex_max_id; _asm_vertex_id++)
            {
                node1 = _asm_vertex_id / 2;
                rev1 = _asm_vertex_id % 2;
                bool found = false;
                for(map<uint32_t, BBG_Edge_t>::iterator it = (*_asm_bbgraph)[node1].edges[rev1].begin(); it != (*_asm_bbgraph)[node1].edges[rev1].end(); it++)
                {
                    if(it->second.flag != _asm_flag_visited)
                    {
                        _asm_it = it;
                        found = true;
                        break;
                    }
                }
                if(found) break;
            }
        }
    }
    else
    {
        ret_val = false;
    }
    pthread_mutex_unlock(&_asm_coord_lock);
    return ret_val;
}

void* asm_calc_edge_coordinates_ST(void *args_ptr)
{
    tuple<int, FILE*> args = *(tuple<int, FILE*> *)args_ptr;
    int th_id = std::get<0>(args);
    FILE *fp_log = std::get<1>(args);
    pair<uint32_t, map<uint32_t, BBG_Edge_t>::iterator> p1;
    while(asm_get_next_edge(p1))
    {
        uint32_t node_id = p1.first / 2;
        uint32_t node_rev = p1.first % 2;
        fprintf(fp_log, "calc_coords th_id:%d %u:%c -> %u:%c\n", th_id, node_id, "+-"[node_rev], p1.second->first >> 1,  "+-"[p1.second->first & 1]);
        fflush(fp_log);
        asm_calc_single_edge_coordinates(node_id, node_rev, p1.second->first >> 1, p1.second->first & 1, fp_log);
    }
    return NULL;
}

void asm_calc_edge_coordinates_MT(vector<BBG_Node_t> *bbgraph, Contig_List_t *contig_list, Longread_List_t *lr_list, vector<vector<Align_Seq_t*>> *compact_lr_list, string log_path)
{
    _asm_bbgraph = bbgraph;
    _asm_contig_list = contig_list;
    _asm_lr_list = lr_list;
    _asm_compact_lr_list = compact_lr_list;
    // 
    if(asm_init_vertex_index(11) == false) return;
    // init arguments
    FILE *fp_log = file_open_write(log_path.c_str());
    vector<tuple<int, FILE*>> args_MT(gopt.num_threads);
    for(uint32_t i = 0; i < gopt.num_threads; i++)
    {
        std::get<0>(args_MT[i]) = i;
        std::get<1>(args_MT[i]) = fp_log;
    }
    // spawn threads and join
    pthread_t *coord_threads = (pthread_t*)malloc(gopt.num_threads * sizeof(pthread_t));
    for(uint32_t i = 0; i < gopt.num_threads; i++)
        pthread_create(coord_threads + i, NULL, asm_calc_edge_coordinates_ST, &args_MT[i]);
    // 
    for (uint32_t i = 0; i < gopt.num_threads; i++)
        pthread_join(coord_threads[i], NULL);
    free(coord_threads);
}

void asm_calc_single_cns_seq(uint32_t node1, uint32_t rev1, uint32_t node2, uint32_t rev2, FILE* fp_log)
{
    uint32_t i;
    // 
    uint32_t to_node1 = (node2 << 1) | rev2;
    if((*_asm_bbgraph)[node1].edges[rev1].count(to_node1) == 0)
    {
        fprintf(stderr, "[ERROR] (Assemble::asm_calc_single_cns_seq) outgoing edge not found! %u:%c ===> %u:%c\n", node1, "+-"[rev1], node2, "+-"[rev2]);
        return;
    }
    BBG_Edge_t &edge1 = (*_asm_bbgraph)[node1].edges[rev1][to_node1];
    // 
    uint32_t to_node2 = (node1 << 1) | (1 - rev1);
    if((*_asm_bbgraph)[node2].edges[1-rev2].count(to_node2) == 0)
    {
        fprintf(stderr, "[ERROR] (Assemble::asm_calc_single_cns_seq) outgoing edge not found! %u:%c ===> %u:%c\n", node2, "+-"[1-rev2], node1, "+-"[1-rev1]);
        return;
    }
    BBG_Edge_t &edge2 = (*_asm_bbgraph)[node2].edges[1-rev2][to_node2];
    // generate the consensus sequence
    auto alignment_engine = spoa::createAlignmentEngine(static_cast<spoa::AlignmentType>(SPOA_ALGORITHM), SPOA_SCORE_MATCH, SPOA_SCORE_MIS, SPOA_SCORE_GAP);
    auto graph_cns = spoa::createGraph();
    fprintf(fp_log, "[shared_region] head_end:%u\ttail_beg:%u\n", edge1.head_end, edge1.tail_beg);
    uint32_t cnt_non_empty = 0;
    for(i = 0; i < edge1.cns_supp.size(); i++)
    {
        uint32_t rid = edge1.cns_supp[i].lr_id;
        uint32_t rstrand = edge1.cns_supp[i].lr_strand;
        string rseq = get_uncompressed_dna(_asm_lr_list->reads[rid].comp_seq, _asm_lr_list->reads[rid].len, _asm_lr_list->reads[rid].comp_len);
        string rseq_rev = reverseComplement(rseq);
        fprintf(fp_log, "        [debug] lr_id:%u lr_len:%u region_start:%u region_end:%u subseq_len:%u\n",
            rid, _asm_lr_list->reads[rid].len, edge1.cns_supp[i].spos, edge1.cns_supp[i].epos, edge1.cns_supp[i].epos - edge1.cns_supp[i].spos + 1);
        // fflush(fp_log);
        // if(edge1.cns_supp[i].spos > edge1.cns_supp[i].epos + 1)
        // {
        //     for(uint32_t ii = 0; ii < edge1.edge_supp.size(); ii++)
        //     {
        //         fprintf(fp_log, "        [edge_supp] lr_id:%u lr_rev:%c cmp_head_id:%u cmp_tail_id:%u\n",
        //             edge1.edge_supp[ii].lr_id, "+-"[edge1.edge_supp[ii].lr_strand], edge1.edge_supp[ii].cmp_head_id, edge1.edge_supp[ii].cmp_tail_id);
        //         uint32_t rid2 = edge1.edge_supp[ii].lr_id;
        //         uint32_t cmp1 = edge1.edge_supp[ii].cmp_head_id;
        //         uint32_t cmp2 = edge1.edge_supp[ii].cmp_tail_id;
        //         fprintf(fp_log, "\tsupp_detail head\t%u-%u:%u:%c:%u-%u\ttail\t%u-%u:%u:%c:%u-%u\n", 
        //             (*_asm_compact_lr_list)[rid2][cmp1]->q_start, (*_asm_compact_lr_list)[rid2][cmp1]->q_end, (*_asm_compact_lr_list)[rid2][cmp1]->t_id, "+-"[(*_asm_compact_lr_list)[rid2][cmp1]->is_rev], (*_asm_compact_lr_list)[rid2][cmp1]->t_start, (*_asm_compact_lr_list)[rid2][cmp1]->t_end, 
        //             (*_asm_compact_lr_list)[rid2][cmp2]->q_start, (*_asm_compact_lr_list)[rid2][cmp2]->q_end, (*_asm_compact_lr_list)[rid2][cmp2]->t_id, "+-"[(*_asm_compact_lr_list)[rid2][cmp2]->is_rev], (*_asm_compact_lr_list)[rid2][cmp2]->t_start, (*_asm_compact_lr_list)[rid2][cmp2]->t_end);
        //         fflush(fp_log);
        //     }
        //     // fflush(fp_log);
        // }
        string lr_subseq;
        if(rstrand == 0)
            lr_subseq = rseq.substr(edge1.cns_supp[i].spos, edge1.cns_supp[i].epos - edge1.cns_supp[i].spos + 1);
        else
            lr_subseq = rseq_rev.substr(edge1.cns_supp[i].spos, edge1.cns_supp[i].epos - edge1.cns_supp[i].spos + 1);
        // 
        fprintf(fp_log, ">%u %c %u %u %u\n", rid, (rstrand ? '-' : '+'), edge1.cns_supp[i].spos, edge1.cns_supp[i].epos, edge1.cns_supp[i].epos - edge1.cns_supp[i].spos + 1);
        fprintf(fp_log, "%s\n", lr_subseq.c_str());
        // fflush(fp_log);
        if(lr_subseq.size() > 0)
        {
            auto alignment = alignment_engine->align_sequence_with_graph(lr_subseq, graph_cns);
            graph_cns->add_alignment(alignment, lr_subseq);
            cnt_non_empty++;
        }
    } // END for(i = 0; i < edge1.cns_supp.size(); i++)
    if(cnt_non_empty == 0)
    {
        edge1.cns_seq = "";
        edge2.cns_seq = "";
        fprintf(fp_log, ">CONSENSUS\n");
        fprintf(fp_log, "%s\n", edge1.cns_seq.c_str());
        // fflush(fp_log);
    }
    else
    {
        edge1.cns_seq = graph_cns->generate_consensus();
        edge2.cns_seq = reverseComplement(edge1.cns_seq);
        fprintf(fp_log, ">CONSENSUS\n");
        fprintf(fp_log, "%s\n", edge1.cns_seq.c_str());
        // fflush(fp_log);
    }
}

void* asm_cal_cns_seq_ST(void *args_ptr)
{
    tuple<int, FILE*> args = *(tuple<int, FILE*> *)args_ptr;
    int th_id = std::get<0>(args);
    FILE *fp_log = std::get<1>(args);
    // 
    pair<uint32_t, map<uint32_t, BBG_Edge_t>::iterator> p1;
    while(asm_get_next_edge(p1))
    {
        uint32_t node_id = p1.first / 2;
        uint32_t node_rev = p1.first % 2;
        fprintf(fp_log, "calc_cns th_id:%u %u:%c -> %u:%c\n", th_id, node_id, "+-"[node_rev], p1.second->first >> 1,  "+-"[p1.second->first & 1]);
        fflush(fp_log);
        asm_calc_single_cns_seq(node_id, node_rev, p1.second->first >> 1, p1.second->first & 1, fp_log);
    }
    return NULL;
}

void asm_cal_cns_seq_MT(vector<BBG_Node_t> *bbgraph, Longread_List_t *lr_list, string log_path)
{
    _asm_bbgraph = bbgraph;
    _asm_lr_list = lr_list;
    // 
    if(asm_init_vertex_index(12) == false) return;
    // 
    // init arguments
    FILE *fp_log = file_open_write(log_path.c_str());
    vector<tuple<int, FILE*>> args_MT(gopt.num_threads);
    for(uint32_t i = 0; i < gopt.num_threads; i++)
    {
        std::get<0>(args_MT[i]) = i;
        std::get<1>(args_MT[i]) = fp_log;
    }
    // spawn threads and join
    // int TH_ID[gopt.num_threads]; for(uint32_t i = 0; i < gopt.num_threads; i++) TH_ID[i] = i;
    pthread_t *cns_threads = (pthread_t*)malloc(gopt.num_threads * sizeof(pthread_t));
    for(uint32_t i = 0; i < gopt.num_threads; i++)
        pthread_create(cns_threads + i, NULL, asm_cal_cns_seq_ST, &args_MT[i]);
    // 
    for (uint32_t i = 0; i < gopt.num_threads; i++)
        pthread_join(cns_threads[i], NULL);
    free(cns_threads);
    fclose(fp_log);
}

void asm_find_simple_path_from_source(vector<BBG_Node_t> &graph, uint32_t src_node, uint32_t src_strand, map<uint32_t, BBG_Edge_t>::iterator it, deque<id_strand2_t> &simple_path)
{
    simple_path.clear();
    simple_path.push_back({src_strand, src_node});
    uint32_t curr_node   = it->first >> 1;
    uint32_t curr_strand = it->first  & 1;
    while(true)
    {
        simple_path.push_back({curr_strand, curr_node});
        if(graph[curr_node].edges[curr_strand].size() == 0) break; // not extendable
        if(graph[curr_node].edges[curr_strand].size() > 1 || graph[curr_node].edges[1 - curr_strand].size() > 1) break; // end of simple path
        it = graph[curr_node].edges[curr_strand].begin();
        curr_node   = it->first >> 1;
        curr_strand = it->first  & 1;
    }
}

void asm_assemble_single_path(deque<id_strand2_t> &path, vector<BBG_Node_t> &graph, int &nb_ctg, FILE *fp_asm, FILE *fp_ann, FILE *fp_log)
{
    if(path.size() == 1)
    {
        uint32_t contig1 = path.front().id;
        uint32_t strand1 = path.front().strand;
        string contig1_str = get_uncompressed_dna(_asm_contig_list->contigs[graph[contig1].contig_id].comp_seq, _asm_contig_list->contigs[graph[contig1].contig_id].len, _asm_contig_list->contigs[graph[contig1].contig_id].comp_len);
        fprintf(fp_log, ">%d from:%u:%c to:%u:%c\n%s\n\n", nb_ctg, contig1, "+-"[strand1], contig1, "+-"[strand1], contig1_str.c_str());
        fprintf(fp_asm, ">%d from:%u:%c to:%u:%c\n%s\n", nb_ctg, contig1, "+-"[strand1], contig1, "+-"[strand1], contig1_str.c_str());
        nb_ctg++;        
        return;
    }
    // if(path.size() == 1)
    // {
    //     uint32_t singleton_id = path[0].id;
    //     uint32_t singleton_strand = path[0].strand;
    //     string singleton_str = get_uncompressed_dna(_asm_contig_list->contigs[singleton_id].comp_seq, _asm_contig_list->contigs[singleton_id].len, _asm_contig_list->contigs[singleton_id].comp_len);
    //     fprintf(stdout, ">%d from:%u:%c to:%u:%c singletone\n%s\n\n", nb_ctg, singleton_id, "+-"[singleton_strand], singleton_id, "+-"[singleton_strand], singleton_str.c_str());
    //     fflush(stdout);
    //     fprintf(fp_asm, ">%d from:%u:%c to:%u:%c singletone\n%s\n", nb_ctg, singleton_id, "+-"[singleton_strand], singleton_id, "+-"[singleton_strand], singleton_str.c_str());
    //     nb_ctg++;
    //     return;
    // }
    // 
    // for(uint32_t i = 0; i < path.size() - 1; i++)
    // {
    //     uint32_t node1 = path[i].id;
    //     uint32_t rev1 = path[i].strand;
    //     uint32_t node2 = path[i+1].id;
    //     uint32_t rev2 = path[i+1].strand;
    //     uint32_t to_node1 = (node2 << 1) | rev2;
    //     uint32_t to_node2 = (node1 << 1) | (1 - rev1);
    //     graph[node1].edges[rev1][to_node1].flag = _asm_flag_visited;
    //     graph[node2].edges[1-rev2][to_node2].flag = _asm_flag_visited;
    // }
    // 
    string assembled = "";
    uint32_t source_contig = path[0].id;
    uint32_t source_strand = path[0].strand;
    uint32_t contig1_start = (source_strand == 0 ? 0 : _asm_contig_list->contigs[graph[source_contig].contig_id].len - 1);
    uint32_t target_contig = path[path.size() - 1].id;
    uint32_t target_strand = path[path.size() - 1].strand;
    uint32_t i;
    for(i = 0; i < path.size() - 1; i++)
    {
        uint32_t contig1 = path[i].id;
        uint32_t strand1 = path[i].strand;
        string contig1_str = get_uncompressed_dna(_asm_contig_list->contigs[graph[contig1].contig_id].comp_seq, _asm_contig_list->contigs[graph[contig1].contig_id].len, _asm_contig_list->contigs[graph[contig1].contig_id].comp_len);
        uint32_t contig2 = path[i+1].id;
        uint32_t strand2 = path[i+1].strand;
        // update flags
        uint32_t to_node1 = (contig2 << 1) | strand2;
        // uint32_t to_node2 = (contig1 << 1) | (1 - strand1);
        BBG_Edge_t &edge1 = graph[contig1].edges[strand1][to_node1];
        // BBG_Edge_t &edge2 = graph[contig2].edges[1-strand2][to_node2];
        // edge1.flag = _asm_flag_visited;
        // edge2.flag = _asm_flag_visited;
        string prefix = "";
        if(edge1.cns_supp.size() == 0) // break the assembly
        {
            fprintf(fp_log, "[breaking] contig1_len:%zu    contig1_start:%u    prev_end:%u     next_beg:%u\n", contig1_str.size(), contig1_start, edge1.head_end, edge1.tail_beg);
            if(strand1 == 0)
            {
                prefix = contig1_str.substr(contig1_start);
                fprintf(fp_ann, "%d\t%zu\t%zu\tctg\t+\t%u\t%zu\t%u\t%zu\n", nb_ctg, assembled.size(), assembled.size() + prefix.size(), contig1, contig1_str.size(), contig1_start, contig1_str.size());
            }
            else
            {
                prefix = contig1_str.substr(0, contig1_start + 1);
                fprintf(fp_ann, "%d\t%zu\t%zu\tctg\t-\t%u\t%zu\t%u\t%u\n", nb_ctg, assembled.size(), assembled.size() + prefix.size(), contig1, contig1_str.size(), 0, contig1_start + 1);
                prefix = reverseComplement(prefix);
            }
            assembled += prefix;
            fprintf(fp_log, ">%d from:%u:%c to:%u:%c\n%s\n\n", nb_ctg, source_contig, "+-"[source_strand], contig1, "+-"[strand1], assembled.c_str());
            // fflush(fp_log);
            fprintf(fp_asm, ">%d from:%u:%c to:%u:%c\n%s\n", nb_ctg, source_contig, "+-"[source_strand], contig1, "+-"[strand1], assembled.c_str());
            nb_ctg++;
            assembled = "";
            source_contig = contig2;
            source_strand = strand2;
            contig1_start = (source_strand == 0 ? 0 : _asm_contig_list->contigs[graph[source_contig].contig_id].len - 1);
            fprintf(stderr, "[WARNING] breaking assembly for path %u:%c --> %u:%c between anchors %u:%c --> %u:%c\n", source_contig, "+-"[source_strand], target_contig, "+-"[target_strand], contig1, "+-"[strand1], contig2, "+-"[strand2]);
        }
        else
        {
            fprintf(fp_log, "[stitching] contig1_len:%zu    contig1_start:%u    prev_end:%u     next_beg:%u\n", contig1_str.size(), contig1_start, edge1.head_end, edge1.tail_beg);
            if (strand1 == 0)
            {
                prefix = contig1_str.substr(contig1_start, edge1.head_end - contig1_start + 1);
                fprintf(fp_ann, "%d\t%zu\t%zu\tctg\t+\t%u\t%zu\t%u\t%zu\n", nb_ctg, assembled.size(), assembled.size() + prefix.size(), contig1, contig1_str.size(), contig1_start, contig1_start + prefix.size());
            }
            else
            {
                prefix = contig1_str.substr(edge1.head_end, contig1_start - edge1.head_end + 1);
                fprintf(fp_ann, "%d\t%zu\t%zu\tctg\t-\t%u\t%zu\t%u\t%zu\n", nb_ctg, assembled.size(), assembled.size() + prefix.size(), contig1, contig1_str.size(), edge1.head_end, edge1.head_end + prefix.size());
                prefix = reverseComplement(prefix);
            }
            assembled += prefix;
            fprintf(fp_ann, "%d\t%zu\t%zu\tcns\t%zu\t%zu\n", nb_ctg, assembled.size(), assembled.size() + edge1.cns_seq.size(), edge1.cns_seq.size(), edge1.cns_supp.size());
            assembled += edge1.cns_seq;
            if(strand2 == 0)
            {
                contig1_start = edge1.tail_beg;
            }
            else
            {
                contig1_start = edge1.tail_beg;
            }
        }
    }
    // concatenate the last contig suffix
    uint32_t contig2 = path[i].id;
    uint32_t strand2 = path[i].strand;
    string contig2_str = get_uncompressed_dna(_asm_contig_list->contigs[graph[contig2].contig_id].comp_seq, _asm_contig_list->contigs[graph[contig2].contig_id].len, _asm_contig_list->contigs[graph[contig2].contig_id].comp_len);
    string suffix = "";
    if (strand2 == 0)
    {
        suffix = contig2_str.substr(contig1_start);
        fprintf(fp_ann, "%d\t%zu\t%zu\tctg\t+\t%u\t%zu\t%u\t%zu\n", nb_ctg, assembled.size(), assembled.size() + suffix.size(), contig2, contig2_str.size(), contig1_start, contig2_str.size());
    }
    else
    {
        suffix = contig2_str.substr(0, contig1_start + 1);
        fprintf(fp_ann, "%d\t%zu\t%zu\tctg\t-\t%u\t%zu\t%u\t%u\n", nb_ctg, assembled.size(), assembled.size() + suffix.size(), contig2, contig2_str.size(), 0, contig1_start + 1);
        suffix = reverseComplement(suffix);
    }
    assembled += suffix;
    fprintf(fp_log, ">%d from:%u:%c to:%u:%c\n%s\n\n", nb_ctg, source_contig, "+-"[source_strand], contig2, "+-"[strand2], assembled.c_str());
    // fflush(fp_log);
    fprintf(fp_asm, ">%d from:%u:%c to:%u:%c\n%s\n", nb_ctg, source_contig, "+-"[source_strand], contig2, "+-"[strand2], assembled.c_str());
    nb_ctg++;
}

void asm_extract_all_simple_paths(vector<BBG_Node_t> &graph, vector<deque<id_strand2_t>> &path_list)
{
    for(uint32_t i = 0; i < graph.size(); i++)
    {
        if(graph[i].edges[0].size() == 1 && graph[i].edges[1].size() == 1)
        {
            // cannot be the start of a path
            continue;
        }
        if(graph[i].edges[0].size() > 1 && graph[i].edges[1].size() > 1)
        {
            deque<id_strand2_t> path;
            path.push_back({0, i});
            path_list.push_back(path);
            // continue;
        }
        for(uint32_t rev = 0; rev < 2; rev++)
        {
            for(map<uint32_t, BBG_Edge_t>::iterator it = graph[i].edges[rev].begin(); it != graph[i].edges[rev].end(); it++)
            {
                if(it->second.flag != _asm_flag_visited) // an unvisited simple path from here
                {
                    deque<id_strand2_t> path;
                    asm_find_simple_path_from_source(graph, i, rev, it, path);
                    // 
                    uint32_t node1, node2, rev1, rev2, to_node1, to_node2;
                    for(uint32_t j = 0; j < path.size() - 1; j++)
                    {
                        node1 = path[j].id;
                        rev1 = path[j].strand;
                        node2 = path[j+1].id;
                        rev2 = path[j+1].strand;
                        to_node1 = (node2 << 1) | rev2;
                        to_node2 = (node1 << 1) | (1 - rev1);
                        graph[node1].edges[rev1][to_node1].flag = _asm_flag_visited;
                        graph[node2].edges[1-rev2][to_node2].flag = _asm_flag_visited;
                    }
                    // trim the front and back if necessary
                    node1 = path.front().id;
                    rev1 = path.front().strand;
                    if(graph[node1].edges[rev1].size() > 1) path.pop_front();
                    node2 = path.back().id;
                    rev2 = path.back().strand;
                    if(graph[node2].edges[1-rev2].size() > 1) path.pop_back();
                    // 
                    if(path.size() > 0)
                    {
                        path_list.push_back(path);
                    }
                }
            }
        }
    }
}

uint32_t asm_get_shared_supp(vector<Edge_Supp_t> &supp1, vector<Edge_Supp_t> &supp2)
{
    set<uint32_t> lrs1;
    for(uint32_t i = 0; i < supp1.size(); i++)
        lrs1.insert(supp1[i].lr_id);
    set<uint32_t> lrs2;
    for(uint32_t i = 0; i < supp2.size(); i++)
        lrs2.insert(supp2[i].lr_id);
    set<uint32_t> shared_lrs;
    std::set_intersection(lrs1.begin(), lrs1.end(), lrs2.begin(), lrs2.end(), std::inserter(shared_lrs, shared_lrs.begin()));
    return shared_lrs.size();
}

void asm_connect_paths(vector<deque<id_strand2_t>> &path_list, uint32_t middle_path,
    map<uint32_t, BBG_Edge_t>::iterator &in_edge, map<uint32_t, BBG_Edge_t>::iterator &out_edge,
    map<pair<uint32_t, uint32_t>, pair<uint32_t, uint32_t>> &tails, vector<uint8_t> &deleted_paths,
    bool delete_middle)
{
    map<pair<uint32_t, uint32_t>, pair<uint32_t, uint32_t>>::iterator it_in, it_out;
    deque<id_strand2_t> path_merged;
    it_in = tails.find(make_pair(uint32_t(in_edge->first >> 1), uint32_t(in_edge->first & 1)));
    if(it_in == tails.end())
    {
        fprintf(stdout, "[ERROR] here1\tin_edge not found\n");
        fflush(stdout);
        exit(EXIT_FAILURE);
    }
    it_out = tails.find(make_pair(uint32_t(out_edge->first >> 1), uint32_t(out_edge->first & 1)));
    if(it_out == tails.end())
    {
        fprintf(stdout, "[ERROR] here2\tout_edge not found\n");
        fflush(stdout);
        exit(EXIT_FAILURE);
    }
    uint32_t pid1 = it_in->second.first;
    uint32_t pid2 = it_out->second.first;
    if(pid1 == pid2)
    {
        fprintf(stdout, "[warning] weird case!\n");
        fflush(stdout);
        // exit(EXIT_FAILURE);
        deleted_paths[middle_path] = 1;
        return;
    }
    // append ingoing path
    if(it_in->second.second == 0) // front; reverse of the path
    {
        for(deque<id_strand2_t>::reverse_iterator z = path_list[pid1].rbegin(); z != path_list[pid1].rend(); z++)
            path_merged.push_back({uint32_t(1 - z->strand), z->id});
    }
    else // back; same path
    {
        for(deque<id_strand2_t>::iterator z = path_list[pid1].begin(); z != path_list[pid1].end(); z++)
            path_merged.push_back({z->strand, z->id});
    }
    // append middle path
    for(deque<id_strand2_t>::iterator z = path_list[middle_path].begin(); z != path_list[middle_path].end(); z++)
        path_merged.push_back({z->strand, z->id});
    // append outgoing path
    if(it_out->second.second == 0) // front; same path
    {
        for(deque<id_strand2_t>::iterator z = path_list[pid2].begin(); z != path_list[pid2].end(); z++)
            path_merged.push_back({z->strand, z->id});
    }
    else // back; reverse of the path
    {
        for(deque<id_strand2_t>::reverse_iterator z = path_list[pid2].rbegin(); z != path_list[pid2].rend(); z++)
            path_merged.push_back({uint32_t(1 - z->strand), z->id});
    }
    // erase old tails
    tails.erase(make_pair((uint32_t)path_list[pid1].front().id, (uint32_t)path_list[pid1].front().strand));
    tails.erase(make_pair((uint32_t)path_list[pid1].back().id,  (uint32_t)(1 - path_list[pid1].back().strand)));
    tails.erase(make_pair((uint32_t)path_list[pid2].front().id, (uint32_t)path_list[pid2].front().strand));
    tails.erase(make_pair((uint32_t)path_list[pid2].back().id,  (uint32_t)(1 - path_list[pid2].back().strand)));
    if(delete_middle == true)
    {
        tails.erase(make_pair((uint32_t)path_list[middle_path].front().id, (uint32_t)path_list[middle_path].front().strand));
        tails.erase(make_pair((uint32_t)path_list[middle_path].back().id,  (uint32_t)(1 - path_list[middle_path].back().strand)));
    }
    // add new tails
    tails[make_pair((uint32_t)path_merged.front().id, (uint32_t)path_merged.front().strand)] = make_pair(pid1, 0); // front
    tails[make_pair((uint32_t)path_merged.back().id,  (uint32_t)(1 - path_merged.back().strand))]  = make_pair(pid1, 1); // back
    // update paths
    path_list[pid1] = path_merged;
    deleted_paths[pid2] = 1;
    if(delete_middle == true)
    {
        deleted_paths[middle_path] = 1;
    }
    fprintf(stdout, "\tmerged:");
    for(uint32_t ii = 0; ii < path_merged.size(); ii++)
        fprintf(stdout, "\t%u:%c", path_merged[ii].id, "+-"[path_merged[ii].strand]);
    fprintf(stdout, "\n");
    fflush(stdout);
}

void asm_resolve_4way_nodes(vector<BBG_Node_t> &graph, vector<deque<id_strand2_t>> &path_list, vector<uint8_t> &deleted_paths)
{
    map<pair<uint32_t, uint32_t>, pair<uint32_t, uint32_t>> tails;
    for(uint32_t i = 0; i < path_list.size(); i++)
    {
        tails[make_pair((uint32_t)path_list[i].front().id, (uint32_t)path_list[i].front().strand)] = make_pair(i, 0); // front
        tails[make_pair((uint32_t)path_list[i].back().id,  (uint32_t)(1 - path_list[i].back().strand))]  = make_pair(i, 1); // back
    }
    // 
    for(uint32_t i = 0; i < path_list.size(); i++)
    {
        uint32_t node1, strand1, node2, strand2;
        node1 = path_list[i].front().id;
        strand1 = path_list[i].front().strand;
        node2 = path_list[i].back().id;
        strand2 = path_list[i].back().strand;
        if(graph[node2].edges[strand2].size() == 2 && graph[node1].edges[1-strand1].size() == 2)
        {
            fprintf(stdout, "from:%u:%c\tto:%u:%c\n", node1, "+-"[strand1], node2, "+-"[strand2]);
            map<uint32_t, BBG_Edge_t>::iterator in1, in2, out1, out2;
            in1 = in2 = graph[node1].edges[1-strand1].begin();
            in2++;
            out1 = out2 = graph[node2].edges[strand2].begin();
            out2++;
            fprintf(stdout, "\tin1 %u:%c supp:%zu\n", in1->first >> 1, "+-"[in1->first & 1], in1->second.edge_supp.size());
            fprintf(stdout, "\tin2 %u:%c supp:%zu\n", in2->first >> 1, "+-"[in2->first & 1], in2->second.edge_supp.size());
            fprintf(stdout, "\tout1 %u:%c supp:%zu\n", out1->first >> 1, "+-"[out1->first & 1], out1->second.edge_supp.size());
            fprintf(stdout, "\tout2 %u:%c supp:%zu\n", out2->first >> 1, "+-"[out2->first & 1], out2->second.edge_supp.size());
            uint32_t sup11, sup12, sup21, sup22;
            sup11 = asm_get_shared_supp(in1->second.edge_supp, out1->second.edge_supp); // in1-out1
            sup12 = asm_get_shared_supp(in1->second.edge_supp, out2->second.edge_supp); // in1-out2
            sup21 = asm_get_shared_supp(in2->second.edge_supp, out1->second.edge_supp); // in2-out1
            sup22 = asm_get_shared_supp(in2->second.edge_supp, out2->second.edge_supp); // in2-out2
            fprintf(stdout, "\tin1-out1 %u\n", sup11);
            fprintf(stdout, "\tin1-out2 %u\n", sup12);
            fprintf(stdout, "\tin2-out1 %u\n", sup21);
            fprintf(stdout, "\tin2-out2 %u\n", sup22);
            if((sup11 > 2 * sup12 && !(sup21 > 2 * sup22)) || (sup22 > 2 * sup21 && !(sup12 > 2 *sup11)))
            {
                // connect in1-out1
                asm_connect_paths(path_list, i, in1, out1, tails, deleted_paths, false);
                // connect in2-out2
                asm_connect_paths(path_list, i, in2, out2, tails, deleted_paths, true);
            }
            if((sup12 > 2 * sup11 && !(sup22 > 2 * sup21)) || (sup21 > 2 * sup22 && !(sup11 > 2 *sup12)))
            {
                // connect in1-out2
                asm_connect_paths(path_list, i, in1, out2, tails, deleted_paths, false);
                // connect in2-out1
                asm_connect_paths(path_list, i, in2, out1, tails, deleted_paths, true);
            }
        }
    }
}

void asm_identify_unused_longreads(vector<BBG_Node_t> &graph, vector<deque<id_strand2_t>> &path_list)
{
    uint32_t i, j;
    vector<uint8_t> unused(_asm_lr_list->reads_size, 1);
    // mark reads used for uniq assembly
    for(i = 0; i < path_list.size(); i++)
    {
        for(j = 0; j < path_list[i].size(); j++)
        {
            uint32_t contig1 = path_list[i][j].id;
            for(uint32_t rev1 = 0; rev1 < 2; rev1++)
            {
                for(map<uint32_t, BBG_Edge_t>::iterator it = graph[contig1].edges[rev1].begin(); it != graph[contig1].edges[rev1].end(); it++)
                {
                    for(uint32_t x = 0; x < it->second.edge_supp.size(); x++)
                    {
                        unused[it->second.edge_supp[x].lr_id] = 0;
                    }
                }
            }
        }
    }
    // include reads covering two extremeties
    for(i = 0; i < path_list.size(); i++)
    {
        uint32_t contig1 = path_list[i].front().id;
        for(uint32_t rev1 = 0; rev1 < 2; rev1++)
        {
            for(map<uint32_t, BBG_Edge_t>::iterator it = graph[contig1].edges[rev1].begin(); it != graph[contig1].edges[rev1].end(); it++)
            {
                for(uint32_t x = 0; x < it->second.edge_supp.size(); x++)
                {
                    unused[it->second.edge_supp[x].lr_id] = 2;
                }
            }
        }
        // 
        contig1 = path_list[i].back().id;
        for(uint32_t rev1 = 0; rev1 < 2; rev1++)
        {
            for(map<uint32_t, BBG_Edge_t>::iterator it = graph[contig1].edges[rev1].begin(); it != graph[contig1].edges[rev1].end(); it++)
            {
                for(uint32_t x = 0; x < it->second.edge_supp.size(); x++)
                {
                    unused[it->second.edge_supp[x].lr_id] = 2;
                }
            }
        }
    }
    // 

    FILE *fp = file_open_write(gopt.out_dir + "/lr.unused.fasta");
    if(fp == NULL)
    {
        fprintf(stderr, "[ERROR] could not open file: %s\n", (gopt.out_dir + "/lr.unused.fasta").c_str());
        exit(EXIT_FAILURE);
    }
    for(i = 0; i < _asm_lr_list->reads_size; i++)
    {
        if(unused[i])
        {
            string tmp_seq = get_uncompressed_dna(_asm_lr_list->reads[i].comp_seq, _asm_lr_list->reads[i].len, _asm_lr_list->reads[i].comp_len);
            fprintf(fp, ">u%u %s\n%s\n", i, (unused[i] == 2 ? "tail" : ""), tmp_seq.c_str());
        }
    }
    fclose(fp);
    // 
    // fp = fopen((gopt.outDir + "/unused.txt").c_str(), "wt");
    // if(fp == NULL)
    // {
    //     fprintf(stderr, "[ERROR] could not open file: %s\n", (gopt.outDir + "/unused.txt").c_str());
    //     exit(EXIT_FAILURE);
    // }
    // // for(i = 0; i < lr_list.reads_size; i++)
    // for(i = 0; i < compact_lr_list.size(); i++)
    // {
    //     // fprintf(fp, "%u\t%u\n", i, (unused[i] ? 1 : 0));
    //     fprintf(fp, "%u\t%u\n", i, unused[i]);
    // }
    // fclose(fp);
}

void asm_get_assembly(vector<BBG_Node_t> &graph)
{
    // string asm_path = gopt.out_dir + "/asm.backbone.fa";
    // string ann_path = gopt.out_dir + "/asm.backbone.ann";
    // string log_path = gopt.out_dir + "/log_asmbackbone.txt";
    string asm_path = gopt.out_dir + "/asm.final.fa";
    string ann_path = gopt.out_dir + "/asm.final.ann";
    string log_path = gopt.out_dir + "/log_asmfinal.txt";
    FILE *fp_asm = file_open_write(asm_path);
    FILE *fp_ann = file_open_write(ann_path);
    FILE *fp_log = file_open_write(log_path);

    _asm_flag_visited = 21;
    vector<deque<id_strand2_t>> path_list;
    asm_extract_all_simple_paths(graph, path_list);
    for(uint32_t i = 0; i < path_list.size(); i++)
    {
        fprintf(fp_log, "simple_path %u size:%zu\tfrom:%u:%c\tto:%u:%c\n", i, path_list[i].size(), path_list[i].front().id, "+-"[path_list[i].front().strand], path_list[i].back().id, "+-"[path_list[i].back().strand]);
        // for(uint32_t j = 0; j < path_list[i].size(); j++)
        //     fprintf(stdout, "%u:%c\t", path_list[i][j].id, "+-"[path_list[i][j].strand]);
        // fprintf(stdout, "\n");
    }
    // asm_identify_unused_longreads(graph, path_list);
    // 
    int nb_ctg = 0;
    for(uint32_t i = 0; i < path_list.size(); i++)
    {
        asm_assemble_single_path(path_list[i], graph, nb_ctg, fp_asm, fp_ann, fp_log);
    }
    // 
    fclose(fp_log);
    fclose(fp_ann);
    fclose(fp_asm);
    // //////////////////////////////////////////////////////////////// 
    // asm_path = gopt.out_dir + "/asm.final.fa";
    // ann_path = gopt.out_dir + "/asm.final.ann";
    // log_path = gopt.out_dir + "/log_asmfinal.txt";
    // fp_asm = file_open_write(asm_path);
    // fp_ann = file_open_write(ann_path);
    // fp_log = file_open_write(log_path);

    // vector<uint8_t> deleted_paths(path_list.size(), 0);
    // asm_resolve_4way_nodes(graph, path_list, deleted_paths);
    // // 
    // for(uint32_t i = 0; i < path_list.size(); i++)
    // {
    //     if(deleted_paths[i] == 0)
    //     {
    //         fprintf(fp_log, "final_path size:%zu\tfrom:%u:%c\tto:%u:%c\n", path_list[i].size(), path_list[i].front().id, "+-"[path_list[i].front().strand], path_list[i].back().id, "+-"[path_list[i].back().strand]);
    //         // for(uint32_t j = 0; j < path_list[i].size(); j++)
    //         //     fprintf(stdout, "%u:%c\t", path_list[i][j].id, "+-"[path_list[i][j].strand]);
    //         // fprintf(stdout, "\n");
    //     }
    // }
    // // 
    // nb_ctg = 0;
    // for(uint32_t i = 0; i < path_list.size(); i++)
    // {
    //     if(deleted_paths[i] == 0)
    //     {
    //         asm_assemble_single_path(path_list[i], graph, nb_ctg, fp_asm, fp_ann, fp_log);
    //     }
    // }
    // // 
    // fclose(fp_log);
    // fclose(fp_ann);
    // fclose(fp_asm);
}
