// /*****************************************************
//  * Author: Ehsan Haghshenas (ehaghshe AT sfu DOT ca) *
//  *****************************************************/

#include "Graph_repeat.hpp"
#include "Common.hpp"

void lcs_alignment(vector<vector<Align_Seq_t*> > &compact_lr_list, uint32_t lr1, uint32_t lr2, deque<int32_t> &aln1, deque<int32_t> &aln2, int &score)
{
    int LCS_MATCH = 3;
    int LCS_INDEL = -1;
    int i, j;
    vector<Align_Seq_t*> &lr1_comp = compact_lr_list[lr1];
    vector<Align_Seq_t*> &lr2_comp = compact_lr_list[lr2];
    int m = lr1_comp.size();
    int n = lr2_comp.size();
    // fprintf(stderr, "[lcs_alignment] compact_lr_list: %zu\t lr1: %u\t lr1_size: %zu\t lr2: %u\t lr2_size: %zu\n", compact_lr_list.size(), lr1, lr1_comp.size(), lr2, lr2_comp.size());
    int lcs[1000][1000];
    char bt[1000][1000]; // D diagonal, L left, U up
    // init the first column, do not penalize beginning gaps
    for(i = 0; i <= m; i++)
    {
        lcs[i][0] = 0;
        bt[i][0] = 'U'; // up
    }
    // init the first row, do not penalize beginning gaps
    for(j = 0; j <= n; j++)
    {
        lcs[0][j] = 0;
        bt[0][j] = 'L'; // left
    }
    // DP
    for(i = 1; i <= m; i++)
    {
        for(j = 1; j <= n; j++)
        {
            if(lr1_comp[i-1]->t_id == lr2_comp[j-1]->t_id && lr1_comp[i-1]->is_rev == lr2_comp[j-1]->is_rev)
            {
                lcs[i][j] = lcs[i-1][j-1] + LCS_MATCH;
                bt[i][j] = 'D'; // diagonal
            }
            else
            {
                if(lcs[i-1][j] > lcs[i][j-1])
                {
                    lcs[i][j] = lcs[i-1][j] + LCS_INDEL;
                    bt[i][j] = 'U'; // up
                }
                else
                {
                    lcs[i][j] = lcs[i][j-1] + LCS_INDEL;
                    bt[i][j] = 'L'; // lef
                }
            }
        }
    }
    // do not penalize end gaps
    for(i = 0; i < m; i++)
    {
        if(lcs[i][n] > lcs[i+1][n])
        {
            lcs[i+1][n] = lcs[i][n];
            bt[i+1][n] = 'U';
        }
    }
    // do not penalize end gaps
    for(j = 0; j < n; j++)
    {
        if(lcs[m][j] > lcs[m][j+1])
        {
            lcs[m][j+1] = lcs[m][j];
            bt[m][j+1] = 'L';
        }
    }
    // BT
    // deque<pair<int32_t, int8_t>> aln1, aln2;
    score = lcs[m][n];
    aln1.clear();
    aln2.clear();
    i = m;
    j = n;
    while(i > 0 || j > 0)
    {
        if(bt[i][j] == 'L') // left
        {
            aln1.push_front(-1);
            aln2.push_front(j-1);
            j--;
        }
        else if(bt[i][j] == 'U') // top
        {
            aln1.push_front(i-1);
            aln2.push_front(-1);
            i--;
        }
        else // diagonal
        {
            aln1.push_front(i-1);
            aln2.push_front(j-1);
            i--;
            j--;
        }
    }
    
    // // print dp matrix
    // fprintf(stdout, "$$\t$$\t");
    // for(j = 0; j < n; j++)
    //     fprintf(stdout, "%u\t", lr2_comp[j]->t_id);
    // fprintf(stdout, "\n");
    // for(i = 0; i <= m; i++)
    // {
    //     if(i == 0)
    //         fprintf(stdout, "$$\t");
    //     else
    //         fprintf(stdout, "%u\t", lr1_comp[i-1]->t_id);
    //     for(j = 0; j <= n; j++)
    //     {
    //         fprintf(stdout, "%2d:%c\t", lcs[i][j], bt[i][j]);
    //     }
    //     fprintf(stdout, "\n");
    // }
}

void lcs_alignment_reverse(vector<vector<Align_Seq_t*> > &compact_lr_list, uint32_t lr1, uint32_t lr2, deque<int32_t> &aln1, deque<int32_t> &aln2, int &score)
{
    int LCS_MATCH = 3;
    int LCS_INDEL = -1;
    int i, j;
    vector<Align_Seq_t*> &lr1_comp = compact_lr_list[lr1];
    vector<Align_Seq_t*> &lr2_comp = compact_lr_list[lr2];
    int m = lr1_comp.size();
    int n = lr2_comp.size();
    int lcs[1000][1000];
    char bt[1000][1000]; // D diagonal, L left, U up
    // init the first column, do not penalize beginning gaps
    for(i = 0; i <= m; i++)
    {
        lcs[i][0] = 0;
        bt[i][0] = 'U'; // up
    }
    // init the first row, do not penalize beginning gaps
    for(j = 0; j <= n; j++)
    {
        lcs[0][j] = 0;
        bt[0][j] = 'L'; // left
    }
    // DP
    for(i = 1; i <= m; i++)
    {
        for(j = 1; j <= n; j++)
        {
            if(lr1_comp[m-i]->t_id == lr2_comp[j-1]->t_id && lr1_comp[m-i]->is_rev != lr2_comp[j-1]->is_rev)
            {
                lcs[i][j] = lcs[i-1][j-1] + LCS_MATCH;
                bt[i][j] = 'D'; // diagonal
            }
            else
            {
                if(lcs[i-1][j] > lcs[i][j-1])
                {
                    lcs[i][j] = lcs[i-1][j] + LCS_INDEL;
                    bt[i][j] = 'U'; // up
                }
                else
                {
                    lcs[i][j] = lcs[i][j-1] + LCS_INDEL;
                    bt[i][j] = 'L'; // lef
                }
            }
        }
    }
    // do not penalize end gaps
    for(i = 0; i < m; i++)
    {
        if(lcs[i][n] > lcs[i+1][n])
        {
            lcs[i+1][n] = lcs[i][n];
            bt[i+1][n] = 'U';
        }
    }
    // do not penalize end gaps
    for(j = 0; j < n; j++)
    {
        if(lcs[m][j] > lcs[m][j+1])
        {
            lcs[m][j+1] = lcs[m][j];
            bt[m][j+1] = 'L';
        }
    }
    // BT
    // deque<pair<int32_t, int8_t>> aln1, aln2;
    score = lcs[m][n];
    aln1.clear();
    aln2.clear();
    i = m;
    j = n;
    while(i > 0 || j > 0)
    {
        if(bt[i][j] == 'L') // left
        {
            aln1.push_front(-1);
            aln2.push_front(j-1);
            j--;
        }
        else if(bt[i][j] == 'U') // top
        {
            aln1.push_front(m-i);
            aln2.push_front(-1);
            i--;
        }
        else // diagonal
        {
            aln1.push_front(m-i);
            aln2.push_front(j-1);
            i--;
            j--;
        }
    }
    
    // // print dp matrix
    // fprintf(stdout, "$$\t$$\t");
    // for(j = 0; j < n; j++)
    //     fprintf(stdout, "%u\t", lr2_comp[j]->t_id);
    // fprintf(stdout, "\n");
    // for(i = 0; i <= m; i++)
    // {
    //     if(i == 0)
    //         fprintf(stdout, "$$\t");
    //     else
    //         fprintf(stdout, "%u\t", lr1_comp[m-i]->t_id);
    //     for(j = 0; j <= n; j++)
    //     {
    //         fprintf(stdout, "%2d:%c\t", lcs[i][j], bt[i][j]);
    //     }
    //     fprintf(stdout, "\n");
    // }
}

// int asm_is_alignment_spurious(vector<vector<Align_Seq_t*> > &compact_lr_list, uint32_t lr1, uint32_t lr2, deque<int32_t> &aln1, deque<int32_t> &aln2)
// {
//     int gap_num = 0;
//     int gap_bp = 0;
//     vector<Align_Seq_t*> &lr1_comp = compact_lr_list[lr1];
//     vector<Align_Seq_t*> &lr2_comp = compact_lr_list[lr2];
//     for(int32_t i = 0; i < (int32_t)aln1.size(); i++)
//     {
//         if(aln1[i] == -1)
//         {
//             gap_num++;
//             // gap_bp += lr2_comp[aln2[i]]->n_block;
//             gap_bp += (lr2_comp[aln2[i]]->q_end - lr2_comp[aln2[i]]->q_start);
//         }
//         if(aln2[i] == -1)
//         {
//             gap_num++;
//             // gap_bp += lr1_comp[aln1[i]]->n_block;
//             gap_bp += (lr1_comp[aln1[i]]->q_end - lr1_comp[aln1[i]]->q_start);
//         }
//     }
//     for(int32_t i = 0; i < (int32_t)aln1.size() && aln1[i] == -1; i++)
//     {
//         gap_num--;
//         // gap_bp -= lr2_comp[aln2[i]]->n_block;
//         gap_bp -= (lr2_comp[aln2[i]]->q_end - lr2_comp[aln2[i]]->q_start);
//     }
//     for(int32_t i = 0; i < (int32_t)aln2.size() && aln2[i] == -1; i++)
//     {
//         gap_num--;
//         // gap_bp -= lr1_comp[aln1[i]]->n_block;
//         gap_bp -= (lr1_comp[aln1[i]]->q_end - lr1_comp[aln1[i]]->q_start);
//     }
//     for(int32_t i = aln1.size() - 1; i >= 0 && aln1[i] == -1; i--)
//     {
//         gap_num--;
//         // gap_bp -= lr2_comp[aln2[i]]->n_block;
//         gap_bp -= (lr2_comp[aln2[i]]->q_end - lr2_comp[aln2[i]]->q_start);
//     }
//     for(int32_t i = aln2.size() - 1; i >= 0 && aln2[i] == -1; i--)
//     {
//         gap_num--;
//         // gap_bp -= lr1_comp[aln1[i]]->n_block;
//         gap_bp -= (lr1_comp[aln1[i]]->q_end - lr1_comp[aln1[i]]->q_start);
//     }
    
//     fprintf(stdout, "gap_num: %d\n", gap_num);
//     fprintf(stdout, "gap_bp:  %d\n", gap_bp);
//     // if(gap_num > 4 || gap_bp > 7000)
//     if(gap_bp > 8000)
//         return 1;
//     else
//         return 0;
// }

// #include "needleman_wunsch.h"
// int asm_is_overlap_spurious(uint32_t lr1, uint32_t lr2)
// {
//     bool no_start_gap_penalty = true, no_end_gap_penalty = true;
//     bool no_gaps_in_a = false, no_gaps_in_b = false;
//     bool no_mismatches = false, case_sensitive = true;
//     int match = 2;
//     int mismatch = -4;
//     int gap_open = -4;
//     int gap_extend = -1;

//     scoring_t scoring;
//     scoring_init(&scoring, match, mismatch, gap_open, gap_extend,
//         no_start_gap_penalty, no_end_gap_penalty, 
//         no_gaps_in_a, no_gaps_in_b,
//         no_mismatches, case_sensitive);

//     // Alignment results stored here
//     nw_aligner_t *nw = needleman_wunsch_new();
//     alignment_t *aln = alignment_create(256);

//     string aseq = get_uncompressed_dna(_asm_lr_list->reads[lr1].comp_seq, _asm_lr_list->reads[lr1].len, _asm_lr_list->reads[lr1].comp_len);
//     string bseq = get_uncompressed_dna(_asm_lr_list->reads[lr2].comp_seq, _asm_lr_list->reads[lr2].len, _asm_lr_list->reads[lr2].comp_len);

//     needleman_wunsch_align(aseq.c_str(), bseq.c_str(), &scoring, nw, aln);

//     int n_match = 0;
//     for(uint32_t i=0; i < aln->length; i++)
//     {
//         // cout<< "[debug] " << aln->result_a[i] << " " << aln->result_b[i] << endl;
//         n_match += (aln->result_a[i] == aln->result_b[i] ? 1 : 0);
//     }

//     fprintf(stdout, "%s\n", aln->result_a);
//     fprintf(stdout, "%s\n", aln->result_b);
//     fprintf(stdout, "score: %d\n", aln->score);
//     fprintf(stdout, "match: %d\n", n_match);
//     fprintf(stdout, "  sim: %lf\n", (double) n_match / aln->length);

//     needleman_wunsch_free(nw);
//     alignment_free(aln);

//     return 0;
// }

#include "minimap.h"
#include "mmpriv.h"
int asm_is_overlap_spurious(Longread_List_t &lr_list, uint32_t qid, uint32_t tid, mm_idx_t *mi, mm_mapopt_t *mopt, mm_tbuf_t *tbuf, int t_rev, int ov_type)
{
    int par_min_aln_len = 3000;
    int par_max_end_gap = 50;
    //   t_rev: whether forward strand (0) or reverse strand (1) of target should overlap with the query
    // ov_type: +1 =======              -1 =======
    //                =======>          =======> 
    // fprintf(stdout, "in function asm_is_overlap_spurious()\n");
    int32_t tlen = mi->seq[0].len;
    int32_t qlen = lr_list.reads[qid].len;
    string qseq = get_uncompressed_dna(lr_list.reads[qid].comp_seq, qlen, lr_list.reads[qid].comp_len);
    // 
    // map the query
    mm_reg1_t *regs;
    int n_regs;
    regs = mm_map(mi, qlen, qseq.c_str(), &n_regs, tbuf, mopt, 0); // get all hits for the query
    // for (int i = 0; i < n_regs; i++) // traverse hits and print them out
    // {
    //     mm_reg1_t *r = &regs[i];
    //     fprintf(stdout, "%d\t%u\t%d\t%d\t%c\t", qid, qlen, r->qs, r->qe, "+-"[r->rev]);
    //     fprintf(stdout, "%d\t%u\t%d\t%d\t%d\t%d\t%d\tcg:Z:", tid, tlen, r->rs, r->re, r->mlen, r->blen, r->mapq);
    //     if(mopt->flag & MM_F_CIGAR)
    //     {
    //         assert(r->p); // with MM_F_CIGAR, this should not be NULL
    //         for (uint32_t z = 0; z < r->p->n_cigar; z++) // IMPORTANT: this gives the CIGAR in the aligned regions. NO soft/hard clippings!
    //             fprintf(stdout, "%d%c", r->p->cigar[z]>>4, "MIDNSH"[r->p->cigar[z]&0xf]);
    //         free(r->p);
    //     }
    //     fprintf(stdout, "\ts0:%d\n", r->score0);
    // }
    // return 0;
    int is_spurious = 1; // yes unless otherwise is proven below
    if(n_regs > 0)
    {
        mm_reg1_t *r = &regs[0];
        if(r->blen >= par_min_aln_len)
        {
            if(t_rev == 0 && ov_type > 0)
            {
                if(tlen - r->re < par_max_end_gap && r->qs < par_max_end_gap)
                {
                    is_spurious = 0;
                }
            }
            else if(t_rev == 0 && ov_type < 0)
            {
                if(r->rs < par_max_end_gap && qlen - r->qe < par_max_end_gap)
                {
                    is_spurious = 0;
                }
            }
            else if(t_rev == 1 && ov_type > 0)
            {
                if(r->rs < par_max_end_gap && r->qs < par_max_end_gap)
                {
                    is_spurious = 0;
                }
            }
            else // if(t_rev == 1 && ov_type < 0)
            {
                if(tlen - r->re < par_max_end_gap && qlen - r->qe < par_max_end_gap)
                {
                    is_spurious = 0;
                }
            }
        }
    }
    free(regs);
    // fprintf(stdout, "return value: %d\n", is_spurious);
    return is_spurious;
}

int asm_is_alignment_ok(vector<vector<Align_Seq_t*> > &compact_lr_list, uint32_t lr1, uint32_t lr2, deque<int32_t> &aln1, deque<int32_t> &aln2, int t_rev)
{
    //   t_rev: whether forward strand (0) or reverse strand (1) of target should align with the query
    int gap_bp = 0;
    vector<Align_Seq_t*> &lr1_comp = compact_lr_list[lr1];
    vector<Align_Seq_t*> &lr2_comp = compact_lr_list[lr2];

    int last = -1;
    int gap_start = -1;
    for(int32_t i = 0; i < (int32_t)aln1.size(); i++)
    {
        if(aln1[i] != -1) // match
        {
            if(last == -1 && gap_start != -1)
            {
                // calculate number of bases in the gap
                gap_bp += lr2_comp[aln2[i-1]]->q_end - lr2_comp[aln2[gap_start]]->q_start;
            }
            gap_start = -1;
        }
        else // gap
        {
            if(last != -1)
            {
                gap_start = i;
            }
        }
        last = aln1[i];
    }
    // fprintf(stdout, "in function asm_is_alignment_ok() gap_bp:  %d\n", gap_bp);
    // 
    last = -1;
    gap_start = -1;
    for(int32_t i = 0; i < (int32_t)aln2.size(); i++)
    {
        if(aln2[i] != -1) // match
        {
            if(last == -1 && gap_start != -1)
            {
                // calculate number of bases in the gap
                if(t_rev)
                    gap_bp += lr1_comp[aln1[gap_start]]->q_end - lr1_comp[aln1[i-1]]->q_start;
                else
                    gap_bp += lr1_comp[aln1[i-1]]->q_end - lr1_comp[aln1[gap_start]]->q_start;
            }
            gap_start = -1;
        }
        else // gap
        {
            if(last != -1)
            {
                gap_start = i;
            }
        }
        last = aln2[i];
    }

    // fprintf(stdout, "in function asm_is_alignment_ok() gap_bp:  %d\n", gap_bp);
    return gap_bp;
}

// int asm_is_lr_contained(uint32_t qid, uint32_t tid, mm_idx_t *mi, mm_mapopt_t *mopt, mm_tbuf_t *tbuf, int check_target)
// {
//     // int par_min_aln_len = 3000;
//     int par_max_end_gap = 50;
//     // check_target: is target sequence contained or the query sequence
//     fprintf(stdout, "in function asm_is_lr_contained() check_target:%d\n", check_target);
//     int32_t tlen = mi->seq[0].len;
//     int32_t qlen = _asm_lr_list->reads[qid].len;
//     string qseq = get_uncompressed_dna(_asm_lr_list->reads[qid].comp_seq, qlen, _asm_lr_list->reads[qid].comp_len);
//     // 
//     // map the query
//     mm_reg1_t *regs;
//     int n_regs;
//     regs = mm_map(mi, qlen, qseq.c_str(), &n_regs, tbuf, mopt, 0); // get all hits for the query
//     for (int i = 0; i < n_regs; i++) // traverse hits and print them out
//     {
//         mm_reg1_t *r = &regs[i];
//         fprintf(stdout, "%d\t%u\t%d\t%d\t%c\t", qid, qlen, r->qs, r->qe, "+-"[r->rev]);
//         fprintf(stdout, "%d\t%u\t%d\t%d\t%d\t%d\t%d\tcg:Z:", tid, tlen, r->rs, r->re, r->mlen, r->blen, r->mapq);
//         if(mopt->flag & MM_F_CIGAR)
//         {
//             assert(r->p); // with MM_F_CIGAR, this should not be NULL
//             for (uint32_t z = 0; z < r->p->n_cigar; z++) // IMPORTANT: this gives the CIGAR in the aligned regions. NO soft/hard clippings!
//                 fprintf(stdout, "%d%c", r->p->cigar[z]>>4, "MIDNSH"[r->p->cigar[z]&0xf]);
//             free(r->p);
//         }
//         fprintf(stdout, "\n");
//     }
//     // return 0;
//     int is_contained = 0; // no unless otherwise is proven below
//     if(n_regs > 0)
//     {
//         mm_reg1_t *r = &regs[0];
//         if(check_target)
//         {
//             // if(r->blen >= par_min_aln_len && r->rs < par_max_end_gap && tlen - r->re < par_max_end_gap)
//             if(r->rs < par_max_end_gap && tlen - r->re < par_max_end_gap)
//             {
//                 is_contained = 1;
//             }
//         }
//         else
//         {
//             // if(r->blen >= par_min_aln_len && r->qs < par_max_end_gap && qlen - r->qe < par_max_end_gap)
//             if(r->qs < par_max_end_gap && qlen - r->qe < par_max_end_gap)
//             {
//                 is_contained = 1;
//             }
//         }
//     }
//     free(regs);
//     return is_contained;
// }

void asm_get_overlap_type(deque<int32_t> &aln1, deque<int32_t> &aln2, int &front_type, int &back_type)
{
    // check front
    front_type = 0;
    if(aln1.front() == -1)
        front_type = -1; // before
    else if(aln2.front() == -1)
        front_type = +1; // after
    // check back
    back_type = 0;
    if(aln1.back() == -1)
        back_type = +1; // after
    else if(aln2.back() == -1)
        back_type = -1; // before
}

void asm_ovgrpah_add_edge(vector<Ov_Node_t> &ov_graph, uint32_t lr1, uint8_t rev1, uint32_t lr2, uint8_t rev2, deque<int32_t> &aln1, deque<int32_t> &aln2, deque<int32_t> &aln1_invert, deque<int32_t> &aln2_invert)
{
    if(ov_graph[lr1].is_contained || ov_graph[lr2].is_contained)
        return;
    // 
    uint32_t to_node1, to_node2;
    if(rev1 == 0 && rev2 == 0)
    {
        // add edge +:n:+
        to_node1 = (lr2 << 1) | 0;
        ov_graph[lr1].out[to_node1].is_trasitive = 0;
        ov_graph[lr1].out[to_node1].aln_rev1 = 0;
        ov_graph[lr1].out[to_node1].aln_rev2 = 0;
        ov_graph[lr1].out[to_node1].aln_cmp1 = aln1;
        ov_graph[lr1].out[to_node1].aln_cmp2 = aln2;
        // add twin -:n:-
        to_node2 = (lr1 << 1) | 1;
        ov_graph[lr2].out_rev[to_node2].is_trasitive = 0;
        ov_graph[lr2].out_rev[to_node2].aln_rev1 = 1;
        ov_graph[lr2].out_rev[to_node2].aln_rev2 = 1;
        ov_graph[lr2].out_rev[to_node2].aln_cmp1 = aln2_invert;
        ov_graph[lr2].out_rev[to_node2].aln_cmp2 = aln1_invert;
    }
    else if(rev1 == 0 && rev2 == 1)
    {
        // add edge +:n:-
        to_node1 = (lr2 << 1) | 1;
        ov_graph[lr1].out[to_node1].is_trasitive = 0;
        ov_graph[lr1].out[to_node1].aln_rev1 = 0;
        ov_graph[lr1].out[to_node1].aln_rev2 = 1;
        ov_graph[lr1].out[to_node1].aln_cmp1 = aln1;
        ov_graph[lr1].out[to_node1].aln_cmp2 = aln2;
        // add twin +:n:-
        to_node2 = (lr1 << 1) | 1;
        ov_graph[lr2].out[to_node2].is_trasitive = 0;
        ov_graph[lr2].out[to_node2].aln_rev1 = 0;
        ov_graph[lr2].out[to_node2].aln_rev2 = 1;
        ov_graph[lr2].out[to_node2].aln_cmp1 = aln2_invert;
        ov_graph[lr2].out[to_node2].aln_cmp2 = aln1_invert;
    }
    else if(rev1 == 1 && rev2 == 0)
    {
        // add edge -:n:+
        to_node1 = (lr2 << 1) | 0;
        ov_graph[lr1].out_rev[to_node1].is_trasitive = 0;
        ov_graph[lr1].out_rev[to_node1].aln_rev1 = 1;
        ov_graph[lr1].out_rev[to_node1].aln_rev2 = 0;
        ov_graph[lr1].out_rev[to_node1].aln_cmp1 = aln1;
        ov_graph[lr1].out_rev[to_node1].aln_cmp2 = aln2;
        // add twin -:n:+
        to_node2 = (lr1 << 1) | 0;
        ov_graph[lr2].out_rev[to_node2].is_trasitive = 0;
        ov_graph[lr2].out_rev[to_node2].aln_rev1 = 1;
        ov_graph[lr2].out_rev[to_node2].aln_rev2 = 0;
        ov_graph[lr2].out_rev[to_node2].aln_cmp1 = aln2_invert;
        ov_graph[lr2].out_rev[to_node2].aln_cmp2 = aln1_invert;
    }
    else // if(rev1 == 1 && rev2 == 1)
    {
        // add edge -:n:-
        to_node1 = (lr2 << 1) | 1;
        ov_graph[lr1].out_rev[to_node1].is_trasitive = 0;
        ov_graph[lr1].out_rev[to_node1].aln_rev1 = 1;
        ov_graph[lr1].out_rev[to_node1].aln_rev2 = 1;
        ov_graph[lr1].out_rev[to_node1].aln_cmp1 = aln1;
        ov_graph[lr1].out_rev[to_node1].aln_cmp2 = aln2;
        // add twin +:n:+
        to_node2 = (lr1 << 1) | 0;
        ov_graph[lr2].out[to_node2].is_trasitive = 0;
        ov_graph[lr2].out[to_node2].aln_rev1 = 0;
        ov_graph[lr2].out[to_node2].aln_rev2 = 0;
        ov_graph[lr2].out[to_node2].aln_cmp1 = aln2_invert;
        ov_graph[lr2].out[to_node2].aln_cmp2 = aln1_invert;
    }
}

void asm_ovgraph_print_reduced_gfa(vector<Ov_Node_t> &ov_graph, string path)
{
    // print the reduced graph
    FILE *fp = fopen(path.c_str(), "w");
    if(fp == NULL)
    {
        fprintf(stderr, "[ERROR] cannot open file: %s\n", path.c_str());
        exit(EXIT_FAILURE);
    }
    // print only nodes that appear on some edges
    uint32_t i;
    set<uint32_t> to_print;
    for(i = 0; i < ov_graph.size(); i++)
    {
        // fprintf(stdout, "[check] out: %zu\n", ov_graph[i].out.size());
        if(ov_graph[i].out.size() > 0)
        {
            for(map<uint32_t, Ov_edge_t>::iterator it = ov_graph[i].out.begin(); it != ov_graph[i].out.end(); it++)
            {
                // fprintf(stdout, "[check] adding out %u\n", i);
                to_print.insert(i);
                to_print.insert((it->first >> 1));
            }
        }
        // fprintf(stdout, "[check] out_rev: %zu\n", ov_graph[i].out_rev.size());
        if(ov_graph[i].out_rev.size() > 0)
        {
            for(map<uint32_t, Ov_edge_t>::iterator it = ov_graph[i].out_rev.begin(); it != ov_graph[i].out_rev.end(); it++)
            {
                // fprintf(stdout, "[check] adding out_rev %u\n", i);
                to_print.insert(i);
                to_print.insert((it->first >> 1));
            }
        }
    }
    for(set<uint32_t>::iterator it = to_print.begin(); it != to_print.end(); it++)
    {
        fprintf(fp, "S\t%u\t%s\n", *it, string(100, 'N').c_str());
    }

    // print links
    for(i = 0; i < ov_graph.size(); i++)
    {
        if(ov_graph[i].out.size() > 0)
        {
            for(map<uint32_t, Ov_edge_t>::iterator it = ov_graph[i].out.begin(); it != ov_graph[i].out.end(); it++)
            {
                fprintf(fp, "L\t%u\t+\t%u\t%c\t0M\n", i, (it->first >> 1), ((it->first & 1) ? '-' : '+'));
            }
        }
        if(ov_graph[i].out_rev.size() > 0)
        {
            for(map<uint32_t, Ov_edge_t>::iterator it = ov_graph[i].out_rev.begin(); it != ov_graph[i].out_rev.end(); it++)
            {
                fprintf(fp, "L\t%u\t-\t%u\t%c\t0M\n", i, (it->first >> 1), ((it->first & 1) ? '-' : '+'));
            }
        }
    }

    fclose(fp);
}

// prints the overlap graph in GFA format
void asm_ovgraph_print_gfa(vector<Ov_Node_t> &ov_graph, string path)
{
    FILE *fp = fopen(path.c_str(), "w");
    if(fp == NULL)
    {
        fprintf(stderr, "[ERROR] cannot open file: %s\n", path.c_str());
        exit(EXIT_FAILURE);
    }

    // print only nodes that appear on some edges
    set<uint32_t> to_print;
    for(uint32_t i = 0; i < ov_graph.size(); i++)
    {
        map<uint32_t, Ov_edge_t>::iterator it;
        if(ov_graph[i].out.size() > 0)
        {
            for(it = ov_graph[i].out.begin(); it != ov_graph[i].out.end(); it++)
            {
                if(it->second.is_trasitive == 0)
                {
                    if(ov_graph[i].is_contained == false)
                        to_print.insert(i);
                    if(ov_graph[(it->first >> 1)].is_contained == false)
                        to_print.insert((it->first >> 1));
                }
            }
        }
        if(ov_graph[i].out_rev.size() > 0)
        {
            for(it = ov_graph[i].out_rev.begin(); it != ov_graph[i].out_rev.end(); it++)
            {
                if(it->second.is_trasitive == 0)
                {
                    if(ov_graph[i].is_contained == false)
                        to_print.insert(i);
                    if(ov_graph[(it->first >> 1)].is_contained == false)
                        to_print.insert((it->first >> 1));
                }
            }
        }
    }
    for(set<uint32_t>::iterator it = to_print.begin(); it != to_print.end(); it++)
    {
        fprintf(fp, "S\t%u\t%s\n", *it, string(100, 'N').c_str());
    }

    // print links
    for(uint32_t i = 0; i < ov_graph.size(); i++)
    {
        map<uint32_t, Ov_edge_t>::iterator it;
        if(ov_graph[i].out.size() > 0)
        {
            for(it = ov_graph[i].out.begin(); it != ov_graph[i].out.end(); it++)
            {
                if(it->second.is_trasitive == 0 && ov_graph[i].is_contained == false && ov_graph[(it->first >> 1)].is_contained == false)
                    fprintf(fp, "L\t%u\t+\t%u\t%c\t0M\n", i, (it->first >> 1), ((it->first & 1) ? '-' : '+'));
            }
        }
        if(ov_graph[i].out_rev.size() > 0)
        {
            for(it = ov_graph[i].out_rev.begin(); it != ov_graph[i].out_rev.end(); it++)
            {
                if(it->second.is_trasitive == 0 && ov_graph[i].is_contained == false && ov_graph[(it->first >> 1)].is_contained == false)
                    fprintf(fp, "L\t%u\t-\t%u\t%c\t0M\n", i, (it->first >> 1), ((it->first & 1) ? '-' : '+'));
            }
        }
    }

    fclose(fp);
}

void asm_ovgraph_transitive_reduction(vector<Ov_Node_t> &ov_graph)
{
    uint32_t i;
    map<uint32_t, Ov_edge_t>::iterator it_incoming, it_outgoing;
    for(i = 0; i < ov_graph.size(); i++)
    {
        for(it_incoming = ov_graph[i].out_rev.begin(); it_incoming != ov_graph[i].out_rev.end(); it_incoming++)
        {
            for(it_outgoing = ov_graph[i].out.begin(); it_outgoing != ov_graph[i].out.end(); it_outgoing++)
            {
                // NOTE: I'm using '+' and '-' only to simplify things
                // a -> b -> c
                uint32_t a_id;
                // char a_from, a_to;
                char a_from;
                a_id   = (it_incoming->first >> 1);
                a_from = ((it_incoming->first & 1) ? '+' : '-');
                // a_to   = '+';
                uint32_t c_id;
                // char c_from, c_to;
                char c_to;
                c_id   = (it_outgoing->first >> 1);
                // c_from = '+';
                c_to   = ((it_outgoing->first & 1) ? '-' : '+');
                // check edge a -> c
                uint32_t to_node1, to_node2;
                if(a_from == '+' && c_to == '+')
                {
                    to_node1 = (c_id << 1) | 0;
                    to_node2 = (a_id << 1) | 1; // since a_from is '+' it's twin is '-'
                    if(ov_graph[a_id].out.count(to_node1) > 0) // transitive edge exists
                        ov_graph[a_id].out[to_node1].is_trasitive = 1; // mark as transitive
                    if(ov_graph[c_id].out_rev.count(to_node2) > 0) // transitive edge exists
                        ov_graph[c_id].out_rev[to_node2].is_trasitive = 1; // mark as transitive
                }
                else if(a_from == '+' && c_to == '-')
                {
                    to_node1 = (c_id << 1) | 1;
                    to_node2 = (a_id << 1) | 1; // since a_from is '+' it's twin is '-'
                    if(ov_graph[a_id].out.count(to_node1) > 0) // transitive edge exists
                        ov_graph[a_id].out[to_node1].is_trasitive = 1; // mark as transitive
                    if(ov_graph[c_id].out.count(to_node2) > 0) // transitive edge exists
                        ov_graph[c_id].out[to_node2].is_trasitive = 1; // mark as transitive
                }
                else if(a_from == '-' && c_to == '+')
                {
                    to_node1 = (c_id << 1) | 0;
                    to_node2 = (a_id << 1) | 0; // since a_from is '-' it's twin is '+'
                    if(ov_graph[a_id].out_rev.count(to_node1) > 0) // transitive edge exists
                        ov_graph[a_id].out_rev[to_node1].is_trasitive = 1; // mark as transitive
                    if(ov_graph[c_id].out_rev.count(to_node2) > 0) // transitive edge exists
                        ov_graph[c_id].out_rev[to_node2].is_trasitive = 1; // mark as transitive
                }
                else // if(a_from == '-' && c_to == '-')
                {
                    to_node1 = (c_id << 1) | 1;
                    to_node2 = (a_id << 1) | 0; // since a_from is '-' it's twin is '+'
                    if(ov_graph[a_id].out_rev.count(to_node1) > 0) // transitive edge exists
                        ov_graph[a_id].out_rev[to_node1].is_trasitive = 1; // mark as transitive
                    if(ov_graph[c_id].out.count(to_node2) > 0) // transitive edge exists
                        ov_graph[c_id].out[to_node2].is_trasitive = 1; // mark as transitive
                }
            }
        }
    }
}

// removes transitive edges as well as edges connected to the contained reads
void asm_ovgraph_remove_transitive(vector<Ov_Node_t> &ov_graph)
{
    uint32_t i;
    for(i = 0; i < ov_graph.size(); i++)
    {
        if(ov_graph[i].is_contained == true)
        {
            ov_graph[i].out.clear();
            ov_graph[i].out_rev.clear();
        }
        else
        {
            uint32_t node2, to_node2;
            uint8_t rev2;
            for(map<uint32_t, Ov_edge_t>::iterator it = ov_graph[i].out.begin(); it != ov_graph[i].out.end(); )
            {
                if(it->second.is_trasitive == 1 || ov_graph[(it->first >> 1)].is_contained == true)
                {
                    node2 = (it->first >> 1);
                    rev2 = (it->first & 1);
                    // remove edges
                    to_node2 = (i << 1) | 1;
                    it = ov_graph[i].out.erase(it);
                    if(rev2 == 0)
                    {
                        ov_graph[node2].out_rev.erase(to_node2); // NOTE: here we are removing the dual edge
                    }
                    else
                    {
                        ov_graph[node2].out.erase(to_node2); // NOTE: here we are removing the dual edge
                    }
                }
                else
                {
                    it++;
                }
            }
            // 
            for(map<uint32_t, Ov_edge_t>::iterator it = ov_graph[i].out_rev.begin(); it != ov_graph[i].out_rev.end(); )
            {
                if(it->second.is_trasitive == 1 || ov_graph[(it->first >> 1)].is_contained == true)
                {
                    node2 = (it->first >> 1);
                    rev2 = (it->first & 1);
                    // remove edges
                    to_node2 = (i << 1) | 0;
                    it = ov_graph[i].out_rev.erase(it);
                    if(rev2 == 0)
                    {
                        ov_graph[node2].out_rev.erase(to_node2); // NOTE: here we are removing the dual edge
                    }
                    else
                    {
                        ov_graph[node2].out.erase(to_node2); // NOTE: here we are removing the dual edge
                    }
                }
                else
                {
                    it++;
                }
            }
        }
    }
}

bool asm_ovgraph_find_next_node(vector<Ov_Node_t> &ov_graph, uint32_t curr_node, uint32_t curr_strand, uint32_t &next_node, uint32_t &next_strand)
{
    if(curr_strand == 0)
    {
        if(ov_graph[curr_node].out.size() == 1)
        {
            next_node = (ov_graph[curr_node].out.begin()->first >> 1);
            next_strand = (ov_graph[curr_node].out.begin()->first & 1);
            return true;
        }
        else
        {
            return false;
        }
    }
    else // curr_strand == 1
    {
        if(ov_graph[curr_node].out_rev.size() == 1)
        {
            next_node = (ov_graph[curr_node].out_rev.begin()->first >> 1);
            next_strand = (ov_graph[curr_node].out_rev.begin()->first & 1);
            return true;
        }
        else
        {
            return false;
        }
    }
}

void asm_ovgraph_get_paths(vector<Ov_Node_t> &ov_graph, vector<uint8_t> &unused, vector<vector<Align_Seq_t*> > &compact_lr_list, vector<vector<path_elem_t>> &ov_paths, vector<vector<uint32_t>> &covering_lr)
{
    fprintf(stdout, "in function asm_ovgraph_get_paths()\n");
    // for(uint32_t i = 0; i < ov_graph.size(); i++)
    // {
    //     fprintf(stdout, "node %u\n", i);
    //     if(unused[i])
    //     {
    //         if(ov_graph[i].out.size() == 0 && ov_graph[i].out_rev.size() == 0)
    //         {
    //             if(ov_graph[i].is_contained == 0)
    //                 fprintf(stdout, "\tsingletone to be printed\n");
    //         }
    //         else
    //         {
    //             map<uint32_t, uint8_t>::iterator it;
    //             if(ov_graph[i].out.size() > 0)
    //             {
    //                 for(it = ov_graph[i].out.begin(); it != ov_graph[i].out.end(); it++)
    //                 {
    //                     if(it->second == 0 && ov_graph[i].is_contained == 0 && ov_graph[(it->first >> 1)].is_contained == 0)
    //                     {
    //                         fprintf(stdout, "\tto_check1 %u\n", i);
    //                         fprintf(stdout, "\tto_check2 %u\n", (it->first >> 1));
    //                     }
    //                 }
    //             }
    //             if(ov_graph[i].out_rev.size() > 0)
    //             {
    //                 for(it = ov_graph[i].out_rev.begin(); it != ov_graph[i].out_rev.end(); it++)
    //                 {
    //                     if(it->second == 0 && ov_graph[i].is_contained == 0 && ov_graph[(it->first >> 1)].is_contained == 0)
    //                     {
    //                         fprintf(stdout, "\tto_check1 %u\n", i);
    //                         fprintf(stdout, "\tto_check2 %u\n", (it->first >> 1));
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }
    uint32_t i, j;
    // bool is_singleton;
    bool traverse_outgoing;
    bool traverse_incoming;
    vector<bool> visited(ov_graph.size(), false);
    while(true)
    {
        traverse_outgoing = false;
        traverse_incoming = false;
        for(i = 0; i < ov_graph.size(); i++)
        {
            // TODO: what if out-degree is more than 1? should I traverse?
            if (ov_graph[i].out.size() == 1 && ov_graph[i].out_rev.size() == 0 && visited[i] == false) // node is a tip (leaf) with outgoing edges
            {
                traverse_outgoing = true;
                break;
            }
            // TODO: what if in-degree is more than 1? should I traverse?
            else if (ov_graph[i].out.size() == 0 && ov_graph[i].out_rev.size() == 1 && visited[i] == false) // node is a tip (leaf) with incoming edges
            {
                traverse_incoming = true;
                break;
            }
            else if (ov_graph[i].out.size() == 1 && ov_graph[i].out_rev.size() == 1 && visited[i] == false) // node is on a simple path
            {
                traverse_outgoing = true;
                traverse_incoming = true;
                break;
            }
        }
        if(i == ov_graph.size()) // no other node to start traversing from
            break;
        fprintf(stdout, "new_traversal node:%u outgoing:%zu incoming:%zu\n", i, ov_graph[i].out.size(), ov_graph[i].out_rev.size());
        visited[i] = true;
        vector<id_strand_t> path_outgoing;
        if(traverse_outgoing == true)
        {
            uint32_t curr_node = i;
            uint32_t curr_strand = 0;
            uint32_t next_node, next_strand;
            bool found_next;
            while(true)
            {
                found_next = asm_ovgraph_find_next_node(ov_graph, curr_node, curr_strand, next_node, next_strand);
                if(found_next == false || visited[next_node] == true)
                    break;
                curr_node = next_node;
                curr_strand = next_strand;
                visited[curr_node] = true;
                path_outgoing.push_back({curr_strand, 0, curr_node});
            }
        }
        // fprintf(stdout, "==== path_outgoing ");
        // for(size_t j = 0; j < path_outgoing.size(); j++)
        //     fprintf(stdout, "%u:%c ", path_outgoing[j].id, "+-"[path_outgoing[j].strand]);
        // fprintf(stdout, "\n");
        // 
        vector<id_strand_t> path_incoming;
        if(traverse_incoming == true)
        {
            uint32_t curr_node = i;
            uint32_t curr_strand = 1;
            uint32_t next_node, next_strand;
            bool found_next;
            while(true)
            {
                found_next = asm_ovgraph_find_next_node(ov_graph, curr_node, curr_strand, next_node, next_strand);
                if(found_next == false || visited[next_node] == true)
                    break;
                curr_node = next_node;
                curr_strand = next_strand;
                visited[curr_node] = true;
                path_incoming.push_back({curr_strand, 0, curr_node});
            }
        }
        // fprintf(stdout, "==== path_incoming ");
        // for(size_t j = 0; j < path_incoming.size(); j++)
        //     fprintf(stdout, "%u:%c ", path_incoming[j].id, "+-"[path_incoming[j].strand]);
        // fprintf(stdout, "\n");
        // 
        vector<id_strand_t> path;
        // add the incoming path
        for (vector<id_strand_t>::reverse_iterator it = path_incoming.rbegin(); it != path_incoming.rend(); ++it)
            path.push_back({(uint32_t)1 - it->strand, 0, it->id});
        // add the node
        path.push_back({0, 0, i});
        // add the outgoing path
        for (vector<id_strand_t>::iterator it = path_outgoing.begin(); it != path_outgoing.end(); ++it)
            path.push_back({it->strand, 0, it->id});
        // print the path
        if(path.size() > 1)
            fprintf(stdout, "==== path %u:%c -> %u:%c\n", path.front().id, "+-"[path.front().strand], path.back().id, "+-"[path.back().strand]);
        else
            fprintf(stdout, "==== path %u:%c\n", path.front().id, "+-"[path.front().strand]);
        fprintf(stdout, "==== ");
        for(j = 0; j < path.size(); j++)
            fprintf(stdout, "%u:%c ", path[j].id, "+-"[path[j].strand]);
        fprintf(stdout, "\n");
        fflush(stdout);
        // // 
        // if(path.size() > 0)
        //     _asm_all_paths.push_back(path);
        if(path.size() < 2) continue;

        uint32_t p_node, p_strand;
        uint32_t p_node2, p_strand2;
        deque<int32_t> *aln1 = NULL;
        deque<int32_t> *aln2 = NULL;
        uint8_t aln1_rev = 0, aln2_rev = 0;
        vector<path_elem_t> sr_contig_path;
        int32_t start_index = 0;
        int32_t aln_index = 0;
        int32_t last_non_gap = -1;
        for(j = 0; j < path.size() - 1; j++)
        {
            p_node = path[j].id;
            p_strand = path[j].strand;
            p_node2 = path[j+1].id;
            p_strand2 = path[j+1].strand;
            fprintf(stdout, "%u:%c\n", p_node, "+-"[p_strand]);
            fflush(stdout);
            if(p_strand == 0)
            {
                // fprintf(stdout, "next %u:%c (# out: %zu)\n", (ov_graph[p_node].out.begin()->first >> 1), "+-"[ov_graph[p_node].out.begin()->first & 1], ov_graph[p_node].out.size());
                fprintf(stdout, "next %u:%c (# out: %zu)\n", p_node2, "+-"[p_strand2], ov_graph[p_node].out.size());
                fflush(stdout);
                uint32_t node_key = (p_node2 << 1) | p_strand2;
                aln1 = &(ov_graph[p_node].out[node_key].aln_cmp1);
                aln1_rev = ov_graph[p_node].out[node_key].aln_rev1;
                aln2 = &(ov_graph[p_node].out[node_key].aln_cmp2);
                aln2_rev = ov_graph[p_node].out[node_key].aln_rev2;
            }
            else
            {
                // fprintf(stdout, "next %u:%c (# out_rev: %zu)\n", (ov_graph[p_node].out_rev.begin()->first >> 1), "+-"[ov_graph[p_node].out_rev.begin()->first & 1], ov_graph[p_node].out_rev.size());
                fprintf(stdout, "next %u:%c (# out_rev: %zu)\n", p_node2, "+-"[p_strand2], ov_graph[p_node].out_rev.size());
                fflush(stdout);
                uint32_t node_key = (p_node2 << 1) | p_strand2;
                aln1 = &(ov_graph[p_node].out_rev[node_key].aln_cmp1);
                aln1_rev = ov_graph[p_node].out_rev[node_key].aln_rev1;
                aln2 = &(ov_graph[p_node].out_rev[node_key].aln_cmp2);
                aln2_rev = ov_graph[p_node].out_rev[node_key].aln_rev2;
            }
            for(size_t z = 0; z < aln1->size(); z++)
            {
                int32_t ind = aln1->operator[](z);
                fprintf(stdout, "%8d ", ind);
            }
            fprintf(stdout, "\n");
            for(size_t z = 0; z < aln1->size(); z++)
            {
                int32_t ind = aln1->operator[](z);
                if(ind == -1)
                {
                    fprintf(stdout, "%8s ", "---");
                }
                else
                {
                    fprintf(stdout, "%6d:%c ", compact_lr_list[p_node][ind]->t_id, "+-"[compact_lr_list[p_node][ind]->is_rev ^ aln1_rev]);
                }
            }
            fprintf(stdout, "\n");
            for(size_t z = 0; z < aln2->size(); z++)
            {
                int32_t ind = aln2->operator[](z);
                if(ind == -1)
                {
                    fprintf(stdout, "%8s ", "---");
                }
                else
                {
                    fprintf(stdout, "%6d:%c ", compact_lr_list[p_node2][ind]->t_id, "+-"[compact_lr_list[p_node2][ind]->is_rev ^ aln2_rev]);
                }
            }
            fprintf(stdout, "\n");
            for(size_t z = 0; z < aln2->size(); z++)
            {
                int32_t ind = aln2->operator[](z);
                fprintf(stdout, "%8d ", ind);
            }
            fprintf(stdout, "\n");
            fflush(stdout);
            // 
            if(j == 0)
            {
                aln_index = 0;
            }
            else
            {
                for(size_t z = 0; z < aln1->size(); z++)
                {
                    if(aln1->operator[](z) == start_index)
                        aln_index = z;
                }
            }
            // 
            for(size_t z = aln_index; z < aln1->size(); z++)
            {
                int32_t ind = aln1->operator[](z);
                if(ind != -1)
                {
                    last_non_gap = z;
                    sr_contig_path.push_back({uint32_t(compact_lr_list[p_node][ind]->is_rev ^ aln1_rev), compact_lr_list[p_node][ind]->t_id, set<uint32_t>(), map<uint32_t, int32_t>()});
                }
            }
            while(aln2->operator[](++last_non_gap) == -1);
            start_index = aln2->operator[](last_non_gap);
        }
        p_node = path[j].id;
        p_strand = path[j].strand;
        fprintf(stdout, "%u:%c\n", p_node, "+-"[p_strand]);
        fprintf(stdout, "\n");
        fflush(stdout);
        // 
        for(size_t z = last_non_gap; z < aln2->size(); z++)
        {
            int32_t ind = aln2->operator[](z);
            if(ind != -1)
            {
                sr_contig_path.push_back({uint32_t(compact_lr_list[p_node][ind]->is_rev ^ aln2_rev), compact_lr_list[p_node][ind]->t_id, set<uint32_t>(), map<uint32_t, int32_t>()});
            }
        }
        fprintf(stdout, "#### ");
        for(uint32_t z = 0; z < sr_contig_path.size(); z++)
        {
            fprintf(stdout, "%u:%c ", sr_contig_path[z].id, "+-"[sr_contig_path[z].strand]);
        }
        fprintf(stdout, "\n\n");
        ov_paths.push_back(sr_contig_path);
        uint32_t ov_paths_size = ov_paths.size();
        covering_lr.resize(ov_paths_size);
        fprintf(stdout, "\n\nReads that fall in this path:\n");
        for(size_t z = 0; z < path.size(); z++)
        {
            covering_lr[ov_paths_size - 1].push_back(path[z].id);
            fprintf(stdout, "%u\t", path[z].id);
            for(size_t x = 0; x < ov_graph[path[z].id].containing_list.size(); x++)
            {
                fprintf(stdout, "%u\t", ov_graph[path[z].id].containing_list[x]);
                covering_lr[ov_paths_size - 1].push_back(ov_graph[path[z].id].containing_list[x]);
            }
        }
        fprintf(stdout, "\n\n");
    }
}

void asm_build_ovgraph_from_unused_lrs(Longread_List_t &lr_list, vector<vector<Align_Seq_t*>> &compact_lr_list, vector<vector<id_strand_t>> &map_contig2lr, vector<vector<path_elem_t>> &simple_paths, vector<vector<id_strand_t>> &map_contig2lr_uniq, vector<uint8_t> &unused, vector<vector<path_elem_t>> &ov_paths, vector<vector<uint32_t>> &covering_lr)
{
    uint32_t i, j;
    // build the overlap graph
    vector<Ov_Node_t> ov_graph;
    ov_graph.resize(compact_lr_list.size()); // FIXME: we can save on memory by allocating only the number of unused long reads. This requires extracting a list of unused long reads.
    for(i = 0; i < compact_lr_list.size(); i++)
        ov_graph[i].is_contained = 0;

    vector<uint32_t> lrs_to_align;
    for(i = 0; i < compact_lr_list.size(); i++)
    {
        if(unused[i]) lrs_to_align.push_back(i);
    }

    for(i = 0; i < lrs_to_align.size(); i++)
    {
        // 
        uint32_t target_id = lrs_to_align[i];
        string target_seq = get_uncompressed_dna(lr_list.reads[target_id].comp_seq, lr_list.reads[target_id].len, lr_list.reads[target_id].comp_len);
        mm_idxopt_t iopt;
        mm_mapopt_t mopt;
        mm_idxopt_init(&iopt);
        mm_mapopt_init(&mopt);
        mm_set_opt("ava-pb", &iopt, &mopt);
        iopt.flag |= MM_I_HPC;  // compress homopolymers
        iopt.k = 17;
        // iopt.w = 3;
        mopt.flag |= MM_F_CIGAR; // perform alignment
        mopt.flag &= ~MM_F_NO_DUAL; // 
        mopt.mid_occ_frac = 0.1;
        mopt.min_chain_score = 200;
        mopt.max_gap = 2000;
        // index the target
        const char *dummy_name = "NA";
        const char *dummy_seq  = target_seq.c_str();
        mm_idx_t *mi;
        mi = mm_idx_str(iopt.w, iopt.k, iopt.flag & MM_I_HPC, iopt.bucket_bits, 1, (const char**)&dummy_seq, (const char**)&dummy_name);
        mm_mapopt_update(&mopt, mi);
        mm_tbuf_t *tbuf = mm_tbuf_init();
        // 
        for(j = i+1; j < lrs_to_align.size(); j++)
        {
            fprintf(stdout, "\nLCS between long reads %u and %u\n", lrs_to_align[i], lrs_to_align[j]);
            // fflush(stdout);
            int lcs_score, lcs_score_rev;
            deque<int32_t> aln1, aln2; // contain alignment between lr1 and lr2
            deque<int32_t> aln1_rev, aln2_rev; // contain alignment between lr1_rev and lr2
            lcs_alignment(compact_lr_list, lrs_to_align[i], lrs_to_align[j], aln1, aln2, lcs_score);
            lcs_alignment_reverse(compact_lr_list, lrs_to_align[i], lrs_to_align[j], aln1_rev, aln2_rev, lcs_score_rev);

            if(lcs_score > lcs_score_rev)
            {
                if(lcs_score < 7)
                {
                    fprintf(stdout, "Reject: not enough matching contigs\n");
                    continue;
                }
                int front_type, back_type;
                asm_get_overlap_type(aln1, aln2, front_type, back_type);
                // 
                fprintf(stdout, "lr1 lr2 score: %d\n", lcs_score);
                fprintf(stdout, "front_type:%d\tback_type:%d\n", front_type, back_type);
                // 
                // for(uint32_t z = 0; z < aln1.size(); z++)
                //     if(aln1[z] == -1)
                //         fprintf(stdout, "%d\t", aln1[z]);
                //     else
                //         fprintf(stdout, "%u:%u\t", compact_lr_list[lrs_to_align[i]][aln1[z]]->t_id, compact_lr_list[lrs_to_align[i]][aln1[z]]->is_rev);
                // fprintf(stdout, "\n");
                // for(uint32_t z = 0; z < aln2.size(); z++)
                //     if(aln2[z] == -1)
                //         fprintf(stdout, "%d\t", aln2[z]);
                //     else
                //         fprintf(stdout, "%u:%u\t", compact_lr_list[lrs_to_align[j]][aln2[z]]->t_id, compact_lr_list[lrs_to_align[j]][aln2[z]]->is_rev);
                // fprintf(stdout, "\n");
                // 
                if(front_type > 0)
                {
                    if(back_type <= 0)
                    {
                        // if(asm_is_lr_contained(lrs_to_align[j], target_id, mi, &mopt, tbuf, 0))
                        if(asm_is_alignment_ok(compact_lr_list, lrs_to_align[i], lrs_to_align[j], aln1, aln2, 0) < 10000)
                        {
                            ov_graph[lrs_to_align[j]].is_contained = 1;
                            ov_graph[lrs_to_align[i]].containing_list.push_back(lrs_to_align[j]);
                            fprintf(stdout, "#### %u contained ####\n", lrs_to_align[j]);
                        }
                    }
                    else // back_type > 0
                    {
                        if(asm_is_overlap_spurious(lr_list, lrs_to_align[j], target_id, mi, &mopt, tbuf, 0, +1)) continue;
                        // deque<int32_t> cmp_aln1, cmp_aln2;
                        deque<int32_t> cmp_aln1_invert, cmp_aln2_invert;
                        cmp_aln1_invert.resize(aln1.size());
                        reverse_copy(aln1.begin(), aln1.end(), cmp_aln1_invert.begin());
                        cmp_aln2_invert.resize(aln2.size());
                        reverse_copy(aln2.begin(), aln2.end(), cmp_aln2_invert.begin());
                        // asm_ovgraph_get_cmp_alignments(compact_lr_list, lrs_to_align[i], lrs_to_align[j], aln1, aln2, false, cmp_aln1, cmp_aln2, cmp_aln1_invert, cmp_aln2_invert);
                        asm_ovgrpah_add_edge(ov_graph, lrs_to_align[i], 0, lrs_to_align[j], 0, aln1, aln2, cmp_aln1_invert, cmp_aln2_invert);
                    }
                }
                else if(front_type < 0)
                {
                    if(back_type >= 0)
                    {
                        // if(asm_is_lr_contained(lrs_to_align[j], target_id, mi, &mopt, tbuf, 1))
                        if(asm_is_alignment_ok(compact_lr_list, lrs_to_align[i], lrs_to_align[j], aln1, aln2, 0) < 10000)
                        {
                            ov_graph[lrs_to_align[i]].is_contained = 1;
                            ov_graph[lrs_to_align[j]].containing_list.push_back(lrs_to_align[i]);
                            fprintf(stdout, "#### %u contained ####\n", lrs_to_align[i]);
                        }
                    }
                    else // back_type < 0
                    {
                        if(asm_is_overlap_spurious(lr_list, lrs_to_align[j], target_id, mi, &mopt, tbuf, 0, -1)) continue;
                        // deque<cmp_aln_t> cmp_aln1, cmp_aln2;
                        deque<int32_t> cmp_aln1_invert, cmp_aln2_invert;
                        cmp_aln1_invert.resize(aln1.size());
                        reverse_copy(aln1.begin(), aln1.end(), cmp_aln1_invert.begin());
                        cmp_aln2_invert.resize(aln2.size());
                        reverse_copy(aln2.begin(), aln2.end(), cmp_aln2_invert.begin());
                        // asm_ovgraph_get_cmp_alignments(compact_lr_list, lrs_to_align[i], lrs_to_align[j], aln1, aln2, false, cmp_aln1, cmp_aln2, cmp_aln1_invert, cmp_aln2_invert);
                        asm_ovgrpah_add_edge(ov_graph, lrs_to_align[i], 1, lrs_to_align[j], 1, cmp_aln1_invert, cmp_aln2_invert, aln1, aln2);
                    }
                }
                else // front_type == 0
                {
                    if(back_type < 0)
                    {
                        // if(asm_is_lr_contained(lrs_to_align[j], target_id, mi, &mopt, tbuf, 0))
                        if(asm_is_alignment_ok(compact_lr_list, lrs_to_align[i], lrs_to_align[j], aln1, aln2, 0) < 10000)
                        {
                            ov_graph[lrs_to_align[j]].is_contained = 1;
                            ov_graph[lrs_to_align[i]].containing_list.push_back(lrs_to_align[j]);
                            fprintf(stdout, "#### %u contained ####\n", lrs_to_align[j]);
                        }
                    }
                    else if(back_type == 0)
                    {
                        if(lrs_to_align[i] < lrs_to_align[j])
                        {
                            // if(asm_is_lr_contained(lrs_to_align[j], target_id, mi, &mopt, tbuf, 0))
                            if(asm_is_alignment_ok(compact_lr_list, lrs_to_align[i], lrs_to_align[j], aln1, aln2, 0) < 10000)
                            {
                                ov_graph[lrs_to_align[j]].is_contained = 1;
                                ov_graph[lrs_to_align[i]].containing_list.push_back(lrs_to_align[j]);
                                fprintf(stdout, "#### %u contained ####\n", lrs_to_align[j]);
                            }
                        }
                        else
                        {
                            // if(asm_is_lr_contained(lrs_to_align[j], target_id, mi, &mopt, tbuf, 1))
                            if(asm_is_alignment_ok(compact_lr_list, lrs_to_align[i], lrs_to_align[j], aln1, aln2, 0) < 10000)
                            {
                                ov_graph[lrs_to_align[i]].is_contained = 1;
                                ov_graph[lrs_to_align[j]].containing_list.push_back(lrs_to_align[i]);
                                fprintf(stdout, "#### %u contained ####\n", lrs_to_align[i]);
                            }
                        }
                    }
                    else // back_type > 0
                    {
                        // if(asm_is_lr_contained(lrs_to_align[j], target_id, mi, &mopt, tbuf, 1))
                        if(asm_is_alignment_ok(compact_lr_list, lrs_to_align[i], lrs_to_align[j], aln1, aln2, 0) < 10000)
                        {
                            ov_graph[lrs_to_align[i]].is_contained = 1;
                            ov_graph[lrs_to_align[j]].containing_list.push_back(lrs_to_align[i]);
                            fprintf(stdout, "#### %u contained ####\n", lrs_to_align[i]);
                        }
                    }
                }
            }
            else
            {
                if(lcs_score_rev < 7)
                {
                    fprintf(stdout, "Reject: not enough matching contigs\n");
                    continue;
                }
                int front_type, back_type;
                asm_get_overlap_type(aln1_rev, aln2_rev, front_type, back_type);
                // 
                fprintf(stdout, "lr1_rev lr2 score: %d\n", lcs_score_rev);
                fprintf(stdout, "front_type:%d\tback_type:%d\n", front_type, back_type);
                // // 
                // for(uint32_t z = 0; z < aln1_rev.size(); z++)
                //     if(aln1_rev[z] == -1)
                //         fprintf(stdout, "%d\t", aln1_rev[z]);
                //     else
                //         fprintf(stdout, "%u:%u\t", compact_lr_list[lrs_to_align[i]][aln1_rev[z]]->t_id, 1 - compact_lr_list[lrs_to_align[i]][aln1_rev[z]]->is_rev);
                // fprintf(stdout, "\n");
                // for(uint32_t z = 0; z < aln2_rev.size(); z++)
                //     if(aln2_rev[z] == -1)
                //         fprintf(stdout, "%d\t", aln2_rev[z]);
                //     else
                //         fprintf(stdout, "%u:%u\t", compact_lr_list[lrs_to_align[j]][aln2_rev[z]]->t_id, compact_lr_list[lrs_to_align[j]][aln2_rev[z]]->is_rev);
                // fprintf(stdout, "\n");
                // 
                if(front_type > 0)
                {
                    if(back_type <= 0)
                    {
                        // if(asm_is_lr_contained(lrs_to_align[j], target_id, mi, &mopt, tbuf, 0))
                        if(asm_is_alignment_ok(compact_lr_list, lrs_to_align[i], lrs_to_align[j], aln1_rev, aln2_rev, 1) < 10000)
                        {
                            ov_graph[lrs_to_align[j]].is_contained = 1;
                            ov_graph[lrs_to_align[i]].containing_list.push_back(lrs_to_align[j]);
                            fprintf(stdout, "#### %u contained ####\n", lrs_to_align[j]);
                        }
                    }
                    else // back_type > 0
                    {
                        if(asm_is_overlap_spurious(lr_list, lrs_to_align[j], target_id, mi, &mopt, tbuf, 1, +1)) continue;
                        // deque<cmp_aln_t> cmp_aln1, cmp_aln2;
                        deque<int32_t> cmp_aln1_invert, cmp_aln2_invert;
                        cmp_aln1_invert.resize(aln1_rev.size());
                        reverse_copy(aln1_rev.begin(), aln1_rev.end(), cmp_aln1_invert.begin());
                        cmp_aln2_invert.resize(aln2_rev.size());
                        reverse_copy(aln2_rev.begin(), aln2_rev.end(), cmp_aln2_invert.begin());
                        // asm_ovgraph_get_cmp_alignments(compact_lr_list, lrs_to_align[i], lrs_to_align[j], aln1_rev, aln2_rev, true, cmp_aln1, cmp_aln2, cmp_aln1_invert, cmp_aln2_invert);
                        asm_ovgrpah_add_edge(ov_graph, lrs_to_align[i], 1, lrs_to_align[j], 0, aln1_rev, aln2_rev, cmp_aln1_invert, cmp_aln2_invert);
                    }
                }
                else if(front_type < 0)
                {
                    if(back_type >= 0)
                    {
                        // if(asm_is_lr_contained(lrs_to_align[j], target_id, mi, &mopt, tbuf, 1))
                        if(asm_is_alignment_ok(compact_lr_list, lrs_to_align[i], lrs_to_align[j], aln1_rev, aln2_rev, 1) < 10000)
                        {
                            ov_graph[lrs_to_align[i]].is_contained = 1;
                            ov_graph[lrs_to_align[j]].containing_list.push_back(lrs_to_align[i]);
                            fprintf(stdout, "#### %u contained ####\n", lrs_to_align[i]);
                        }
                    }
                    else // back_type < 0
                    {
                        if(asm_is_overlap_spurious(lr_list, lrs_to_align[j], target_id, mi, &mopt, tbuf, 1, -1)) continue;
                        // deque<cmp_aln_t> cmp_aln1, cmp_aln2;
                        deque<int32_t> cmp_aln1_invert, cmp_aln2_invert;
                        cmp_aln1_invert.resize(aln1_rev.size());
                        reverse_copy(aln1_rev.begin(), aln1_rev.end(), cmp_aln1_invert.begin());
                        cmp_aln2_invert.resize(aln2_rev.size());
                        reverse_copy(aln2_rev.begin(), aln2_rev.end(), cmp_aln2_invert.begin());
                        // asm_ovgraph_get_cmp_alignments(compact_lr_list, lrs_to_align[i], lrs_to_align[j], aln1_rev, aln2_rev, true, cmp_aln1, cmp_aln2, cmp_aln1_invert, cmp_aln2_invert);
                        asm_ovgrpah_add_edge(ov_graph, lrs_to_align[i], 0, lrs_to_align[j], 1, cmp_aln1_invert, cmp_aln2_invert, aln1_rev, aln2_rev);
                    }
                }
                else // front_type == 0
                {
                    if(back_type < 0)
                    {
                        // if(asm_is_lr_contained(lrs_to_align[j], target_id, mi, &mopt, tbuf, 0))
                        if(asm_is_alignment_ok(compact_lr_list, lrs_to_align[i], lrs_to_align[j], aln1_rev, aln2_rev, 1) < 10000)
                        {
                            ov_graph[lrs_to_align[j]].is_contained = 1;
                            ov_graph[lrs_to_align[i]].containing_list.push_back(lrs_to_align[j]);
                            fprintf(stdout, "#### %u contained ####\n", lrs_to_align[j]);
                        }
                    }
                    else if(back_type == 0)
                    {
                        if(lrs_to_align[i] < lrs_to_align[j])
                        {
                            // if(asm_is_lr_contained(lrs_to_align[j], target_id, mi, &mopt, tbuf, 0))
                            if(asm_is_alignment_ok(compact_lr_list, lrs_to_align[i], lrs_to_align[j], aln1_rev, aln2_rev, 1) < 10000)
                            {
                                ov_graph[lrs_to_align[j]].is_contained = 1;
                                ov_graph[lrs_to_align[i]].containing_list.push_back(lrs_to_align[j]);
                                fprintf(stdout, "#### %u contained ####\n", lrs_to_align[j]);
                            }
                        }
                        else
                        {
                            // if(asm_is_lr_contained(lrs_to_align[j], target_id, mi, &mopt, tbuf, 1))
                            if(asm_is_alignment_ok(compact_lr_list, lrs_to_align[i], lrs_to_align[j], aln1_rev, aln2_rev, 1) < 10000)
                            {
                                ov_graph[lrs_to_align[i]].is_contained = 1;
                                ov_graph[lrs_to_align[j]].containing_list.push_back(lrs_to_align[i]);
                                fprintf(stdout, "#### %u contained ####\n", lrs_to_align[i]);
                            }
                        }
                    }
                    else // back_type > 0
                    {
                        // if(asm_is_lr_contained(lrs_to_align[j], target_id, mi, &mopt, tbuf, 1))
                        if(asm_is_alignment_ok(compact_lr_list, lrs_to_align[i], lrs_to_align[j], aln1_rev, aln2_rev, 1) < 10000)
                        {
                            ov_graph[lrs_to_align[i]].is_contained = 1;
                            ov_graph[lrs_to_align[j]].containing_list.push_back(lrs_to_align[i]);
                            fprintf(stdout, "#### %u contained ####\n", lrs_to_align[i]);
                        }
                    }
                }
            }
        }
        mm_tbuf_destroy(tbuf);
        mm_idx_destroy(mi);
    }
    // print the overlap graph
    asm_ovgraph_transitive_reduction(ov_graph);
    asm_ovgraph_print_gfa(ov_graph, gopt.outDir + "/ov_graph.gfa");
    asm_ovgraph_print_reduced_gfa(ov_graph, gopt.outDir + "/ov_graph_all.gfa");
    asm_ovgraph_remove_transitive(ov_graph);
    asm_ovgraph_print_reduced_gfa(ov_graph, gopt.outDir + "/ov_graph_reduced.gfa");
    asm_ovgraph_get_paths(ov_graph, unused, compact_lr_list, ov_paths, covering_lr);
    // fprintf(stdout, "\n#### PATHS OF OVERLAP GRAPH ####\n");
    // for(i = 0; i < ov_paths.size(); i++)
    // {
    //     for(uint32_t z = 0; z < ov_paths[i].size(); z++)
    //     {
    //         fprintf(stdout, "%u:%c ", ov_paths[i][z].id, "+-"[ov_paths[i][z].strand]);
    //     }
    //     fprintf(stdout, "\n");
    // }
}
