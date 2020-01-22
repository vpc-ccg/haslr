// /*****************************************************
//  * Author: Ehsan Haghshenas (ehaghshe AT sfu DOT ca) *
//  *****************************************************/

#include "Align_LR2path.hpp"
#include <deque>
#include <tuple>

int  **_ap_lcs;
char **_ap_bt; // D diagonal, L left, U up

//   GOAL: store information about how long reads align onto the path
//  INPUT: 1) an extracted simple path 2) compacted long reads
// OUTPUT: data structure that for each consecutive pair of nodes in path stores which alignments should be used to extract the sub-sequence for every long read aligning to this pair of nodes
//   TODO: currently we only look at the pair of contigs to decide whether a long read supports both contigs or not. One improvement on this is to actually align compacted long read onto the simple path
void path_uniq_align_lrs(vector<vector<path_elem_t>> &simple_paths, vector<vector<Align_Seq_t*>> &compact_lr_list, vector<vector<id_strand_t>> &map_contig2lr)
{
    // initialize _asm_all_supp_lrs
    // vector<vector<vector<pair<Align_Seq_t*, Align_Seq_t*>>>> _asm_all_supp_lrs;
    // _asm_all_supp_lrs.resize(simple_paths.size());
    // 
    uint32_t i, j;
    for(i = 0; i < simple_paths.size(); i++)
    {
        vector<vector<pair<Align_Seq_t*, Align_Seq_t*> > > supp_lrs(simple_paths[i].size() - 1);
        for(j = 0; j < simple_paths[i].size() - 1; j++)
        {
            uint32_t contig1 = simple_paths[i][j].id;
            uint32_t strand1 = simple_paths[i][j].strand;
            uint32_t contig2 = simple_paths[i][j+1].id;
            uint32_t strand2 = simple_paths[i][j+1].strand;
            // fprintf(stdout, "[debug] %u:%c -> %u:%c\n", contig1, "+-"[strand1], contig2, "+-"[strand2]);
            for(uint32_t x = 0; x < map_contig2lr[contig1].size(); x++)
            {
                uint32_t lr_id = map_contig2lr[contig1][x].id;
                uint32_t lr_strand = map_contig2lr[contig1][x].strand;
                uint32_t count1 = 0;
                uint32_t count2 = 0;
                for(uint32_t y = 0; y < compact_lr_list[lr_id].size(); y++)
                {
                    if(compact_lr_list[lr_id][y]->t_id == contig1)
                    {
                        count1++;
                    }
                    if(compact_lr_list[lr_id][y]->t_id == contig2)
                    {
                        count2++;
                    }
                }
                if(count2 < 1) continue; // no alignment to the contig2
                if(count1 > 1 || count2 > 1) // multiple alignments to contig1 or contig2
                {
                    fprintf(stdout, "    duplicate alignment: contigs %u or %u on long read %u\n", contig1, contig2, lr_id);
                    continue;
                }
                // TODO: if there is only one alignment to contig1 and only one alignment to contig2, then one loop is enough (look below)
                int32_t comp1_id = -1;
                int32_t comp2_id = -1;
                if(lr_strand == strand1)
                {
                    for(int32_t y = compact_lr_list[lr_id].size() - 1; y >= 0; y--)
                    {
                        if(compact_lr_list[lr_id][y]->t_id == contig1)
                        {
                            comp1_id = y;
                            break;
                        }
                    }
                    for(int32_t y = comp1_id + 1; y < (int32_t)compact_lr_list[lr_id].size(); y++)
                    {
                        if(compact_lr_list[lr_id][y]->t_id == contig2)
                        {
                            comp2_id = y;
                            break;
                        }
                    }
                }
                else // lr_strand != strand1
                {
                    for(int32_t y = 0; y < (int32_t)compact_lr_list[lr_id].size(); y++)
                    {
                        if(compact_lr_list[lr_id][y]->t_id == contig1)
                        {
                            comp1_id = y;
                            break;
                        }
                    }
                    for(int32_t y = comp1_id - 1; y >= 0; y--)
                    {
                        if(compact_lr_list[lr_id][y]->t_id == contig2)
                        {
                            comp2_id = y;
                            break;
                        }
                    }
                }
                if(comp1_id != -1 && comp2_id != -1)
                {
                    if(compact_lr_list[lr_id][comp1_id]->is_rev == strand1 && compact_lr_list[lr_id][comp2_id]->is_rev == strand2
                        && compact_lr_list[lr_id][comp1_id]->q_end <= compact_lr_list[lr_id][comp2_id]->q_start)
                    {
                        // supp_lrs[j].push_back({compact_lr_list[lr_id][comp1_id], compact_lr_list[lr_id][comp2_id]});
                        // fprintf(stdout, "    +++ support %u\n", lr_id);
                        simple_paths[i][j].lrs.insert(lr_id);
                        simple_paths[i][j].alns[lr_id] = comp1_id;
                        simple_paths[i][j+1].lrs.insert(lr_id);
                        simple_paths[i][j+1].alns[lr_id] = comp2_id;
                    }
                    else if(compact_lr_list[lr_id][comp1_id]->is_rev != strand1 && compact_lr_list[lr_id][comp2_id]->is_rev != strand2
                        && compact_lr_list[lr_id][comp2_id]->q_end <= compact_lr_list[lr_id][comp1_id]->q_start)
                    {
                        // supp_lrs[j].push_back({compact_lr_list[lr_id][comp1_id], compact_lr_list[lr_id][comp2_id]});
                        // fprintf(stdout, "    +++ support %u\n", lr_id);
                        simple_paths[i][j].lrs.insert(lr_id);
                        simple_paths[i][j].alns[lr_id] = comp1_id;
                        simple_paths[i][j+1].lrs.insert(lr_id);
                        simple_paths[i][j+1].alns[lr_id] = comp2_id;
                    }
                }
            }
        }
        // _asm_all_supp_lrs[i] = supp_lrs;
        // fprintf(stdout, "=============\n");
    }
}

void path_repeat_align_single_lcs(vector<path_elem_t> &ov_path, vector<Align_Seq_t*> &compact_lr, deque<int32_t> &aln1, deque<int32_t> &aln2, int &score)
{
    int LCS_MATCH = 3;
    int LCS_INDEL = -1;
    int i, j;
    int m = compact_lr.size(); // number of rows
    int n = ov_path.size();    // number of columns
    // fprintf(stdout, "[path_repeat_align_single_lcs] ov_path: %zu\t compact_lr: %zu\n", ov_path.size(), compact_lr.size());
    // fflush(stdout);
    // int lcs[1000][1000];
    // char bt[1000][1000]; // D diagonal, L left, U up
    // init the first column
    for(i = 0; i <= m; i++)
    {
        // lcs[i][0] = 0;
        _ap_lcs[i][0] = i * LCS_INDEL;
        _ap_bt[i][0] = 'U'; // up
    }
    // init the first row, do not penalize beginning gaps
    for(j = 0; j <= n; j++)
    {
        _ap_lcs[0][j] = 0;
        _ap_bt[0][j] = 'L'; // left
    }
    // DP
    for(i = 1; i <= m; i++)
    {
        for(j = 1; j <= n; j++)
        {
            if(compact_lr[i-1]->t_id == ov_path[j-1].id && compact_lr[i-1]->is_rev == ov_path[j-1].strand)
            {
                _ap_lcs[i][j] = _ap_lcs[i-1][j-1] + LCS_MATCH;
                _ap_bt[i][j] = 'D'; // diagonal
            }
            else
            {
                if(_ap_lcs[i-1][j] > _ap_lcs[i][j-1])
                {
                    _ap_lcs[i][j] = _ap_lcs[i-1][j] + LCS_INDEL;
                    _ap_bt[i][j] = 'U'; // up
                }
                else
                {
                    _ap_lcs[i][j] = _ap_lcs[i][j-1] + LCS_INDEL;
                    _ap_bt[i][j] = 'L'; // lef
                }
            }
        }
    }
    // // do not penalize end gaps
    // for(i = 0; i < m; i++)
    // {
    //     if(lcs[i][n] > lcs[i+1][n])
    //     {
    //         lcs[i+1][n] = lcs[i][n];
    //         bt[i+1][n] = 'U';
    //     }
    // }
    // do not penalize end gaps (last row)
    for(j = 0; j < n; j++)
    {
        if(_ap_lcs[m][j] > _ap_lcs[m][j+1])
        {
            _ap_lcs[m][j+1] = _ap_lcs[m][j];
            _ap_bt[m][j+1] = 'L';
        }
    }
    // BT
    // deque<pair<int32_t, int8_t>> aln1, aln2;
    score = _ap_lcs[m][n];
    aln1.clear();
    aln2.clear();
    i = m;
    j = n;
    while(i > 0 || j > 0)
    {
        if(_ap_bt[i][j] == 'L') // left
        {
            aln1.push_front(-1);
            aln2.push_front(j-1);
            j--;
        }
        else if(_ap_bt[i][j] == 'U') // top
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
    //         fprintf(stdout, "%u\t", compact_lr[i-1]->t_id);
    //     for(j = 0; j <= n; j++)
    //     {
    //         fprintf(stdout, "%2d:%c\t", lcs[i][j], bt[i][j]);
    //     }
    //     fprintf(stdout, "\n");
    // }
}

void path_repeat_align_single_lcs_rev(vector<path_elem_t> &ov_path, vector<Align_Seq_t*> &compact_lr, deque<int32_t> &aln1, deque<int32_t> &aln2, int &score)
{
    int LCS_MATCH = 3;
    int LCS_INDEL = -1;
    int i, j;
    int m = compact_lr.size(); // number of rows
    int n = ov_path.size();   // number of columns
    // fprintf(stdout, "[path_repeat_align_single_lcs_rev] ov_path: %zu\t compact_lr: %zu\n", ov_path.size(), compact_lr.size());
    // fflush(stdout);
    // int lcs[1000][1000];
    // char bt[1000][1000]; // D diagonal, L left, U up
    // init the first column
    for(i = 0; i <= m; i++)
    {
        _ap_lcs[i][0] = i * LCS_INDEL;
        _ap_bt[i][0] = 'U'; // up
    }
    // init the first row, do not penalize beginning gaps
    for(j = 0; j <= n; j++)
    {
        _ap_lcs[0][j] = 0;
        _ap_bt[0][j] = 'L'; // left
    }
    // DP
    for(i = 1; i <= m; i++)
    {
        for(j = 1; j <= n; j++)
        {
            if(compact_lr[m-i]->t_id == ov_path[j-1].id && compact_lr[m-i]->is_rev != ov_path[j-1].strand)
            {
                _ap_lcs[i][j] = _ap_lcs[i-1][j-1] + LCS_MATCH;
                _ap_bt[i][j] = 'D'; // diagonal
            }
            else
            {
                if(_ap_lcs[i-1][j] > _ap_lcs[i][j-1])
                {
                    _ap_lcs[i][j] = _ap_lcs[i-1][j] + LCS_INDEL;
                    _ap_bt[i][j] = 'U'; // up
                }
                else
                {
                    _ap_lcs[i][j] = _ap_lcs[i][j-1] + LCS_INDEL;
                    _ap_bt[i][j] = 'L'; // lef
                }
            }
        }
    }
    // // do not penalize end gaps
    // for(i = 0; i < m; i++)
    // {
    //     if(_ap_lcs[i][n] > _ap_lcs[i+1][n])
    //     {
    //         _ap_lcs[i+1][n] = _ap_lcs[i][n];
    //         _ap_bt[i+1][n] = 'U';
    //     }
    // }
    // do not penalize end gaps
    for(j = 0; j < n; j++)
    {
        if(_ap_lcs[m][j] > _ap_lcs[m][j+1])
        {
            _ap_lcs[m][j+1] = _ap_lcs[m][j];
            _ap_bt[m][j+1] = 'L';
        }
    }
    // BT
    // deque<pair<int32_t, int8_t>> aln1, aln2;
    score = _ap_lcs[m][n];
    aln1.clear();
    aln2.clear();
    i = m;
    j = n;
    while(i > 0 || j > 0)
    {
        if(_ap_bt[i][j] == 'L') // left
        {
            aln1.push_front(-1);
            aln2.push_front(j-1);
            j--;
        }
        else if(_ap_bt[i][j] == 'U') // top
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

void path_repeat_align_all_lrs(vector<vector<path_elem_t>> &ov_paths, vector<vector<uint32_t>> &covering_lr, vector<vector<Align_Seq_t*>> &compact_lr_list)
{
    _ap_lcs = (int**)malloc(1000 * sizeof(int*));
    for(size_t i = 0; i < 1000; i++) _ap_lcs[i] = (int*)malloc(10000 * sizeof(int));
    _ap_bt = (char**)malloc(1000 * sizeof(char*));
    for(size_t i = 0; i < 1000; i++) _ap_bt[i] = (char*)malloc(10000 * sizeof(char));
    // TODO: realloc if needed
    fprintf(stdout, "\n#### PATHS OF OVERLAP GRAPH ####\n");
    for(size_t i = 0; i < ov_paths.size(); i++)
    {
        fprintf(stdout, "\n$$$ ");
        for(uint32_t z = 0; z < ov_paths[i].size(); z++)
        {
            fprintf(stdout, "%u:%c ", ov_paths[i][z].id, "+-"[ov_paths[i][z].strand]);
        }
        fprintf(stdout, "\n+++ ");
        for(uint32_t z = 0; z < covering_lr[i].size(); z++)
        {
            fprintf(stdout, "%u\t", covering_lr[i][z]);
        }
        fprintf(stdout, "\n");
        fflush(stdout);
        // 
        for(size_t j = 0; j < covering_lr[i].size(); j++)
        {
            uint32_t lr_id = covering_lr[i][j];
            // fprintf(stdout, "\naligning read: %u\n", lr_id);
            // fflush(stdout);
            deque<int32_t> aln1, aln2;
            deque<int32_t> aln1_rev, aln2_rev;
            int score, score_rev;
            path_repeat_align_single_lcs(ov_paths[i], compact_lr_list[lr_id], aln1, aln2, score);
            path_repeat_align_single_lcs_rev(ov_paths[i], compact_lr_list[lr_id], aln1_rev, aln2_rev, score_rev);
            if(score > score_rev)
            {
                // // 
                // for(uint32_t z = 0; z < aln1.size(); z++)
                //     if(aln1[z] == -1)
                //         fprintf(stdout, "%d\t", aln1[z]);
                //     else
                //         fprintf(stdout, "%u:%u\t", compact_lr_list[lr_id][aln1[z]]->t_id, compact_lr_list[lr_id][aln1[z]]->is_rev);
                // fprintf(stdout, "\n");
                // for(uint32_t z = 0; z < aln2.size(); z++)
                //     if(aln2[z] == -1)
                //         fprintf(stdout, "%d\t", aln2[z]);
                //     else
                //         fprintf(stdout, "%u:%u\t", ov_paths[i][aln2[z]].id, ov_paths[i][aln2[z]].strand);
                // fprintf(stdout, "\n");
                for(uint32_t z = 0; z < aln1.size(); z++)
                {
                    if(aln1[z] != -1 && aln2[z] != -1) // match
                    {
                        ov_paths[i][aln2[z]].lrs.insert(lr_id);
                        ov_paths[i][aln2[z]].alns[lr_id] = aln1[z];
                    }
                }
            }
            else
            {
                // // 
                // for(uint32_t z = 0; z < aln1_rev.size(); z++)
                //     if(aln1_rev[z] == -1)
                //         fprintf(stdout, "%d\t", aln1_rev[z]);
                //     else
                //         fprintf(stdout, "%u:%u\t", compact_lr_list[lr_id][aln1_rev[z]]->t_id, 1 - compact_lr_list[lr_id][aln1_rev[z]]->is_rev);
                // fprintf(stdout, "\n");
                // for(uint32_t z = 0; z < aln2_rev.size(); z++)
                //     if(aln2_rev[z] == -1)
                //         fprintf(stdout, "%d\t", aln2_rev[z]);
                //     else
                //         fprintf(stdout, "%u:%u\t", ov_paths[i][aln2_rev[z]].id, ov_paths[i][aln2_rev[z]].strand);
                // fprintf(stdout, "\n");
                for(uint32_t z = 0; z < aln1_rev.size(); z++)
                {
                    if(aln1_rev[z] != -1 && aln2_rev[z] != -1) // match
                    {
                        ov_paths[i][aln2_rev[z]].lrs.insert(lr_id);
                        ov_paths[i][aln2_rev[z]].alns[lr_id] = aln1_rev[z];
                    }
                }
            }
        }
        // for(uint32_t j = 0; j < ov_paths[i].size(); j++)
        // {
        //     fprintf(stdout, "=== %u:%c \n", ov_paths[i][j].id, "+-"[ov_paths[i][j].strand]);
        //     for(set<uint32_t>::iterator it = ov_paths[i][j].lrs.begin(); it != ov_paths[i][j].lrs.end(); it++)
        //     {
        //         fprintf(stdout, "%u\t", *it);
        //     }
        //     fprintf(stdout, "\n");
        // }
        // break;
    }
}

int32_t bridge_find_next_node_forward(vector<path_elem_t> &ov_path, int32_t curr_index, int32_t last_index)
{
    uint32_t best_size = 0;
    int32_t best_index = -1;
    set<uint32_t> &curr_lrs = ov_path[curr_index].lrs;
    int32_t next_index = curr_index + 1;
    while(next_index <= last_index)
    {
        // fprintf(stdout, "[debug] curr_index: %u \tnext_index: %u\n", curr_index, next_index);
        // fflush(stdout);
        vector<uint32_t> shared_lrs;
        std::set_intersection(curr_lrs.begin(), curr_lrs.end(), ov_path[next_index].lrs.begin(), ov_path[next_index].lrs.end(), std::inserter(shared_lrs, shared_lrs.begin()));
        // fprintf(stdout, "[debug] intersection: %zu\n", shared_lrs.size());
        // fflush(stdout);
        if(shared_lrs.size() == 0) // there is no shared long read between curr_index and nodes starting next_index and onward
        {
            break;
        }
        else if(shared_lrs.size() > best_size)
        {
            best_index = next_index;
            best_size = shared_lrs.size();
        }
        next_index++;
    }
    return best_index;
}

int32_t bridge_find_next_node_backward(vector<path_elem_t> &ov_path, int32_t curr_index, int32_t first_index)
{
    uint32_t best_size = 0;
    int32_t best_index = -1;
    set<uint32_t> &curr_lrs = ov_path[curr_index].lrs;
    int32_t next_index = curr_index - 1;
    while(next_index >= first_index)
    {
        // fprintf(stdout, "[debug] curr_index: %u \tnext_index: %u\n", curr_index, next_index);
        // fflush(stdout);
        vector<uint32_t> shared_lrs;
        std::set_intersection(curr_lrs.begin(), curr_lrs.end(), ov_path[next_index].lrs.begin(), ov_path[next_index].lrs.end(), std::inserter(shared_lrs, shared_lrs.begin()));
        // fprintf(stdout, "[debug] intersection: %zu\n", shared_lrs.size());
        // fflush(stdout);
        if(shared_lrs.size() == 0) // there is no shared long read between curr_index and nodes starting next_index and onward
        {
            break;
        }
        else if(shared_lrs.size() > best_size)
        {
            best_index = next_index;
            best_size = shared_lrs.size();
        }
        next_index--;
    }
    return best_index;
}

// TODO: generalize the idea of connecting simple paths; relax the condition of nodes that can be connected (right now only tails can be connected)
void bridge_between_simple_paths(vector<vector<path_elem_t>> &simple_paths, vector<vector<path_elem_t>> &ov_paths, vector<vector<Align_Seq_t*>> &compact_lr_list)
{
    // tail_node_id,   <path_id, index in path>
    map<uint32_t, pair<uint32_t, uint32_t>> tails;
    for(uint32_t j = 0; j < simple_paths.size(); j++)
    {
        tails[simple_paths[j].front().id] = {j, 0}; // head
        tails[simple_paths[j].back().id]  = {j, simple_paths[j].size() - 1}; // tail
    }
    vector<uint8_t> paths_to_delete(simple_paths.size(), 0);
    // for(uint32_t i = 8; i <= 8; i++)
    for(uint32_t i = 0; i < ov_paths.size(); i++)
    {
        //          <ov_path_index, simple_path_id, simple_path_index>
        vector<tuple<uint32_t,      uint32_t,       uint32_t>> candid; // paths that share the tail
        for(uint32_t j = 0; j < ov_paths[i].size(); j++)
        {
            if(tails.count(ov_paths[i][j].id) > 0) // found in the map
            {
                candid.push_back(make_tuple(j, tails[ov_paths[i][j].id].first, tails[ov_paths[i][j].id].second));
            }
        }
        fprintf(stdout, "\n$$$ ");
        for(uint32_t z = 0; z < ov_paths[i].size(); z++)
        {
            fprintf(stdout, "%u:%c ", ov_paths[i][z].id, "+-"[ov_paths[i][z].strand]);
        }
        fprintf(stdout, "\ncandidate simple paths\n");
        for(uint32_t j = 0; j < candid.size(); j++)
        {
            uint32_t ov_path_index, simple_path_id, simple_path_index;
            tie(ov_path_index, simple_path_id, simple_path_index) = candid[j];
            for(uint32_t z = 0; z < simple_paths[simple_path_id].size(); z++)
            {
                fprintf(stdout, "%u:%c\t", simple_paths[simple_path_id][z].id, "+-"[simple_paths[simple_path_id][z].strand]);
            }
            fprintf(stdout, "\n");
            fflush(stdout);
        }
        if(candid.size() == 2)
        {
            uint32_t ov_path_index1, simple_path_id1, simple_path_index1;
            tie(ov_path_index1, simple_path_id1, simple_path_index1) = candid[0];
            uint32_t ov_path_index2, simple_path_id2, simple_path_index2;
            tie(ov_path_index2, simple_path_id2, simple_path_index2) = candid[1];
            fprintf(stdout, "path1 - id: %u - size: %zu - index: %u\n", simple_path_id1, simple_paths[simple_path_id1].size(), simple_path_index1);
            fprintf(stdout, "path2 - id: %u - size: %zu - index: %u\n", simple_path_id2, simple_paths[simple_path_id2].size(), simple_path_index2);
            // 
            if(simple_path_id1 != simple_path_id2) // connecting two different paths
            {
                fprintf(stdout, "CASE 1: connecting two different paths!\n");
                fflush(stdout);
                // 1) find the heaviset sequence of SRCs from the ov_path
                vector<path_elem_t> heaviest_path;
                uint32_t curr_index = ov_path_index1;
                heaviest_path.push_back(ov_paths[i][curr_index]);
                while(curr_index != ov_path_index2)
                {
                    int32_t ret_index = bridge_find_next_node_forward(ov_paths[i], curr_index, ov_path_index2);
                    // fprintf(stdout, "next_index: %d\n", ret_index);
                    // fflush(stdout);
                    if(ret_index != -1)
                    {
                        curr_index = ret_index;
                        // fprintf(stdout, "next %u:%c\n", ov_paths[i][curr_index].id, "+-"[ov_paths[i][curr_index].strand]);
                        // fflush(stdout);
                        heaviest_path.push_back(ov_paths[i][curr_index]);
                    }
                    else
                    {
                        fprintf(stdout, "UNEXPECTED!\n");
                        fflush(stdout);
                        exit(0);
                    }
                }
                fprintf(stdout, "heaviest path\n");
                for(uint32_t z = 0; z < heaviest_path.size(); z++)
                {
                    fprintf(stdout, "%u:%c\t", heaviest_path[z].id, "+-"[heaviest_path[z].strand]);
                }
                fprintf(stdout, "\n");
                fflush(stdout);
                // check correct direction
                // TODO: this automatically won't connect any simple_path of size 1
                bool is_ok = false;
                if(heaviest_path.front().strand == simple_paths[simple_path_id1][simple_path_index1].strand && simple_path_index1 > simple_paths[simple_path_id1].size() / 2)
                    is_ok = true;
                if(heaviest_path.front().strand != simple_paths[simple_path_id1][simple_path_index1].strand && simple_path_index1 < simple_paths[simple_path_id1].size() / 2)
                    is_ok = true;
                if(heaviest_path.back().strand == simple_paths[simple_path_id2][simple_path_index2].strand && simple_path_index2 < simple_paths[simple_path_id2].size() / 2)
                    is_ok = true;
                if(heaviest_path.back().strand != simple_paths[simple_path_id2][simple_path_index2].strand && simple_path_index2 > simple_paths[simple_path_id2].size() / 2)
                    is_ok = true;
                if(is_ok == false) continue;
                // 2) connect paths using heaviest path
                vector<path_elem_t> connected_path;
                // step 1: add the first simple path
                if(simple_paths[simple_path_id1][simple_path_index1].strand == heaviest_path.front().strand)
                {
                    fprintf(stdout, "here1\n");
                    fflush(stdout);
                    for(uint32_t j = 0; j < simple_path_index1; j++)
                    {
                        connected_path.push_back(simple_paths[simple_path_id1][j]);
                    }
                    fprintf(stdout, "here11\n");
                    fflush(stdout);
                }
                else
                {
                    fprintf(stdout, "here2\n");
                    fflush(stdout);
                    for(int32_t j = simple_paths[simple_path_id1].size() - 1; j > (int32_t)simple_path_index1; j--)
                    {
                        connected_path.push_back({(uint32_t)1 - simple_paths[simple_path_id1][j].strand, simple_paths[simple_path_id1][j].id, simple_paths[simple_path_id1][j].lrs, simple_paths[simple_path_id1][j].alns});
                    }
                    fprintf(stdout, "here22\n");
                    fflush(stdout);
                }
                // step 2: add the heaviest path
                // TODO: note that insert function of c++ map does not update the value if key exists. This means we need to be sure that the alignments chosen for simple path and overlap path are the same!
                // -- update set of long reads and alignments on the tails of heaviest_path
                heaviest_path.front().lrs.insert(simple_paths[simple_path_id1][simple_path_index1].lrs.begin(), simple_paths[simple_path_id1][simple_path_index1].lrs.end());
                heaviest_path.front().alns.insert(simple_paths[simple_path_id1][simple_path_index1].alns.begin(), simple_paths[simple_path_id1][simple_path_index1].alns.end());
                heaviest_path.back().lrs.insert(simple_paths[simple_path_id2][simple_path_index2].lrs.begin(), simple_paths[simple_path_id2][simple_path_index2].lrs.end());
                heaviest_path.back().alns.insert(simple_paths[simple_path_id2][simple_path_index2].alns.begin(), simple_paths[simple_path_id2][simple_path_index2].alns.end());
                fprintf(stdout, "here3\n");
                fflush(stdout);
                // -- add path now
                for(uint32_t j = 0; j < heaviest_path.size(); j++)
                    connected_path.push_back(heaviest_path[j]);
                // step 3: add the second simple path
                if(simple_paths[simple_path_id2][simple_path_index2].strand == heaviest_path.back().strand)
                {
                    fprintf(stdout, "here4\n");
                    fflush(stdout);
                    for(uint32_t j = simple_path_index2 + 1; j < simple_paths[simple_path_id2].size(); j++)
                    {
                        connected_path.push_back(simple_paths[simple_path_id2][j]);
                    }
                    fprintf(stdout, "here44\n");
                    fflush(stdout);
                }
                else
                {
                    fprintf(stdout, "here5\n");
                    fflush(stdout);
                    for(int32_t j = simple_path_index2 - 1; j >= 0; j--)
                    {
                        connected_path.push_back({(uint32_t)1 - simple_paths[simple_path_id2][j].strand, simple_paths[simple_path_id2][j].id, simple_paths[simple_path_id2][j].lrs, simple_paths[simple_path_id2][j].alns});
                    }
                    fprintf(stdout, "here55\n");
                    fflush(stdout);
                }
                // 3) add the connected path and mark simple paths for deletion
                // update tails
                tails.erase(simple_paths[simple_path_id1][simple_path_index1].id);
                tails.erase(simple_paths[simple_path_id2][simple_path_index2].id);
                // update path1
                simple_paths[simple_path_id1] = connected_path;
                if(tails.count(simple_paths[simple_path_id1].front().id) > 0)
                    tails[simple_paths[simple_path_id1].front().id] = {simple_path_id1, 0}; // head
                if(tails.count(simple_paths[simple_path_id1].back().id) > 0)
                tails[simple_paths[simple_path_id1].back().id]  = {simple_path_id1, simple_paths[simple_path_id1].size() - 1}; // tail
                // delete path2
                paths_to_delete[simple_path_id2] = 1;
                // 
                fprintf(stdout, "connected path\n");
                fflush(stdout);
                for(uint32_t z = 0; z < connected_path.size(); z++)
                {
                    fprintf(stdout, "%u:%c\t", connected_path[z].id, "+-"[connected_path[z].strand]);
                }
                fprintf(stdout, "\n");
                fflush(stdout);
            }
            else // is it a palindrome chromosome?
            {
                // TODO: deal with palindromic chromosomes
                fprintf(stdout, "CASE 2: Palindromic chromosome? not implemented yet!\n");
                fflush(stdout);
            }
        }
        else if(candid.size() == 1) // extend a simple path to the left/right
        {
            uint32_t ov_path_index, simple_path_id, simple_path_index;
            tie(ov_path_index, simple_path_id, simple_path_index) = candid[0];
            fprintf(stdout, "path1 - id: %u - size: %zu - index: %u\n", simple_path_id, simple_paths[simple_path_id].size(), simple_path_index);
            // TODO: what is a good way to figure out whether we should extend to the left or right?
            if(simple_path_index < simple_paths[simple_path_id].size() / 2) // extend to the left
            {
                if(ov_paths[i][ov_path_index].strand == simple_paths[simple_path_id][simple_path_index].strand) // iterate backward on overlap path
                {
                    fprintf(stdout, "CASE 3: extending to the left, iterating backward\n");
                    fflush(stdout);
                    vector<path_elem_t> connected_path;
                    // 1) add the heaviest path
                    uint32_t curr_index = 0;
                    uint32_t ov_path_last = ov_path_index;
                    connected_path.push_back(ov_paths[i][curr_index]);
                    while(curr_index < ov_path_last)
                    {
                        int32_t ret_index = bridge_find_next_node_forward(ov_paths[i], curr_index, ov_path_last);
                        fprintf(stdout, "next_index: %d\n", ret_index);
                        fflush(stdout);
                        if(ret_index != -1)
                        {
                            curr_index = ret_index;
                            fprintf(stdout, "next %u:%c\n", ov_paths[i][curr_index].id, "+-"[ov_paths[i][curr_index].strand]);
                            fflush(stdout);
                            connected_path.push_back(ov_paths[i][curr_index]);
                        }
                        else
                        {
                            break;
                        }
                    }
                    // 2) add the simple path
                    connected_path.back().lrs.insert(simple_paths[simple_path_id][simple_path_index].lrs.begin(), simple_paths[simple_path_id][simple_path_index].lrs.end());
                    connected_path.back().alns.insert(simple_paths[simple_path_id][simple_path_index].alns.begin(), simple_paths[simple_path_id][simple_path_index].alns.end());
                    for(uint32_t j = simple_path_index + 1; j < simple_paths[simple_path_id].size(); j++)
                    {
                        connected_path.push_back(simple_paths[simple_path_id][j]);
                    }
                    fprintf(stdout, "connected path\n");
                    for(uint32_t z = 0; z < connected_path.size(); z++)
                    {
                        fprintf(stdout, "%u:%c\t", connected_path[z].id, "+-"[connected_path[z].strand]);
                    }
                    fprintf(stdout, "\n");
                    fflush(stdout);
                    // // update tails
                    tails.erase(simple_paths[simple_path_id][simple_path_index].id);
                    simple_paths[simple_path_id] = connected_path;
                    if(tails.count(simple_paths[simple_path_id].back().id) > 0)
                        tails[simple_paths[simple_path_id].back().id]  = {simple_path_id, simple_paths[simple_path_id].size() - 1}; // tail
                }
                else
                {
                    fprintf(stdout, "CASE 4: extending to the left, iterating forward\n");
                    fflush(stdout);
                    vector<path_elem_t> connected_path;
                    // 1) add the heaviest path
                    uint32_t curr_index = ov_paths[i].size() - 1;
                    uint32_t ov_path_first = ov_path_index;
                    connected_path.push_back({(uint32_t)1 - ov_paths[i][curr_index].strand, ov_paths[i][curr_index].id, ov_paths[i][curr_index].lrs, ov_paths[i][curr_index].alns});
                    while(curr_index > ov_path_first)
                    {
                        int32_t ret_index = bridge_find_next_node_backward(ov_paths[i], curr_index, ov_path_first);
                        fprintf(stdout, "next_index: %d\n", ret_index);
                        fflush(stdout);
                        if(ret_index != -1)
                        {
                            curr_index = ret_index;
                            fprintf(stdout, "next %u:%c\n", ov_paths[i][curr_index].id, "+-"[ov_paths[i][curr_index].strand]);
                            fflush(stdout);
                            connected_path.push_back({(uint32_t)1 - ov_paths[i][curr_index].strand, ov_paths[i][curr_index].id, ov_paths[i][curr_index].lrs, ov_paths[i][curr_index].alns});
                        }
                        else
                        {
                            break;
                        }
                    }
                    // 2) add the simple path
                    connected_path.back().lrs.insert(simple_paths[simple_path_id][simple_path_index].lrs.begin(), simple_paths[simple_path_id][simple_path_index].lrs.end());
                    connected_path.back().alns.insert(simple_paths[simple_path_id][simple_path_index].alns.begin(), simple_paths[simple_path_id][simple_path_index].alns.end());
                    for(uint32_t j = simple_path_index + 1; j < simple_paths[simple_path_id].size(); j++)
                    {
                        connected_path.push_back(simple_paths[simple_path_id][j]);
                    }
                    fprintf(stdout, "connected path\n");
                    for(uint32_t z = 0; z < connected_path.size(); z++)
                    {
                        fprintf(stdout, "%u:%c\t", connected_path[z].id, "+-"[connected_path[z].strand]);
                    }
                    fprintf(stdout, "\n");
                    fflush(stdout);
                    // // update tails
                    tails.erase(simple_paths[simple_path_id][simple_path_index].id);
                    simple_paths[simple_path_id] = connected_path;
                    if(tails.count(simple_paths[simple_path_id].back().id) > 0)
                        tails[simple_paths[simple_path_id].back().id]  = {simple_path_id, simple_paths[simple_path_id].size() - 1}; // tail
                }
            }
            else // extend to the right
            {
                if(ov_paths[i][ov_path_index].strand == simple_paths[simple_path_id][simple_path_index].strand) // iterate forward on overlap path
                {
                    fprintf(stdout, "CASE 5: extending to the right, iterating forward\n");
                    fflush(stdout);
                    vector<path_elem_t> connected_path;
                    // 1) add the simple path
                    for(uint32_t j = 0; j <= simple_path_index; j++)
                    {
                        connected_path.push_back(simple_paths[simple_path_id][j]);
                    }
                    connected_path.back().lrs.insert(ov_paths[i][ov_path_index].lrs.begin(), ov_paths[i][ov_path_index].lrs.end());
                    connected_path.back().alns.insert(ov_paths[i][ov_path_index].alns.begin(), ov_paths[i][ov_path_index].alns.end());
                    // 2) add the heaviest path
                    uint32_t curr_index = ov_path_index;
                    uint32_t ov_path_last = ov_paths[i].size() - 1;
                    while(curr_index < ov_path_last)
                    {
                        int32_t ret_index = bridge_find_next_node_forward(ov_paths[i], curr_index, ov_path_last);
                        fprintf(stdout, "next_index: %d\n", ret_index);
                        fflush(stdout);
                        if(ret_index != -1)
                        {
                            curr_index = ret_index;
                            fprintf(stdout, "next %u:%c\n", ov_paths[i][curr_index].id, "+-"[ov_paths[i][curr_index].strand]);
                            fflush(stdout);
                            connected_path.push_back(ov_paths[i][curr_index]);
                        }
                        else
                        {
                            break;
                        }
                    }
                    fprintf(stdout, "connected path\n");
                    for(uint32_t z = 0; z < connected_path.size(); z++)
                    {
                        fprintf(stdout, "%u:%c\t", connected_path[z].id, "+-"[connected_path[z].strand]);
                    }
                    fprintf(stdout, "\n");
                    fflush(stdout);
                    // update tails
                    tails.erase(simple_paths[simple_path_id][simple_path_index].id);
                    simple_paths[simple_path_id] = connected_path;
                    if(tails.count(simple_paths[simple_path_id].front().id) > 0)
                        tails[simple_paths[simple_path_id].front().id] = {simple_path_id, 0}; // head
                }
                else // iterate backward on overlap path
                {
                    fprintf(stdout, "CASE 6: extending to the right, iterating backward\n");
                    fflush(stdout);
                    vector<path_elem_t> connected_path;
                    // 1) add the simple path
                    for(uint32_t j = 0; j <= simple_path_index; j++)
                    {
                        connected_path.push_back(simple_paths[simple_path_id][j]);
                    }
                    connected_path.back().lrs.insert(ov_paths[i][ov_path_index].lrs.begin(), ov_paths[i][ov_path_index].lrs.end());
                    connected_path.back().alns.insert(ov_paths[i][ov_path_index].alns.begin(), ov_paths[i][ov_path_index].alns.end());
                    // 2) add the heaviest path
                    uint32_t curr_index = ov_path_index;
                    uint32_t ov_path_first = 0;
                    while(curr_index > ov_path_first)
                    {
                        int32_t ret_index = bridge_find_next_node_backward(ov_paths[i], curr_index, ov_path_first);
                        fprintf(stdout, "next_index: %d\n", ret_index);
                        fflush(stdout);
                        if(ret_index != -1)
                        {
                            curr_index = ret_index;
                            fprintf(stdout, "next %u:%c\n", ov_paths[i][curr_index].id, "+-"[ov_paths[i][curr_index].strand]);
                            fflush(stdout);
                            connected_path.push_back({(uint32_t)1 - ov_paths[i][curr_index].strand, ov_paths[i][curr_index].id, ov_paths[i][curr_index].lrs, ov_paths[i][curr_index].alns});
                        }
                        else
                        {
                            break;
                        }
                    }
                    fprintf(stdout, "connected path\n");
                    for(uint32_t z = 0; z < connected_path.size(); z++)
                    {
                        fprintf(stdout, "%u:%c\t", connected_path[z].id, "+-"[connected_path[z].strand]);
                    }
                    fprintf(stdout, "\n");
                    fflush(stdout);
                    // update tails
                    tails.erase(simple_paths[simple_path_id][simple_path_index].id);
                    simple_paths[simple_path_id] = connected_path;
                    if(tails.count(simple_paths[simple_path_id].front().id) > 0)
                        tails[simple_paths[simple_path_id].front().id] = {simple_path_id, 0}; // head
                }
            }
        }
        // else // TODO: check cases that path.size() > 2
        // {
        // }
    }
    // 
    vector<vector<path_elem_t> > simple_paths_tmp;
    for(uint32_t j = 0; j < simple_paths.size(); j++)
        if(paths_to_delete[j] == 0)
            simple_paths_tmp.push_back(simple_paths[j]);
    simple_paths = simple_paths_tmp;
}