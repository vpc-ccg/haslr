// /*****************************************************
//  * Author: Ehsan Haghshenas (ehaghshe AT sfu DOT ca) *
//  *****************************************************/

#include "Backbone_graph.hpp"
#include <algorithm>
#include <unordered_map>
#include <queue>

void bbg_add_edge(vector<BBG_Node_t> &graph, uint32_t lr_id, uint32_t lr_strand, vector<Align_Seq_t*> &compact_lr, uint32_t index1, uint32_t index2)
{
    uint32_t node1 = compact_lr[index1]->t_id;
    uint32_t rev1  = compact_lr[index1]->is_rev;
    uint32_t node2 = compact_lr[index2]->t_id;
    uint32_t rev2  = compact_lr[index2]->is_rev;
    // TODO: which one should I use?
    // if(node1 == node2 && rev1 != rev2) return; // most probably palindromic
    // if(node1 == node2) return; // most probably palindromic
    // fprintf(stderr, "[debug] %u:%c -> %u:%c on lr:%u\n", node1, (rev1 ? '-' : '+'), node2, (rev2 ? '-' : '+'), lr_id);
    // 
    uint32_t to_node1 = (node2 << 1) | rev2;
    uint32_t to_node2 = (node1 << 1) | (1 - rev1);
    graph[node1].edges[rev1][to_node1].edge_supp.push_back({lr_id, lr_strand, index1, index2});
    graph[node2].edges[1-rev2][to_node2].edge_supp.push_back({lr_id, 1 - lr_strand, index2, index1});
}

void bbg_add_edge_with_supp(vector<BBG_Node_t> &graph, uint32_t node1, uint8_t rev1, uint32_t node2, uint8_t rev2, vector<Edge_Supp_t> &shared_supp)
{
    // 
    uint32_t to_node1 = (node2 << 1) | rev2;
    uint32_t to_node2 = (node1 << 1) | (1 - rev1);
    for(uint32_t i = 0; i < shared_supp.size(); i++)
    {
        graph[node1].edges[rev1][to_node1].edge_supp.push_back({shared_supp[i].lr_id, shared_supp[i].lr_strand, shared_supp[i].cmp_head_id, shared_supp[i].cmp_tail_id});
        graph[node2].edges[1-rev2][to_node2].edge_supp.push_back({shared_supp[i].lr_id, uint32_t(1 - shared_supp[i].lr_strand), shared_supp[i].cmp_tail_id, shared_supp[i].cmp_head_id});
    }
}

BBG_Edge_t* bbg_get_edge(uint32_t node1, uint8_t rev1, uint32_t node2, uint8_t rev2, vector<BBG_Node_t> &graph)
{
    uint32_t to_node1 = (node2 << 1) | rev2;
    return &(graph[node1].edges[rev1][to_node1]);
}

void bbg_remove_edge(uint32_t node1, uint8_t rev1, uint32_t node2, uint8_t rev2, vector<BBG_Node_t> &graph)
{
    uint32_t to_node1 = (node2 << 1) | rev2;
    uint32_t to_node2 = (node1 << 1) | (1 - rev1);
    graph[node1].edges[rev1].erase(to_node1);
    graph[node2].edges[1-rev2].erase(to_node2);
}

void bbg_add_repetitive_SRCs(vector<BBG_Node_t> &graph, Contig_List_t &contig_list, vector<vector<Align_Seq_t*>> &compact_lr_list)
{
    vector<BBG_Node_t> tmp_bbg;
    tmp_bbg.resize(contig_list.contigs_size);
    for(uint32_t i = 0; i < compact_lr_list.size(); i++)
    {
        if(compact_lr_list[i].size() > 1)
        {
            // add one edge for each pair of consecutive contigs on the long read
            for(uint32_t j = 0; j < compact_lr_list[i].size() - 1; j++)
            {
                bbg_add_edge(tmp_bbg, i, 0, compact_lr_list[i], j, j + 1);
            }
        }
    }
    // 
    bbg_remove_weak_edges(tmp_bbg);
    bbg_print_graph_gfa(tmp_bbg, contig_list, gopt.out_dir + "/backbone.tmp.gfa");
    // 
    // int cnt = 0;
    for(uint32_t i = 0; i < graph.size(); i++)
    {
        for(int rev = 0; rev < 2; rev++)
        {
            for(map<uint32_t, BBG_Edge_t>::iterator it = graph[i].edges[rev].begin(); it != graph[i].edges[rev].end(); it++)
            {
                uint32_t contig1 = i;
                uint32_t strand1 = rev;
                uint32_t contig2 = (it->first >> 1);
                uint32_t strand2 = (it->first & 1);
                uint32_t to_node1 = (contig2 << 1) | strand2;
                uint32_t to_node2 = (contig1 << 1) | (1 - strand1);
                if(tmp_bbg[contig1].edges[strand1].count(to_node1) == 0) // if edge doesn't exist in BBG2
                {
                    // update
                    BBG_Edge_t &edge1 = graph[contig1].edges[strand1][to_node1];
                    fprintf(stdout, "check %u:%c -> %u:%c\tsupp:%zu\n", contig1, "+-"[strand1], contig2, "+-"[strand2], edge1.edge_supp.size());
                    for(uint32_t j = 0; j < edge1.edge_supp.size(); j++)
                    {
                        uint32_t lr_id     = edge1.edge_supp[j].lr_id;
                        uint32_t lr_strand = edge1.edge_supp[j].lr_strand;
                        uint32_t cmp_beg   = edge1.edge_supp[j].cmp_head_id;
                        uint32_t cmp_end   = edge1.edge_supp[j].cmp_tail_id;
                        fprintf(stdout, "debug lr:%u strand:%u cmp_beg:%u cmp_end:%u\n", lr_id, lr_strand, cmp_beg, cmp_end);
                        if(lr_strand == 0) // iterate forward
                        {
                            // fprintf(stdout, "here1\n");
                            for(uint32_t z = cmp_beg; z < cmp_end; z++)
                            {
                                // fprintf(stdout, "here2\n");
                                uint32_t node1 = compact_lr_list[lr_id][z]->t_id;
                                uint32_t rev1 = compact_lr_list[lr_id][z]->is_rev;
                                uint32_t node2 = compact_lr_list[lr_id][z+1]->t_id;
                                uint32_t rev2 = compact_lr_list[lr_id][z+1]->is_rev;
                                uint32_t tonode1 = (node2 << 1) | rev2;
                                if(tmp_bbg[node1].edges[rev1].count(tonode1) > 0)
                                {
                                    fprintf(stdout, "\tremove %u:%c -> %u:%c\tsupp:%zu\n", node1, "+-"[rev1], node2, "+-"[rev2], tmp_bbg[node1].edges[rev1][tonode1].edge_supp.size());
                                    bbg_remove_edge(node1, rev1, node2, rev2, tmp_bbg);
                                }
                            }
                        }
                        else  // iterate backward
                        {
                            // fprintf(stdout, "here3\n");
                            for(uint32_t z = cmp_beg; z > cmp_end; z--)
                            {
                                // fprintf(stdout, "here4\n");
                                uint32_t node1 = compact_lr_list[lr_id][z]->t_id;
                                uint32_t rev1 = compact_lr_list[lr_id][z]->is_rev;
                                uint32_t node2 = compact_lr_list[lr_id][z-1]->t_id;
                                uint32_t rev2 = compact_lr_list[lr_id][z-1]->is_rev;
                                uint32_t tonode1 = (node2 << 1) | (1 - rev2);
                                if(tmp_bbg[node1].edges[1-rev1].count(tonode1) > 0)
                                {
                                    fprintf(stdout, "\tremove %u:%c -> %u:%c\tsupp:%zu\n", node1, "+-"[1-rev1], node2, "+-"[1-rev2], tmp_bbg[node1].edges[1-rev1][tonode1].edge_supp.size());
                                    bbg_remove_edge(node1, 1-rev1, node2, 1-rev2, tmp_bbg);
                                }
                            }
                        }
                    }
                    // 
                    tmp_bbg[contig1].edges[strand1][to_node1] = graph[contig1].edges[strand1][to_node1];
                    tmp_bbg[contig2].edges[1-strand2][to_node2] = graph[contig2].edges[1-strand2][to_node2];
                    // 
                    // bbg_print_graph_gfa(tmp_bbg, contig_list, gopt.out_dir + "/backbone.tmp." + type2str<int>(cnt) + ".gfa");
                    // cnt++;
                    // if(cnt == 1) return;
                }
            }
        }
    }
    graph = tmp_bbg;
}

void bbg_build_graph(vector<BBG_Node_t> &graph, Contig_List_t &contig_list, vector<vector<Align_Seq_t*>> &compact_lr_list)
{
    graph.resize(contig_list.contigs_size);
    for(uint32_t i = 0; i < contig_list.contigs_size; i++) graph[i].contig_id = i;
    for(uint32_t i = 0; i < compact_lr_list.size(); i++)
    {
        if(compact_lr_list[i].size() > 1)
        {
            vector<uint32_t> selected_ids;
            for(uint32_t j = 0; j < compact_lr_list[i].size(); j++)
            {
                uint32_t tid = compact_lr_list[i][j]->t_id;
                if(contig_list.contigs[tid].mean_kmer <= gopt.uniq_freq * (1 + gopt.max_uniq_dev)) // only use unique SRCs for the initial backbone graph
                    selected_ids.push_back(j);
            }
            // add one edge for each pair of consecutive contigs on the long read
            if(selected_ids.size() > 1)
            {
                for(uint32_t j = 0; j < selected_ids.size() - 1; j++)
                    bbg_add_edge(graph, i, 0, compact_lr_list[i], selected_ids[j], selected_ids[j+1]);
            }
        }
    }
}

// also integrate repetitive SRCs into the backbone graph
void bbg_build_graph_3(vector<BBG_Node_t> &graph, Contig_List_t &contig_list, vector<vector<Align_Seq_t*>> &compact_lr_list)
{
    vector<BBG_Node_t> tmp_bbg;
    tmp_bbg.resize(contig_list.contigs_size);
    for(uint32_t i = 0; i < compact_lr_list.size(); i++)
    {
        // fprintf(stderr, "here1 %u\n", i);
        if(compact_lr_list[i].size() > 1)
        {
            // get alignments of unique SRCs
            // fprintf(stderr, "here2\n");
            vector<uint32_t> selected_ids;
            for(uint32_t j = 0; j < compact_lr_list[i].size(); j++)
            {
                uint32_t tid = compact_lr_list[i][j]->t_id;
                if(contig_list.contigs[tid].mean_kmer <= gopt.uniq_freq * (1 + gopt.max_uniq_dev)) // only use unique SRCs for the initial backbone graph
                    selected_ids.push_back(j);
            }
            // add one edge for each pair of consecutive contigs on the long read
            // fprintf(stderr, "here3\n");
            if(selected_ids.size() > 1)
            {
                // fprintf(stderr, "here4\n");
                for(uint32_t j = 0; j < selected_ids.size() - 1; j++)
                    bbg_add_edge(tmp_bbg, i, 0, compact_lr_list[i], selected_ids[j], selected_ids[j+1]);
                uint32_t jfirst = selected_ids.front();
                uint32_t jlast = selected_ids.back();
                // add anything before jfirst
                for(uint32_t j = 0; j < jfirst; j++)
                    bbg_add_edge(tmp_bbg, i, 0, compact_lr_list[i], j, j + 1);
                // add anything after jfirst
                for(uint32_t j = jlast; j < compact_lr_list[i].size() - 1; j++)
                    bbg_add_edge(tmp_bbg, i, 0, compact_lr_list[i], j, j + 1);
            }
            else // either one or zero unique SRCs
            {
                for(uint32_t j = 0; j < compact_lr_list[i].size() - 1; j++)
                    bbg_add_edge(tmp_bbg, i, 0, compact_lr_list[i], j, j + 1);
            }
        }
    }
    // 
    // 
    // int cnt = 0;
    for(uint32_t i = 0; i < graph.size(); i++)
    {
        for(int rev = 0; rev < 2; rev++)
        {
            for(map<uint32_t, BBG_Edge_t>::iterator it = graph[i].edges[rev].begin(); it != graph[i].edges[rev].end(); it++)
            {
                uint32_t contig1 = i;
                uint32_t strand1 = rev;
                uint32_t contig2 = (it->first >> 1);
                uint32_t strand2 = (it->first & 1);
                uint32_t to_node1 = (contig2 << 1) | strand2;
                uint32_t to_node2 = (contig1 << 1) | (1 - strand1);
                if(tmp_bbg[contig1].edges[strand1].count(to_node1) == 0) // if edge doesn't exist in BBG2
                {
                    // update
                    BBG_Edge_t &edge1 = graph[contig1].edges[strand1][to_node1];
                    fprintf(stdout, "check %u:%c -> %u:%c\tsupp:%zu\n", contig1, "+-"[strand1], contig2, "+-"[strand2], edge1.edge_supp.size());
                    for(uint32_t j = 0; j < edge1.edge_supp.size(); j++)
                    {
                        uint32_t lr_id     = edge1.edge_supp[j].lr_id;
                        uint32_t lr_strand = edge1.edge_supp[j].lr_strand;
                        uint32_t cmp_beg   = edge1.edge_supp[j].cmp_head_id;
                        uint32_t cmp_end   = edge1.edge_supp[j].cmp_tail_id;
                        fprintf(stdout, "debug lr:%u strand:%u cmp_beg:%u cmp_end:%u\n", lr_id, lr_strand, cmp_beg, cmp_end);
                        if(lr_strand == 0) // iterate forward
                        {
                            // fprintf(stdout, "here1\n");
                            for(uint32_t z = cmp_beg; z < cmp_end; z++)
                            {
                                // fprintf(stdout, "here2\n");
                                uint32_t node1 = compact_lr_list[lr_id][z]->t_id;
                                uint32_t rev1 = compact_lr_list[lr_id][z]->is_rev;
                                uint32_t node2 = compact_lr_list[lr_id][z+1]->t_id;
                                uint32_t rev2 = compact_lr_list[lr_id][z+1]->is_rev;
                                uint32_t tonode1 = (node2 << 1) | rev2;
                                if(tmp_bbg[node1].edges[rev1].count(tonode1) > 0)
                                {
                                    fprintf(stdout, "\tremove %u:%c -> %u:%c\tsupp:%zu\n", node1, "+-"[rev1], node2, "+-"[rev2], tmp_bbg[node1].edges[rev1][tonode1].edge_supp.size());
                                    bbg_remove_edge(node1, rev1, node2, rev2, tmp_bbg);
                                }
                            }
                        }
                        else  // iterate backward
                        {
                            // fprintf(stdout, "here3\n");
                            for(uint32_t z = cmp_beg; z > cmp_end; z--)
                            {
                                // fprintf(stdout, "here4\n");
                                uint32_t node1 = compact_lr_list[lr_id][z]->t_id;
                                uint32_t rev1 = compact_lr_list[lr_id][z]->is_rev;
                                uint32_t node2 = compact_lr_list[lr_id][z-1]->t_id;
                                uint32_t rev2 = compact_lr_list[lr_id][z-1]->is_rev;
                                uint32_t tonode1 = (node2 << 1) | (1 - rev2);
                                if(tmp_bbg[node1].edges[1-rev1].count(tonode1) > 0)
                                {
                                    fprintf(stdout, "\tremove %u:%c -> %u:%c\tsupp:%zu\n", node1, "+-"[1-rev1], node2, "+-"[1-rev2], tmp_bbg[node1].edges[1-rev1][tonode1].edge_supp.size());
                                    bbg_remove_edge(node1, 1-rev1, node2, 1-rev2, tmp_bbg);
                                }
                            }
                        }
                    }
                    // 
                    tmp_bbg[contig1].edges[strand1][to_node1] = graph[contig1].edges[strand1][to_node1];
                    tmp_bbg[contig2].edges[1-strand2][to_node2] = graph[contig2].edges[1-strand2][to_node2];
                    // 
                    // bbg_print_graph_gfa(tmp_bbg, contig_list, gopt.out_dir + "/backbone.tmp." + type2str<int>(cnt) + ".gfa");
                    // cnt++;
                    // if(cnt == 1) return;
                }
            }
        }
    }
    bbg_remove_weak_edges(tmp_bbg);
    graph = tmp_bbg;
}

// also integrate shorter SRCs into the backbone graph
void bbg_build_graph_2(vector<BBG_Node_t> &graph, uint64_t nb_nodes, vector<vector<Align_Seq_t*>> &compact_lr_list)
{
    graph.resize(nb_nodes);
    for(uint32_t i = 0; i < compact_lr_list.size(); i++)
    {
        // fprintf(stderr, "here1 %u\n", i);
        if(compact_lr_list[i].size() > 1)
        {
            // get alignments of length at least 500 bp
            // fprintf(stderr, "here2\n");
            vector<uint32_t> selected_ids;
            for(uint32_t j = 0; j < compact_lr_list[i].size(); j++)
                if(compact_lr_list[i][j]->n_block >= 500)
                    selected_ids.push_back(j);
            // add one edge for each pair of consecutive contigs on the long read
            // fprintf(stderr, "here3\n");
            if(selected_ids.size() > 1)
            {
                // fprintf(stderr, "here4\n");
                for(uint32_t j = 0; j < selected_ids.size() - 1; j++)
                    bbg_add_edge(graph, i, 0, compact_lr_list[i], selected_ids[j], selected_ids[j+1]);
                // fprintf(stderr, "here8\n");
                uint32_t jfirst = selected_ids.front();
                uint32_t jlast = selected_ids.back();
                selected_ids.clear();
                for(uint32_t j = 0; j < compact_lr_list[i].size(); j++)
                    if(compact_lr_list[i][j]->n_block >= 300)
                        selected_ids.push_back(j);
                // add anything before jfirst
                for(uint32_t j = 0; j < jfirst; j++)
                    bbg_add_edge(graph, i, 0, compact_lr_list[i], selected_ids[j], selected_ids[j+1]);
                // add anything after jfirst
                for(uint32_t j = jlast; j < selected_ids.size() - 1; j++)
                    bbg_add_edge(graph, i, 0, compact_lr_list[i], selected_ids[j], selected_ids[j+1]);
            }
            else
            {
                // fprintf(stderr, "here6\n");
                selected_ids.clear();
                for(uint32_t j = 0; j < compact_lr_list[i].size(); j++)
                    if(compact_lr_list[i][j]->n_block >= 300)
                        selected_ids.push_back(j);
                if(selected_ids.size() > 1)
                {
                    // fprintf(stderr, "here7 %zu\n", selected_ids.size());
                    for(uint32_t j = 0; j < selected_ids.size() - 1; j++)
                        bbg_add_edge(graph, i, 0, compact_lr_list[i], selected_ids[j], selected_ids[j+1]);
                }
            }
        }
    }
}

int bbg_remove_weak_edges(vector<BBG_Node_t> &graph)
{
    int nb_removed = 0;
    for(uint32_t i = 0; i < graph.size(); i++)
    {
        for(int rev1 = 0; rev1 < 2; rev1++) // edge direction: 0 for outgoing and 1 for incoming
        {
            map<uint32_t, BBG_Edge_t>::iterator it;
            for(it = graph[i].edges[rev1].begin(); it != graph[i].edges[rev1].end(); )
            {
                if(it->second.edge_supp.size() < gopt.min_edge_sup)
                {
                    uint32_t node2 = (it->first >> 1);
                    uint8_t rev2 = (it->first & 1);
                    uint32_t to_node2 = (i << 1) | (1 - rev1);
                    it = graph[i].edges[rev1].erase(it); // removing the edge
                    graph[node2].edges[1-rev2].erase(to_node2); // removing the twin edge
                    nb_removed++;
                }
                else
                {
                    it++;
                }
            }
        }
    }
    return nb_removed;
}

// returns false if the simple path is longer than max_depth
bool bbg_find_simple_path_from_source(vector<BBG_Node_t> &graph, uint32_t src_node, uint32_t src_strand, map<uint32_t, BBG_Edge_t>::iterator it, int max_depth, vector<id_strand2_t> &simple_path, float &simple_path_cov)
{
    simple_path.clear();
    simple_path_cov = 0;
    simple_path.push_back({src_strand, src_node});
    uint32_t curr_node   = it->first >> 1;
    uint32_t curr_strand = it->first  & 1;
    int depth = 1;
    while(depth <= max_depth)
    {
        simple_path.push_back({curr_strand, curr_node});
        simple_path_cov += it->second.edge_supp.size();
        if(graph[curr_node].edges[curr_strand].size() == 0) break; // not extendable
        if(graph[curr_node].edges[curr_strand].size() > 1 || graph[curr_node].edges[1 - curr_strand].size() > 1) break; // end of simple path
        it = graph[curr_node].edges[curr_strand].begin();
        curr_node   = it->first >> 1;
        curr_strand = it->first  & 1;
        depth++;
    }
    // 
    if(depth > max_depth) return false;
    // fprintf(stdout, "depth:%d\n", depth);
    simple_path_cov = simple_path_cov / depth;
    return true;
}

bool bbg_find_next_edge(vector<BBG_Node_t> &graph, uint32_t curr_node, uint32_t curr_strand, map<uint32_t, BBG_Edge_t>::iterator &next_it)
{
    if(graph[curr_node].edges[0].size() > 1 || graph[curr_node].edges[1].size() > 1) return false;
    if(curr_strand == 0)
    {
        if(graph[curr_node].edges[0].size() == 1)
        {
            next_it = graph[curr_node].edges[0].begin();
            return true;
        }
        else
        {
            return false;
        }
    }
    else // curr_strand == 1
    {
        if(graph[curr_node].edges[1].size() == 1)
        {
            next_it = graph[curr_node].edges[1].begin();
            return true;
        }
        else
        {
            return false;
        }
    }
}

// this function only extracts simple paths
void bbg_find_simple_paths2(vector<BBG_Node_t> &graph, vector<vector<id_strand2_t>> &simple_paths)
{
    uint32_t i, j;
    uint32_t nb_nodes = graph.size();
    // vector<bool> visited(nb_nodes, false);
    // vector<vector<id_strand2_t>> spaths;
    queue<id_strand2_t> to_explore;
    // 
    for(i = 0; i < nb_nodes; i++)
    {
        // if (graph[i].out_rev.size() == 0 && graph[i].out.size() > 0 && visited[i] == false) // no incoming some outgoing
        if (graph[i].edges[1].size() == 0 && graph[i].edges[0].size() > 0) // no incoming some outgoing
        {
            // fprintf(stdout, "%u outgoing\n", i);
            to_explore.push({0, i});
        }
        // else if (graph[i].out_rev.size() > 0 && graph[i].out.size() == 0 && visited[i] == false) // some incoming no outgoing
        else if (graph[i].edges[1].size() > 0 && graph[i].edges[0].size() == 0) // some incoming no outgoing
        {
            // fprintf(stdout, "%u incoming\n", i);
            to_explore.push({1, i});
        }
    }
    // fprintf(stdout, "\n");
    // fflush(stdout);
    while(to_explore.size() > 0)
    {
        // fprintf(stdout, "\n");
        uint32_t src_node = to_explore.front().id;
        uint32_t src_strand = to_explore.front().strand;
        to_explore.pop();
        // 
        vector<vector<id_strand2_t>> paths_curr; // all simple paths starting from current node
        map<uint32_t, BBG_Edge_t> &edge_list = (src_strand == 0 ? graph[src_node].edges[0] : graph[src_node].edges[1]);
        for(map<uint32_t, BBG_Edge_t>::iterator it = edge_list.begin(); it != edge_list.end(); it++)
        {
            vector<id_strand2_t> path;
            // 
            uint32_t next_node = src_node;
            uint32_t next_strand = src_strand;
            path.push_back({next_strand, next_node});
            // fprintf(stdout, "push %u:%c\n", next_node, "+-"[next_strand]);
            // fflush(stdout);
            map<uint32_t, BBG_Edge_t>::iterator curr_it = it;
            while(true)
            {
                next_node = curr_it->first >> 1;
                next_strand = curr_it->first & 1;
                path.push_back({next_strand, next_node});
                // fprintf(stdout, "push %u:%c\n", next_node, "+-"[next_strand]);
                // fflush(stdout);
                if(bbg_find_next_edge(graph, next_node, next_strand, curr_it) == false) break;
            }
            // 
            if(path.size() > 0)
            {
                paths_curr.push_back(path);
            }
        }
        for(i = 0; i < paths_curr.size(); i++)
        {
            // for(j = 0; j < paths_curr[i].size(); j++)
            // {
            //     fprintf(stdout, "%u:%c\t", paths_curr[i][j].id, "+-"[paths_curr[i][j].strand]);
            // }
            // fprintf(stdout, "\n");
            // fflush(stdout);
            simple_paths.push_back(paths_curr[i]);
            for(j = 0; j < paths_curr[i].size() - 1; j++)
            {
                // fprintf(stdout, "removing %u:%c -> %u:%c\n", node1, "+-"[rev1], node2, "+-"[rev2]);
                // fflush(stdout);
                bbg_remove_edge(paths_curr[i][j].id, paths_curr[i][j].strand, paths_curr[i][j+1].id, paths_curr[i][j+1].strand, graph);
                // fprintf(stdout, "\n");
                // fflush(stdout);
            }
            uint32_t last_node = paths_curr[i].back().id;
            uint32_t last_strand = paths_curr[i].back().strand;
            if(last_strand == 0)
            {
                if(graph[last_node].edges[0].size() > 0 && graph[last_node].edges[1].size() == 0)
                    to_explore.push({last_strand, last_node});
            }
            else
            {
                if(graph[last_node].edges[0].size() == 0 && graph[last_node].edges[1].size() > 0)
                    to_explore.push({last_strand, last_node});
            }
        }
    }
    // 
    fprintf(stdout, "\n");
    for(i = 0; i < simple_paths.size(); i++)
    {
        fprintf(stdout, "simple_path:%u size:%zu\t", i, simple_paths[i].size());
        for(j = 0; j < simple_paths[i].size(); j++)
        {
            fprintf(stdout, "%u:%c\t", simple_paths[i][j].id, "+-"[simple_paths[i][j].strand]);
        }
        fprintf(stdout, "\n");
        fflush(stdout);
    }
    // fclose(flog);
}

// prints the backbone graph in GFA format
void bbg_print_graph_gfa(vector<BBG_Node_t> &graph, Contig_List_t &contig_list, string path)
{
    FILE *fp = fopen(path.c_str(), "w");
    if(fp == NULL)
    {
        fprintf(stderr, "[ERROR] (Backbone_graph::bbg_pring_graph_gfa) cannot open file: %s\n", path.c_str());
        exit(EXIT_FAILURE);
    }

    // // print all nodes
    // for(uint32_t i = 0; i < contig_list.contigs_size; i++)
    // {
    //     fprintf(fp, "S\t%u\t%s\n", i, string(contig_list.contigs[i].len, 'N').c_str());
    // }

    // print only nodes that appear on some edges
    set<uint32_t> to_print;
    for(uint32_t i = 0; i < graph.size(); i++)
    {
        for(int rev = 0; rev < 2; rev++)
        {
            for(map<uint32_t, BBG_Edge_t>::iterator it = graph[i].edges[rev].begin(); it != graph[i].edges[rev].end(); it++)
            {
                to_print.insert(i);
                to_print.insert((it->first >> 1));
            }
        }
    }
    for(set<uint32_t>::iterator it = to_print.begin(); it != to_print.end(); it++)
    {
        uint32_t cid = graph[*it].contig_id;
        string contig_str = get_uncompressed_dna(contig_list.contigs[cid].comp_seq, contig_list.contigs[cid].len, contig_list.contigs[cid].comp_len);
        fprintf(fp, "S\t%u\t%s\tLN:i:%zu\tKC:i:%u\n", *it, contig_str.c_str(), contig_str.size(), contig_list.contigs[cid].kmer_count);
    }

    // print links
    for(uint32_t i = 0; i < graph.size(); i++)
    {
        for(int rev = 0; rev < 2; rev++)
        {
            for(map<uint32_t, BBG_Edge_t>::iterator it = graph[i].edges[rev].begin(); it != graph[i].edges[rev].end(); it++)
            {
                fprintf(fp, "L\t%u\t%c\t%u\t%c\t0M\n", i, "+-"[rev] , (it->first >> 1), ((it->first & 1) ? '-' : '+'));
            }
        }
    }

    fclose(fp);
}

bool sort_tuple_desc(const tuple<uint32_t, uint32_t, uint32_t>& a, const tuple<uint32_t, uint32_t, uint32_t>& b) 
{ 
    return (get<0>(a) > get<0>(b)); 
}

void bbg_general_stats(vector<BBG_Node_t> &graph, Contig_List_t &contig_list, string logpath)
{
    FILE*fp = file_open_write(logpath.c_str());
    // fprintf(stdout, "##### graph statistics: %s\n", graph_name.c_str());
    // 
    uint32_t i;
    uint32_t num = graph.size();
    uint32_t nb_node = 0;
    uint32_t nb_edge = 0;
    for(i = 0; i < num; i++)
    {
        nb_node += (graph[i].edges[0].size() > 0 || graph[i].edges[1].size() > 0);
        nb_edge += graph[i].edges[0].size() + graph[i].edges[1].size();
    }
    fprintf(fp, "nodes: %d\n", nb_node);
    fprintf(fp, "edges: %d\n", nb_edge/2);
    // 
    vector<bool> visited(num, false);
    vector<tuple<uint32_t, uint32_t, uint32_t>> component_list;
    for(i = 0; i < num; i++)
    {
        if(visited[i] == false && (graph[i].edges[0].size() > 0 || graph[i].edges[1].size() > 0))
        {
            uint64_t cc_size = 0;
            uint64_t cc_node = 0;
            queue<uint32_t> q;
            q.push(i);
            cc_node++;
            cc_size += contig_list.contigs[graph[i].contig_id].len;
            visited[i] = true;
            // fprintf(stdout, "new_component:%zu\tstart_node:%u\n", component_list.size(), i);
            while(q.size() > 0)
            {
                uint32_t curr_id = q.front();
                q.pop();
                // fprintf(stdout, "\tpopped %u\n", curr_id);
                for(int rev = 0; rev < 2; rev++)
                {
                    for(map<uint32_t, BBG_Edge_t>::iterator it = graph[curr_id].edges[rev].begin(); it != graph[curr_id].edges[rev].end(); it++)
                    {
                        uint32_t next_id = it->first >> 1;
                        if(visited[next_id] == false)
                        {
                            q.push(next_id);
                            cc_node++;
                            cc_size += contig_list.contigs[graph[next_id].contig_id].len;
                            visited[next_id] = true;
                            // fprintf(stdout, "\t\tpushed %u\n", next_id);
                        }
                    }
                }
            }
            component_list.push_back(make_tuple(cc_size, cc_node, i));
        }
    }
    // 
    sort(component_list.begin(), component_list.end(), sort_tuple_desc);
    fprintf(fp, "connected_components: %zu\n", component_list.size());
    for(i = 0; i < component_list.size(); i++)
    {
        fprintf(fp, "\tcomponent:%u\tsize:%u\tnodes:%u\trepresentative:%u\n", i, get<0>(component_list[i]), get<1>(component_list[i]), get<2>(component_list[i]));
    }
    // fprintf(stdout, "\n");
    fclose(fp);
}

// // prints graph nodes and lists
// void bbg_print_graph_node_list(vector<Node_t> &graph, Contig_List_t &contig_list, string path)
// {
//     FILE *fp = fopen(path.c_str(), "w");
//     if(fp == NULL)
//     {
//         fprintf(stderr, "[ERROR] cannot open file: %s\n", path.c_str());
//         exit(EXIT_FAILURE);
//     }
//     // print all nodes
//     for(uint32_t i = 0; i < contig_list.contigs_size; i++)
//     {
//         fprintf(fp, ">%u\tlen:%u\tlrs:%lu\n", i, contig_list.contigs[i].len, graph[i].lrs.size());
//         for(set<id_strand_t>::iterator it = graph[i].lrs.begin(); it != graph[i].lrs.end(); it++)
//             fprintf(fp, "    %u:%c", it->id, (it->strand ? '-' : '+'));
//         fprintf(fp, "\n");
//     }

//     fclose(fp);
// }

void bbg_report_branching_nodes(vector<BBG_Node_t> &graph, string logpath)
{
    FILE *fp = file_open_write(logpath.c_str());
    // 
    uint32_t num = graph.size();
    uint32_t i;
    for(i = 0; i < num; i++)
    {
        if(graph[i].edges[0].size() >= 2 || graph[i].edges[1].size() >= 2) // cannot be the source of a bubble
            fprintf(fp, "node:%u\tincoming:%zu\toutgoing:%zu\n", i, graph[i].edges[0].size(), graph[i].edges[1].size());
    }
    fclose(fp);
}

// void bbg_identify_unused_longreads(Longread_List_t &lr_list, vector<vector<Align_Seq_t*>> &compact_lr_list, vector<vector<path_elem_t>> &simple_paths, vector<vector<id_strand_t>> &map_contig2lr_uniq, vector<uint8_t> &unused)
// {
//     uint32_t i, j;
//     uint32_t contig1;
//     unused.resize(compact_lr_list.size(), 1);
//     // mark reads used for uniq assembly
//     for(i = 0; i < simple_paths.size(); i++)
//     {
//         for(j = 0; j < simple_paths[i].size(); j++)
//         {
//             contig1 = simple_paths[i][j].id;
//             for(uint32_t x = 0; x < map_contig2lr_uniq[contig1].size(); x++)
//                 unused[map_contig2lr_uniq[contig1][x].id] = 0;
//         }
//         // include reads covering two extremeties
//         contig1 = simple_paths[i].front().id;
//         for(uint32_t x = 0; x < map_contig2lr_uniq[contig1].size(); x++)
//             unused[map_contig2lr_uniq[contig1][x].id] = 2;
//         // 
//         contig1 = simple_paths[i].back().id;
//         for(uint32_t x = 0; x < map_contig2lr_uniq[contig1].size(); x++)
//             unused[map_contig2lr_uniq[contig1][x].id] = 2;
//     }
//     // // mark reads with only one alignment to them
//     // for(i = 0; i < compact_lr_list.size(); i++)
//     // {
//     //     // if(lr_list.reads[i].contig_aln_num <= 1)
//     //     if(compact_lr_list[i].size() <= 1)
//     //     {
//     //         unused[i] = false;
//     //     }
//     // }
//     // 
//     FILE *fp = fopen((gopt.outDir + "/unused.fasta").c_str(), "wt");
//     if(fp == NULL)
//     {
//         fprintf(stderr, "[ERROR] could not open file: %s\n", (gopt.outDir + "/unused.fasta").c_str());
//         exit(EXIT_FAILURE);
//     }
//     // for(i = 0; i < lr_list.reads_size; i++)
//     for(i = 0; i < compact_lr_list.size(); i++)
//     {
//         if(unused[i])
//         {
//             string tmp_seq = get_uncompressed_dna(lr_list.reads[i].comp_seq, lr_list.reads[i].len, lr_list.reads[i].comp_len);
//             fprintf(fp, ">%u %s\n%s\n", i, (unused[i] == 2 ? "tail" : ""), tmp_seq.c_str());
//         }
//     }
//     fclose(fp);
//     // 
//     fp = fopen((gopt.outDir + "/unused.txt").c_str(), "wt");
//     if(fp == NULL)
//     {
//         fprintf(stderr, "[ERROR] could not open file: %s\n", (gopt.outDir + "/unused.txt").c_str());
//         exit(EXIT_FAILURE);
//     }
//     // for(i = 0; i < lr_list.reads_size; i++)
//     for(i = 0; i < compact_lr_list.size(); i++)
//     {
//         // fprintf(fp, "%u\t%u\n", i, (unused[i] ? 1 : 0));
//         fprintf(fp, "%u\t%u\n", i, unused[i]);
//     }
//     fclose(fp);
// }

// bool bbg_find_next_node(vector<Node_t> &graph, uint32_t curr_node, uint32_t curr_strand, uint32_t &next_node, uint32_t &next_strand)
// {
//     if(curr_strand == 0)
//     {
//         if(graph[curr_node].out.size() == 1)
//         {
//             // fprintf(stdout, "here1\n");
//             // fflush(stdout);
//             next_node = (graph[curr_node].out.begin()->first >> 1);
//             next_strand = (graph[curr_node].out.begin()->first & 1);
//             return true;
//         }
//         else
//         {
//             // fprintf(stdout, "here2\n");
//             // fflush(stdout);
//             return false;
//         }
//     }
//     else // curr_strand == 1
//     {
//         if(graph[curr_node].out_rev.size() == 1)
//         {
//             // fprintf(stdout, "here3\n");
//             // fflush(stdout);
//             next_node = (graph[curr_node].out_rev.begin()->first >> 1);
//             next_strand = (graph[curr_node].out_rev.begin()->first & 1);
//             return true;
//         }
//         else
//         {
//             // fprintf(stdout, "here4\n");
//             // fflush(stdout);
//             return false;
//         }
//     }
// }

// // this function only extracts simple paths
// // void asm_traverse_untangling_graph()
// void bbg_find_simple_paths(Contig_List_t &contig_list, vector<Node_t> &graph_uniq, vector<vector<path_elem_t>> &simple_paths)
// {
//     uint32_t i;
//     bool traverse_outgoing;
//     bool traverse_incoming;
//     vector<bool> visited(contig_list.contigs_size, false);
//     // 
//     while (true)
//     {
//         traverse_outgoing = false;
//         traverse_incoming = false;
//         for(i = 0; i < contig_list.contigs_size; i++)
//         {
//             // TODO: what if out-degree is more than 1? should I traverse?
//             if (graph_uniq[i].out.size() == 1 && graph_uniq[i].out_rev.size() == 0 && visited[i] == false) // node is a tip (leaf) with outgoing edges
//             {
//                 traverse_outgoing = true;
//                 break;
//             }
//             // TODO: what if in-degree is more than 1? should I traverse?
//             else if (graph_uniq[i].out.size() == 0 && graph_uniq[i].out_rev.size() == 1 && visited[i] == false) // node is a tip (leaf) with incoming edges
//             {
//                 traverse_incoming = true;
//                 break;
//             }
//             else if (graph_uniq[i].out.size() == 1 && graph_uniq[i].out_rev.size() == 1 && visited[i] == false) // node is on a simple path
//             {
//                 traverse_outgoing = true;
//                 traverse_incoming = true;
//                 break;
//             }
//         }
//         if(i == contig_list.contigs_size) // no other node to start traversing from
//             break;
//         fprintf(stdout, "new_traversal node:%u outgoing:%zu incoming:%zu\n", i, graph_uniq[i].out.size(), graph_uniq[i].out_rev.size());
//         visited[i] = true;
//         vector<id_strand_t> path_outgoing;
//         if(traverse_outgoing == true)
//         {
//             uint32_t curr_node = i;
//             uint32_t curr_strand = 0;
//             uint32_t next_node, next_strand;
//             bool found_next;
//             while(true)
//             {
//                 found_next = bbg_find_next_node(graph_uniq, curr_node, curr_strand, next_node, next_strand);
//                 if(found_next == false || visited[next_node] == true)
//                     break;
//                 curr_node = next_node;
//                 curr_strand = next_strand;
//                 visited[curr_node] = true;
//                 path_outgoing.push_back({curr_strand, 0, curr_node});
//             }
//         }
//         fprintf(stdout, "==== path_outgoing ");
//         for(size_t j = 0; j < path_outgoing.size(); j++)
//             fprintf(stdout, "%u:%c ", path_outgoing[j].id, "+-"[path_outgoing[j].strand]);
//         fprintf(stdout, "\n");
//         vector<id_strand_t> path_incoming;
//         if(traverse_incoming == true)
//         {
//             uint32_t curr_node = i;
//             uint32_t curr_strand = 1;
//             uint32_t next_node, next_strand;
//             bool found_next;
//             while(true)
//             {
//                 found_next = bbg_find_next_node(graph_uniq, curr_node, curr_strand, next_node, next_strand);
//                 if(found_next == false || visited[next_node] == true)
//                     break;
//                 curr_node = next_node;
//                 curr_strand = next_strand;
//                 visited[curr_node] = true;
//                 path_incoming.push_back({curr_strand, 0, curr_node});
//             }
//         }
//         fprintf(stdout, "==== path_incoming ");
//         for(size_t j = 0; j < path_incoming.size(); j++)
//             fprintf(stdout, "%u:%c ", path_incoming[j].id, "+-"[path_incoming[j].strand]);
//         fprintf(stdout, "\n");
//         // 
//         vector<path_elem_t> path;
//         // add the incoming path
//         for (vector<id_strand_t>::reverse_iterator it = path_incoming.rbegin(); it != path_incoming.rend(); ++it)
//             path.push_back({(uint32_t)1 - it->strand, it->id, set<uint32_t>(), map<uint32_t, int32_t>()});
//         // add the node
//         path.push_back({0, i, set<uint32_t>(), map<uint32_t, int32_t>()});
//         // add the outgoing path
//         for (vector<id_strand_t>::iterator it = path_outgoing.begin(); it != path_outgoing.end(); ++it)
//             path.push_back({it->strand, it->id, set<uint32_t>(), map<uint32_t, int32_t>()});
//         // print the path
//         if(path.size() > 1)
//             fprintf(stdout, "==== path %u:%c -> %u:%c\n", path.front().id, "+-"[path.front().strand], path.back().id, "+-"[path.back().strand]);
//         else
//             fprintf(stdout, "==== path %u:%c\n", path.front().id, "+-"[path.front().strand]);
//         for(size_t j = 0; j < path.size(); j++)
//             fprintf(stdout, "%u:%c ", path[j].id, "+-"[path[j].strand]);
//         fprintf(stdout, "\n\n");
//         // 
//         if(path.size() > 0)
//             simple_paths.push_back(path);
//     }
//     // print un-visited nodes
//     fprintf(stdout, "#### total num of paths: %zu\n", simple_paths.size());
//     fprintf(stdout, "#### unvisited nodes: ");
//     for(uint32_t j = 0; j < contig_list.contigs_size; j++)
//         if(visited[j] == false && (graph_uniq[j].out.size() > 0 || graph_uniq[j].out_rev.size() > 0))
//             fprintf(stdout, "%u ", j);
//     fprintf(stdout, "\n\n");
// }