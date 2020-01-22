// /*****************************************************
//  * Author: Ehsan Haghshenas (ehaghshe AT sfu DOT ca) *
//  *****************************************************/

#ifndef __BACKBONE_GRAPH__
#define __BACKBONE_GRAPH__

#include "Common.hpp"
#include "Longread.hpp"
#include <vector>
#include <set>
#include <map>

using namespace std;

// typedef struct
// {
//     uint32_t lr_id;
//     uint32_t index1; // index of the compact element corresponding to from_node
//     uint32_t index2; // index of the compact element corresponding to to_node
// } edge_supp_t;

typedef struct
{
    uint32_t lr_id:31;
    uint32_t lr_strand:1;
    uint32_t cmp_head_id; // index of the compact element corresponding to head node
    uint32_t cmp_tail_id; // index of the compact element corresponding to tail node
} Edge_Supp_t;

typedef struct
{
    uint32_t lr_id:31;
    uint32_t lr_strand:1;
    uint32_t spos;
    uint32_t epos;
} Consensus_Supp_t;

typedef struct // edge between head and tail
{
    uint32_t head_end; // the end position of the shared region on the head node
    uint32_t tail_beg; // the beginning position of the shared region on the tail node
    uint32_t flag:8; // used for traversal
    // uint32_t len_estimate:24;
    string   cns_seq;
    vector<Edge_Supp_t> edge_supp;
    vector<Consensus_Supp_t> cns_supp; // if cns_supp.size() == 0 the path should be broken at this anchor pair
} BBG_Edge_t;

typedef struct
{
    map<uint32_t, BBG_Edge_t> edges[2]; // 0 for outgoing edges and 1 for incoming edges
} BBG_Node_t;

void bbg_add_edge_with_supp(vector<BBG_Node_t> &graph, uint32_t node1, uint8_t rev1, uint32_t node2, uint8_t rev2, vector<Edge_Supp_t> &shared_supp);
void bbg_remove_edge(uint32_t node1, uint8_t rev1, uint32_t node2, uint8_t rev2, vector<BBG_Node_t> &graph);
BBG_Edge_t* bbg_get_edge(uint32_t node1, uint8_t rev1, uint32_t node2, uint8_t rev2, vector<BBG_Node_t> &graph);
void bbg_add_repetitive_SRCs(vector<BBG_Node_t> &graph, Contig_List_t &contig_list, vector<vector<Align_Seq_t*>> &compact_lr_list);
void bbg_build_graph(vector<BBG_Node_t> &graph, Contig_List_t &contig_list, vector<vector<Align_Seq_t*>> &compact_lr_list);
void bbg_build_graph_2(vector<BBG_Node_t> &graph, uint64_t nb_nodes, vector<vector<Align_Seq_t*>> &compact_lr_list);
void bbg_build_graph_3(vector<BBG_Node_t> &graph, Contig_List_t &contig_list, vector<vector<Align_Seq_t*>> &compact_lr_list);
int bbg_remove_weak_edges(vector<BBG_Node_t> &graph);
void bbg_print_graph_gfa(vector<BBG_Node_t> &graph, Contig_List_t &contig_list, string path);
void bbg_general_stats(vector<BBG_Node_t> &graph, Contig_List_t &contig_list, string path);
// void bbg_print_graph_node_list(vector<Node_t> &graph, Contig_List_t &contig_list, string path);
void bbg_report_branching_nodes(vector<BBG_Node_t> &graph, string logpath);
bool bbg_find_simple_path_from_source(vector<BBG_Node_t> &graph, uint32_t src_node, uint32_t src_strand, map<uint32_t, BBG_Edge_t>::iterator it, int max_depth, vector<id_strand2_t> &simple_path, float &simple_path_cov);
// void bbg_find_simple_paths(Contig_List_t &contig_list, vector<Node_t> &graph_uniq, vector<vector<path_elem_t>> &simple_paths);
void bbg_find_simple_paths2(vector<BBG_Node_t> &graph, vector<vector<id_strand2_t>> &simple_paths);
// void bbg_identify_unused_longreads(Longread_List_t &lr_list, vector<vector<Align_Seq_t*>> &compact_lr_list, vector<vector<path_elem_t>> &simple_paths, vector<vector<id_strand_t>> &map_contig2lr_uniq, vector<uint8_t> &unused);

#endif // __BACKBONE_GRAPH__
