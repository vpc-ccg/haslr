// /*****************************************************
//  * Author: Ehsan Haghshenas (ehaghshe AT sfu DOT ca) *
//  *****************************************************/

#ifndef __GRAPH_REPEAT__
#define __GRAPH_REPEAT__

#include "Common.hpp"
#include "Longread.hpp"
#include <vector>
#include <deque>
#include <set>
#include <map>

using namespace std;

struct Ov_edge_t
{
    uint8_t is_trasitive:1;
    uint8_t aln_rev1:1; // alignment transcript is made with reverse complement of the source compressed read 
    uint8_t aln_rev2:1; // alignment transcript is made with reverse complement of the target compressed read 
    deque<int32_t> aln_cmp1; // alignment transcript of the source compressed read
    deque<int32_t> aln_cmp2; // alignment transcript of the target compressed read
};

struct Ov_Node_t // overlap graph node
{
    uint8_t is_contained;
    vector<uint32_t> containing_list;
    map<uint32_t, Ov_edge_t> out;
    map<uint32_t, Ov_edge_t> out_rev;
    // NOTE: the value shows whether the edge is removed during transitive reduction (1) or not (0)
};

void asm_build_ovgraph_from_unused_lrs(Longread_List_t &lr_list, vector<vector<Align_Seq_t*>> &compact_lr_list, vector<vector<id_strand_t>> &map_contig2lr, vector<vector<path_elem_t>> &simple_paths, vector<vector<id_strand_t>> &map_contig2lr_uniq, vector<uint8_t> &unused, vector<vector<path_elem_t>> &ov_paths, vector<vector<uint32_t>> &covering_lr);

#endif // __GRAPH_REPEAT__
