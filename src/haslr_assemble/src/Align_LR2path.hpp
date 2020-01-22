// /*****************************************************
//  * Author: Ehsan Haghshenas (ehaghshe AT sfu DOT ca) *
//  *****************************************************/

#ifndef __ALIGN_LR2PATH__
#define __ALIGN_LR2PATH__

#include "Common.hpp"
#include "Longread.hpp"
#include "Graph_repeat.hpp"
#include <vector>

using namespace std;

void path_uniq_align_lrs(vector<vector<path_elem_t>> &simple_paths, vector<vector<Align_Seq_t*>> &compact_lr_list, vector<vector<id_strand_t>> &map_contig2lr);
void path_repeat_align_all_lrs(vector<vector<path_elem_t>> &ov_paths, vector<vector<uint32_t>> &covering_lr, vector<vector<Align_Seq_t*>> &compact_lr_list);
void bridge_between_simple_paths(vector<vector<path_elem_t>> &simple_paths, vector<vector<path_elem_t>> &ov_paths, vector<vector<Align_Seq_t*>> &compact_lr_list);

#endif // __ALIGN_LR2PATH__
