// /*****************************************************
//  * Author: Ehsan Haghshenas (ehaghshe AT sfu DOT ca) *
//  *****************************************************/

#ifndef __ASSEMBLE__
#define __ASSEMBLE__

#include "Common.hpp"
#include "Longread.hpp"
#include "Cleaning.hpp"

// struct cns_seq_t
// {
//     uint32_t prev_end; // position on the previous conting that consensus sequence starts
//     uint32_t next_start; // position on the next conting that consensus sequence ends
//     uint32_t count; // the number of sequences used to generate the consensus sequence
//     string cns; // the consensus sequence obtained by SPOA
// };

void asm_calc_edge_coordinates_MT(vector<BBG_Node_t> *bbgraph, Contig_List_t *contig_list, Longread_List_t *lr_list, vector<vector<Align_Seq_t*>> *compact_lr_list, string log_path);
void asm_cal_cns_seq_MT(vector<BBG_Node_t> *bbgraph, Longread_List_t *lr_list, string log_path);
void asm_get_assembly(vector<BBG_Node_t> &graph);

// void asm_init(Contig_List_t *contig_list, Longread_List_t *lr_list, vector<vector<Align_Seq_t*>> *compact_lr_list, vector<vector<path_elem_t>> *all_paths);
// void asm_assemble_paths(string asm_path, string ann_path);
// void asm_assemble_paths(Contig_List_t *contig_list, Longread_List_t *lr_list, vector<vector<Align_Seq_t*>> *compact_lr_list, 
//                         vector<vector<id_strand2_t>> *all_paths, vector<vector<cns_info_t>> *cns_list, string asm_path, string ann_path);

#endif // __ASSEMBLE__
