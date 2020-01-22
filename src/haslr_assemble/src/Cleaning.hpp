/*****************************************************
 * Author: Ehsan Haghshenas (ehaghshe AT sfu DOT ca) *
 *****************************************************/

#ifndef __CLEANING__
#define __CLEANING__

#include "Common.hpp"
#include "Backbone_graph.hpp"

// int remove_tips_old(vector<BBG_Node_t> &graph, int max_len, string logpath);
int clean_tips(vector<BBG_Node_t> &graph, int max_depth, string logpath);
int clean_small_bubbles(vector<BBG_Node_t> &graph, string logpath);
int clean_simple_bubbles_old(vector<BBG_Node_t> &graph, int max_depth, string logpath);
int clean_simple_bubbles(vector<BBG_Node_t> &graph, int max_depth, string logpath);
int clean_super_bubbles(vector<BBG_Node_t> &graph, int max_dist, string logpath);

#endif // __CLEANING__
