/*****************************************************
 * Author: Ehsan Haghshenas (ehaghshe AT sfu DOT ca) *
 *****************************************************/

#include "Cleaning.hpp"

int clean_small_bubbles(vector<BBG_Node_t> &graph, string logpath)
{
    FILE *fp = file_open_write(logpath.c_str());
    // 
    uint32_t i;
    int nb_removed = 0;
    map<uint32_t, BBG_Edge_t>::iterator it_out;
    map<uint32_t, BBG_Edge_t>::iterator it_in;
    for(i = 0; i < graph.size(); i++)
    {
        if(graph[i].edges[1].size() > 0 && graph[i].edges[0].size() > 0)
        {
            bool bubble_detected = false;
            for(it_in = graph[i].edges[1].begin(); it_in != graph[i].edges[1].end(); it_in++)
            {
                for(it_out = graph[i].edges[0].begin(); it_out != graph[i].edges[0].end(); it_out++)
                {
                    uint32_t node1 = (it_in->first >> 1);
                    uint32_t rev1  = (it_in->first  & 1);
                    uint32_t to_node = it_out->first;
                    uint32_t node2 = (to_node >> 1);
                    uint32_t rev2  = (to_node & 1);
                    if(graph[node1].edges[1-rev1].count(to_node) > 0) // small bubble detected
                    {
                        double short_path_cov = graph[node1].edges[1-rev1][to_node].edge_supp.size();
                        double long_path_cov  = (it_in->second.edge_supp.size() + it_out->second.edge_supp.size())/2.0;
                        fprintf(fp, "small_bubble cov:%.2lf %u:%c -> %u:%c\n", short_path_cov, node1, "+-"[1-rev1], node2, "+-"[rev2]);
                        fprintf(fp, "             cov:%.2lf %u:%c -> %u:%c -> %u:%c\n", long_path_cov, node1, "+-"[1-rev1], i, "+-"[0], node2, "+-"[rev2]);
                        if(short_path_cov < long_path_cov)
                        {
                            // fprintf(fp, "             removing short path\n");
                            bbg_remove_edge(node1, 1 - rev1, node2, rev2, graph);
                        }
                        else
                        {
                            // fprintf(fp, "             removing long path\n");
                            bbg_remove_edge(node1, 1 - rev1, i, 0, graph);
                            bbg_remove_edge(i, 0, node2, rev2, graph);
                        }
                        nb_removed++;
                        bubble_detected = true;
                    }
                    if(bubble_detected) break;
                }
                if(bubble_detected) break;
            }
        }
    }
    fclose(fp);
    return nb_removed;
}

int clean_tips(vector<BBG_Node_t> &graph, int max_depth, string logpath)
{
    FILE *fp;
    if(max_depth == 1)
        fp = file_open_write(logpath.c_str());
    else
        fp = file_open_append(logpath.c_str());
    // 
    int nb_removed = 0;
    uint32_t num = graph.size();
    uint32_t i, j;
    for(i = 0; i < num; i++)
    {
        uint32_t src_strand;
        if(graph[i].edges[1].size() == 0 && graph[i].edges[0].size() == 1) // no incoming one outgoing
            src_strand = 0;
        else if(graph[i].edges[1].size() == 1 && graph[i].edges[0].size() == 0) // one incoming no outgoing
            src_strand = 1;
        else continue;
        // 
        vector<id_strand2_t> path1;
        float cov1;
        if(bbg_find_simple_path_from_source(graph, i, src_strand, graph[i].edges[src_strand].begin(), max_depth, path1, cov1))
        {
            if(graph[path1.back().id].edges[path1.back().strand].size() == 0) continue;
            fprintf(fp, "tip_len:%zu\t%u:%c -> %u:%c\n", path1.size() - 1, path1.front().id, "+-"[path1.front().strand], path1.back().id, "+-"[path1.back().strand]);
            for(j = 0; j < path1.size() - 1; j++)
                bbg_remove_edge(path1[j].id, path1[j].strand, path1[j+1].id, path1[j+1].strand, graph);
            nb_removed++;
        }
        // else
        // {
        //     fprintf(fp, "tip_src %u:%c %zu long\n", i, "+-"[src_strand], path1.size());
        // }
    }
    fclose(fp);
    return nb_removed;
}

int clean_simple_bubbles_old(vector<BBG_Node_t> &graph, int max_depth, string logpath)
{
    FILE *fp = file_open_write(logpath.c_str());
    // 
    int nb_removed = 0;
    uint32_t num = graph.size();
    uint32_t i, j;
    for(i = 0; i < num; i++)
    {
        if(graph[i].edges[0].size() < 2 && graph[i].edges[1].size() < 2) // cannot be the source of a bubble
            continue;
        if(graph[i].edges[0].size() == 2) // detect bubble at the outgoing edges
        {
            vector<id_strand2_t> path1, path2;
            float cov1, cov2;
            bool found1, found2;
            map<uint32_t, BBG_Edge_t>::iterator it = graph[i].edges[0].begin();
            found1 = bbg_find_simple_path_from_source(graph, i, 0, it, max_depth, path1, cov1);
            it++;
            found2 = bbg_find_simple_path_from_source(graph, i, 0, it, max_depth, path2, cov2);
            if(found1 && found2)
            {
                // check end nodes
                if(path1.back().id == path2.back().id && path1.back().strand == path2.back().strand)
                {
                    fprintf(fp, "simple_bubble cov:%.2lf ", cov1);
                    for(j = 0; j < path1.size(); j++) fprintf(fp, "%u:%c ", path1[j].id, "+-"[path1[j].strand]);
                    fprintf(fp, "\n              cov:%.2lf ", cov2);
                    for(j = 0; j < path2.size(); j++) fprintf(fp, "%u:%c ", path2[j].id, "+-"[path2[j].strand]);
                    fprintf(fp, "\n");
                    if(cov1 < cov2)
                    {
                        for(j = 0; j < path1.size() - 1; j++)
                            bbg_remove_edge(path1[j].id, path1[j].strand, path1[j+1].id, path1[j+1].strand, graph);
                    }
                    else
                    {
                        for(j = 0; j < path2.size() - 1; j++)
                            bbg_remove_edge(path2[j].id, path2[j].strand, path2[j+1].id, path2[j+1].strand, graph);
                    }
                    nb_removed++;
                    // 
                    i--;
                    continue;
                }
            }
        }
        if(graph[i].edges[1].size() == 2) // detect bubble at the incoming edges
        {
            vector<id_strand2_t> path1, path2;
            float cov1, cov2;
            bool found1, found2;
            map<uint32_t, BBG_Edge_t>::iterator it = graph[i].edges[1].begin();
            found1 = bbg_find_simple_path_from_source(graph, i, 1, it, max_depth, path1, cov1);
            it++;
            found2 = bbg_find_simple_path_from_source(graph, i, 1, it, max_depth, path2, cov2);
            if(found1 && found2)
            {
                // check end nodes
                if(path1.back().id == path2.back().id && path1.back().strand == path2.back().strand)
                {
                    fprintf(fp, "simple_bubble cov:%.2lf ", cov1);
                    for(j = 0; j < path1.size(); j++) fprintf(fp, "%u:%c ", path1[j].id, "+-"[path1[j].strand]);
                    fprintf(fp, "\n              cov:%.2lf ", cov2);
                    for(j = 0; j < path2.size(); j++) fprintf(fp, "%u:%c ", path2[j].id, "+-"[path2[j].strand]);
                    fprintf(fp, "\n");
                    if(cov1 < cov2)
                    {
                        for(j = 0; j < path1.size() - 1; j++)
                            bbg_remove_edge(path1[j].id, path1[j].strand, path1[j+1].id, path1[j+1].strand, graph);
                    }
                    else
                    {
                        for(j = 0; j < path2.size() - 1; j++)
                            bbg_remove_edge(path2[j].id, path2[j].strand, path2[j+1].id, path2[j+1].strand, graph);
                    }
                    nb_removed++;
                    // 
                    i--;
                    continue;
                }
            }
        }
    }
    fclose(fp);
    return nb_removed;
}

bool compare_Edge_Supp(const Edge_Supp_t &a, const Edge_Supp_t &b)
{
    return (a.lr_id < b.lr_id);
}

void get_shared_lr_supp(vector<Edge_Supp_t> &edge1_supp, vector<Edge_Supp_t> &edge2_supp, vector<Edge_Supp_t> &shared_supp)
{
    uint32_t i, j;
    // for sanity checking: check that both edge1_supp and edge2_supp are sorted
    // TODO: could be disabled
    for(i = 0; i < edge1_supp.size() - 1; i++)
    {
        if(edge1_supp[i].lr_id > edge1_supp[i+1].lr_id)
        {
            fprintf(stderr, "[ERROR] (cleaning::get_shared_lr_supp) the list of long read supports is not sorted!\n");
            exit(EXIT_FAILURE);
        }
    }
    for(j = 0; j < edge2_supp.size() - 1; j++)
    {
        if(edge2_supp[j].lr_id > edge2_supp[j+1].lr_id)
        {
            fprintf(stderr, "[ERROR] (cleaning::get_shared_lr_supp) the list of long read supports is not sorted!\n");
            exit(EXIT_FAILURE);
        }
    }
    // sort(edge1_supp.begin(), edge1_supp.end(), compare_Edge_Supp);
    // sort(edge2_supp.begin(), edge2_supp.end(), compare_Edge_Supp);
    // 
    i = 0;
    j = 0;
    while(i < edge1_supp.size() && j < edge2_supp.size())
    {
        if(edge1_supp[i].lr_id == edge2_supp[j].lr_id)
        {
            // sanity check
            if(edge1_supp[i].lr_strand != edge2_supp[j].lr_strand)
            {
                fprintf(stderr, "[ERROR] (cleaning::get_shared_lr_supp) same supporting long read has different strand!\n");
                exit(EXIT_FAILURE);
            }
            shared_supp.push_back({edge1_supp[i].lr_id, edge1_supp[i].lr_strand, edge1_supp[i].cmp_head_id, edge2_supp[j].cmp_tail_id});
            i++;
            j++;
        }
        else if(edge1_supp[i].lr_id < edge2_supp[j].lr_id)
        {
            i++;
        }
        else // if(edge1_supp[i].lr_id > edge2_supp[j].lr_id)
        {
            j++;
        }
    }
    sort(shared_supp.begin(), shared_supp.end(), compare_Edge_Supp);
}

int clean_simple_bubbles(vector<BBG_Node_t> &graph, int max_depth, string logpath)
{
    FILE *fp = file_open_write(logpath.c_str());
    // 
    int nb_removed = 0;
    uint32_t num = graph.size();
    uint32_t i, j;
    for(i = 0; i < num; i++)
    {
        if(graph[i].edges[0].size() < 2 && graph[i].edges[1].size() < 2) // cannot be the source of a bubble
            continue;
        if(graph[i].edges[0].size() == 2) // detect bubble at the outgoing edges
        {
            vector<id_strand2_t> path1, path2;
            float cov1, cov2;
            bool found1, found2;
            map<uint32_t, BBG_Edge_t>::iterator it = graph[i].edges[0].begin();
            BBG_Edge_t *edge_start_1 = &(it->second);
            found1 = bbg_find_simple_path_from_source(graph, i, 0, it, max_depth, path1, cov1);
            it++;
            BBG_Edge_t *edge_start_2 = &(it->second);
            found2 = bbg_find_simple_path_from_source(graph, i, 0, it, max_depth, path2, cov2);
            if(found1 && found2)
            {
                // check end nodes
                if(path1.back().id == path2.back().id && path1.back().strand == path2.back().strand)
                {
                    fprintf(fp, "simple_bubble cov:%.2lf ", cov1);
                    for(j = 0; j < path1.size(); j++) fprintf(fp, "%u:%c ", path1[j].id, "+-"[path1[j].strand]);
                    fprintf(fp, "\n              cov:%.2lf ", cov2);
                    for(j = 0; j < path2.size(); j++) fprintf(fp, "%u:%c ", path2[j].id, "+-"[path2[j].strand]);
                    fprintf(fp, "\n");
                    // 
                    BBG_Edge_t *edge_end_1 = bbg_get_edge(path1[path1.size() - 2].id, path1[path1.size() - 2].strand, path1[path1.size() - 1].id, path1[path1.size() - 1].strand, graph);
                    BBG_Edge_t *edge_end_2 = bbg_get_edge(path2[path2.size() - 2].id, path2[path2.size() - 2].strand, path2[path2.size() - 1].id, path2[path2.size() - 1].strand, graph);
                    // 
                    vector<Edge_Supp_t> shared_supp;
                    // for(uint32_t z = 0; z < edge_start_1->edge_supp.size(); z++)
                    //     fprintf(fp, "edge_start_1 %u:%c %u %u\n", edge_start_1->edge_supp[z].lr_id, "+-"[edge_start_1->edge_supp[z].lr_strand], edge_start_1->edge_supp[z].cmp_head_id, edge_start_1->edge_supp[z].cmp_tail_id);
                    // for(uint32_t z = 0; z < edge_end_1->edge_supp.size(); z++)
                    //     fprintf(fp, "edge_end_1 %u:%c %u %u\n", edge_end_1->edge_supp[z].lr_id, "+-"[edge_end_1->edge_supp[z].lr_strand], edge_end_1->edge_supp[z].cmp_head_id, edge_end_1->edge_supp[z].cmp_tail_id);
                    get_shared_lr_supp(edge_start_1->edge_supp, edge_end_1->edge_supp, shared_supp);
                    // for(uint32_t z = 0; z < shared_supp.size(); z++)
                    //     fprintf(fp, "shared %u:%c %u %u\n", shared_supp[z].lr_id, "+-"[shared_supp[z].lr_strand], shared_supp[z].cmp_head_id, shared_supp[z].cmp_tail_id);
                    // for(uint32_t z = 0; z < edge_start_2->edge_supp.size(); z++)
                    //     fprintf(fp, "edge_start_2 %u:%c %u %u\n", edge_start_2->edge_supp[z].lr_id, "+-"[edge_start_2->edge_supp[z].lr_strand], edge_start_2->edge_supp[z].cmp_head_id, edge_start_2->edge_supp[z].cmp_tail_id);
                    // for(uint32_t z = 0; z < edge_end_2->edge_supp.size(); z++)
                    //     fprintf(fp, "edge_end_2 %u:%c %u %u\n", edge_end_2->edge_supp[z].lr_id, "+-"[edge_end_2->edge_supp[z].lr_strand], edge_end_2->edge_supp[z].cmp_head_id, edge_end_2->edge_supp[z].cmp_tail_id);
                    get_shared_lr_supp(edge_start_2->edge_supp, edge_end_2->edge_supp, shared_supp);
                    // for(uint32_t z = 0; z < shared_supp.size(); z++)
                    //     fprintf(fp, "shared %u:%c %u %u\n", shared_supp[z].lr_id, "+-"[shared_supp[z].lr_strand], shared_supp[z].cmp_head_id, shared_supp[z].cmp_tail_id);
                    fprintf(fp, "       shared cov:%zu\n", shared_supp.size());
                    // 
                    if(path1.size() > path2.size()) // path1 is longer
                    {
                        if(cov1 >= cov2)
                        {
                            // keep path1
                            for(j = 0; j < path2.size() - 1; j++)
                                bbg_remove_edge(path2[j].id, path2[j].strand, path2[j+1].id, path2[j+1].strand, graph);
                            fprintf(fp, "           chosen: path1\n");
                        }
                        else
                        {
                            if(shared_supp.size() > cov2)
                            {
                                // keep shared
                                // remove both paths
                                for(j = 0; j < path1.size() - 1; j++)
                                    bbg_remove_edge(path1[j].id, path1[j].strand, path1[j+1].id, path1[j+1].strand, graph);
                                for(j = 0; j < path2.size() - 1; j++)
                                    bbg_remove_edge(path2[j].id, path2[j].strand, path2[j+1].id, path2[j+1].strand, graph);
                                // add a new edge
                                bbg_add_edge_with_supp(graph, path1.front().id, path1.front().strand, path1.back().id, path1.back().strand, shared_supp);
                                fprintf(fp, "           chosen: shared\n");
                            }
                            else
                            {
                                // keep path2
                                for(j = 0; j < path1.size() - 1; j++)
                                    bbg_remove_edge(path1[j].id, path1[j].strand, path1[j+1].id, path1[j+1].strand, graph);
                                fprintf(fp, "           chosen: path2\n");
                            }
                        }
                    }
                    else // path2 is longer
                    {
                        if(cov2 >= cov1)
                        {
                            // keep path2
                            for(j = 0; j < path1.size() - 1; j++)
                                bbg_remove_edge(path1[j].id, path1[j].strand, path1[j+1].id, path1[j+1].strand, graph);
                            fprintf(fp, "           chosen: path2\n");
                        }
                        else
                        {
                            if(shared_supp.size() > cov1)
                            {
                                // keep shared
                                // remove both paths
                                for(j = 0; j < path1.size() - 1; j++)
                                    bbg_remove_edge(path1[j].id, path1[j].strand, path1[j+1].id, path1[j+1].strand, graph);
                                for(j = 0; j < path2.size() - 1; j++)
                                    bbg_remove_edge(path2[j].id, path2[j].strand, path2[j+1].id, path2[j+1].strand, graph);
                                // add a new edge
                                bbg_add_edge_with_supp(graph, path1.front().id, path1.front().strand, path1.back().id, path1.back().strand, shared_supp);
                                fprintf(fp, "           chosen: shared\n");
                            }
                            else
                            {
                                // keep path1
                                for(j = 0; j < path2.size() - 1; j++)
                                    bbg_remove_edge(path2[j].id, path2[j].strand, path2[j+1].id, path2[j+1].strand, graph);
                                fprintf(fp, "           chosen: path1\n");
                            }
                        }
                    }
                    nb_removed++;
                    // 
                    i--;
                    continue;
                }
            }
        }
        if(graph[i].edges[1].size() == 2) // detect bubble at the incoming edges
        {
            vector<id_strand2_t> path1, path2;
            float cov1, cov2;
            bool found1, found2;
            map<uint32_t, BBG_Edge_t>::iterator it = graph[i].edges[1].begin();
            BBG_Edge_t *edge_start_1 = &(it->second);
            found1 = bbg_find_simple_path_from_source(graph, i, 1, it, max_depth, path1, cov1);
            it++;
            BBG_Edge_t *edge_start_2 = &(it->second);
            found2 = bbg_find_simple_path_from_source(graph, i, 1, it, max_depth, path2, cov2);
            if(found1 && found2)
            {
                // check end nodes
                if(path1.back().id == path2.back().id && path1.back().strand == path2.back().strand)
                {
                    fprintf(fp, "simple_bubble cov:%.2lf ", cov1);
                    for(j = 0; j < path1.size(); j++) fprintf(fp, "%u:%c ", path1[j].id, "+-"[path1[j].strand]);
                    fprintf(fp, "\n              cov:%.2lf ", cov2);
                    for(j = 0; j < path2.size(); j++) fprintf(fp, "%u:%c ", path2[j].id, "+-"[path2[j].strand]);
                    fprintf(fp, "\n");
                    // 
                    BBG_Edge_t *edge_end_1 = bbg_get_edge(path1[path1.size() - 2].id, path1[path1.size() - 2].strand, path1[path1.size() - 1].id, path1[path1.size() - 1].strand, graph);
                    BBG_Edge_t *edge_end_2 = bbg_get_edge(path2[path2.size() - 2].id, path2[path2.size() - 2].strand, path2[path2.size() - 1].id, path2[path2.size() - 1].strand, graph);
                    // 
                    vector<Edge_Supp_t> shared_supp;
                    // for(uint32_t z = 0; z < edge_start_1->edge_supp.size(); z++)
                    //     fprintf(fp, "edge_start_1 %u:%c %u %u\n", edge_start_1->edge_supp[z].lr_id, "+-"[edge_start_1->edge_supp[z].lr_strand], edge_start_1->edge_supp[z].cmp_head_id, edge_start_1->edge_supp[z].cmp_tail_id);
                    // for(uint32_t z = 0; z < edge_end_1->edge_supp.size(); z++)
                    //     fprintf(fp, "edge_end_1 %u:%c %u %u\n", edge_end_1->edge_supp[z].lr_id, "+-"[edge_end_1->edge_supp[z].lr_strand], edge_end_1->edge_supp[z].cmp_head_id, edge_end_1->edge_supp[z].cmp_tail_id);
                    get_shared_lr_supp(edge_start_1->edge_supp, edge_end_1->edge_supp, shared_supp);
                    // for(uint32_t z = 0; z < shared_supp.size(); z++)
                    //     fprintf(fp, "shared %u:%c %u %u\n", shared_supp[z].lr_id, "+-"[shared_supp[z].lr_strand], shared_supp[z].cmp_head_id, shared_supp[z].cmp_tail_id);
                    // for(uint32_t z = 0; z < edge_start_2->edge_supp.size(); z++)
                    //     fprintf(fp, "edge_start_2 %u:%c %u %u\n", edge_start_2->edge_supp[z].lr_id, "+-"[edge_start_2->edge_supp[z].lr_strand], edge_start_2->edge_supp[z].cmp_head_id, edge_start_2->edge_supp[z].cmp_tail_id);
                    // for(uint32_t z = 0; z < edge_end_2->edge_supp.size(); z++)
                    //     fprintf(fp, "edge_end_2 %u:%c %u %u\n", edge_end_2->edge_supp[z].lr_id, "+-"[edge_end_2->edge_supp[z].lr_strand], edge_end_2->edge_supp[z].cmp_head_id, edge_end_2->edge_supp[z].cmp_tail_id);
                    get_shared_lr_supp(edge_start_2->edge_supp, edge_end_2->edge_supp, shared_supp);
                    // for(uint32_t z = 0; z < shared_supp.size(); z++)
                    //     fprintf(fp, "shared %u:%c %u %u\n", shared_supp[z].lr_id, "+-"[shared_supp[z].lr_strand], shared_supp[z].cmp_head_id, shared_supp[z].cmp_tail_id);
                    fprintf(fp, "       shared cov:%zu\n", shared_supp.size());
                    //  
                    if(path1.size() > path2.size()) // path1 is longer
                    {
                        if(cov1 >= cov2)
                        {
                            // keep path1
                            for(j = 0; j < path2.size() - 1; j++)
                                bbg_remove_edge(path2[j].id, path2[j].strand, path2[j+1].id, path2[j+1].strand, graph);
                            fprintf(fp, "           chosen: path1\n");
                        }
                        else
                        {
                            if(shared_supp.size() > cov2)
                            {
                                // keep shared
                                // remove both paths
                                for(j = 0; j < path1.size() - 1; j++)
                                    bbg_remove_edge(path1[j].id, path1[j].strand, path1[j+1].id, path1[j+1].strand, graph);
                                for(j = 0; j < path2.size() - 1; j++)
                                    bbg_remove_edge(path2[j].id, path2[j].strand, path2[j+1].id, path2[j+1].strand, graph);
                                // add a new edge
                                bbg_add_edge_with_supp(graph, path1.front().id, path1.front().strand, path1.back().id, path1.back().strand, shared_supp);
                                fprintf(fp, "           chosen: shared\n");
                            }
                            else
                            {
                                // keep path2
                                for(j = 0; j < path1.size() - 1; j++)
                                    bbg_remove_edge(path1[j].id, path1[j].strand, path1[j+1].id, path1[j+1].strand, graph);
                                fprintf(fp, "           chosen: path2\n");
                            }
                        }
                    }
                    else // path2 is longer
                    {
                        if(cov2 >= cov1)
                        {
                            // keep path2
                            for(j = 0; j < path1.size() - 1; j++)
                                bbg_remove_edge(path1[j].id, path1[j].strand, path1[j+1].id, path1[j+1].strand, graph);
                            fprintf(fp, "           chosen: path2\n");
                        }
                        else
                        {
                            if(shared_supp.size() > cov1)
                            {
                                // keep shared
                                // remove both paths
                                for(j = 0; j < path1.size() - 1; j++)
                                    bbg_remove_edge(path1[j].id, path1[j].strand, path1[j+1].id, path1[j+1].strand, graph);
                                for(j = 0; j < path2.size() - 1; j++)
                                    bbg_remove_edge(path2[j].id, path2[j].strand, path2[j+1].id, path2[j+1].strand, graph);
                                // add a new edge
                                bbg_add_edge_with_supp(graph, path1.front().id, path1.front().strand, path1.back().id, path1.back().strand, shared_supp);
                                fprintf(fp, "           chosen: shared\n");
                            }
                            else
                            {
                                // keep path1
                                for(j = 0; j < path2.size() - 1; j++)
                                    bbg_remove_edge(path2[j].id, path2[j].strand, path2[j+1].id, path2[j+1].strand, graph);
                                fprintf(fp, "           chosen: path1\n");
                            }
                        }
                    }
                    nb_removed++;
                    // 
                    i--;
                    continue;
                }
            }
        }
    }
    fclose(fp);
    return nb_removed;
}

// starts at the potential source node, finds the sink node of the super-bubble and returns the best supported path as well as all edges involved in the super bubble
// source node is a branching node
// TODO: implement max_dist
bool detect_super_bubble(vector<BBG_Node_t> &graph, int max_dist, uint32_t src_node, uint32_t src_rev, vector<uint32_t> &best_path, set<pair<uint32_t, uint32_t>> &bubble_edges)
{
    stack<uint32_t> S; // vertices with all "incoming" edges visited
    S.push((src_node << 1) | src_rev);
    unordered_map<uint32_t, int32_t> visited;
    unordered_map<uint32_t, int32_t> gamma; // gamma letter in pseudocode. Number of unvisited incoming edges
    unordered_map<uint32_t, vector<uint32_t>> path; // path with the best support up to this vertext
    unordered_map<uint32_t, uint32_t> support; // sum of edge_supp for the path up to this vertex
    visited[S.top()] = 1;
    path[S.top()].push_back(S.top());
    support[S.top()] = 0;
    int p = 0; // number of visited vertices never added to S
    // fprintf(stdout, "\n");
    while(S.empty() == false)
    {
        uint32_t vertex_v = S.top();
        uint32_t curr_node = vertex_v >> 1;
        uint32_t curr_rev  = vertex_v & 1;
        S.pop();
        // fprintf(stdout, "popped\t%u:%c\tsize:%zu\n", curr_node, "+-"[curr_rev], S.size());
        for(map<uint32_t, BBG_Edge_t>::iterator it = graph[curr_node].edges[curr_rev].begin(); it != graph[curr_node].edges[curr_rev].end(); it++)
        {
            bubble_edges.insert({vertex_v, it->first});
            uint32_t next_node = it->first >> 1;
            uint32_t next_rev = it->first & 1;
            uint32_t next_supp = it->second.edge_supp.size();
            uint32_t vertex_w = (next_node << 1) | next_rev;
            // fprintf(stdout, "\tedge %u:%c -> %u:%c\tcov:%zu\n", curr_node, "+-"[curr_rev], next_node, "+-"[next_rev], it->second.edge_supp.size());
            if(next_node == curr_node) // a circle involving the starting node
            {
                return false;
            }
            if(visited.count(vertex_w) == 0) // not visited before
            {
                gamma[vertex_w] = graph[next_node].edges[1-next_rev].size();
                visited[vertex_w] = 1;
                p++;
            }
            if(support.count(vertex_w) == 0 || double(support[vertex_v] + next_supp) / path[vertex_v].size() > double(support[vertex_w]) / (path[vertex_v].size() - 1))
            {
                support[vertex_w] = support[vertex_v] + next_supp;
                path[vertex_w] = path[vertex_v];
                path[vertex_w].push_back(vertex_w);
            }
            gamma[vertex_w]--;
            if(gamma[vertex_w] == 0)
            {
                if(graph[next_node].edges[next_rev].size() > 0)
                {
                    S.push(vertex_w);
                    p--;
                    // fprintf(stdout, "pushed\t%u:%c\tsize:%zu\n", next_node, "+-"[next_rev], S.size());
                }
            }
        }
        if(S.size() == 1 && p == 0)
        {
            best_path = path[S.top()];
            // uint32_t vertex_sink = S.top();
            // fprintf(stdout, "\tsink_found %u:%c\tcov:%zu:%u:%.2lf\n", vertex_sink >> 1, "+-"[vertex_sink & 1], path[vertex_sink].size(), support[vertex_sink], double(support[vertex_sink]) / (path[vertex_sink].size() - 1));
            // fprintf(stdout, "\tedges:\n");
            // for(set<pair<uint32_t, uint32_t>>::iterator it = bubble_edges.begin(); it != bubble_edges.end(); it++)
            // {
            //     fprintf(stdout, "%u:%c -> %u:%c\n", it->first >> 1, "+-"[it->first & 1], it->second >> 1, "+-"[it->second & 1]);
            // }
            // for(uint32_t j = 0; j < path[vertex_sink].size(); j++)
            // {
            //     fprintf(stdout, "%u:%c ", path[vertex_sink][j] >> 1, "+-"[path[vertex_sink][j] & 1]);
            // }
            // fprintf(stdout, "\n");
            return true;
        }
    }
    return false;
}

// following bubble detection algorithm of miniasm (Algorithm 6 in the paper)
int clean_super_bubbles(vector<BBG_Node_t> &graph, int max_dist, string logpath)
{
    FILE *fp = file_open_write(logpath.c_str());
    // 
    int nb_removed = 0;
    uint32_t num = graph.size();
    uint32_t i, j;
    for(i = 0; i < num; i++)
    {
        // fprintf(stdout, "check:%u\tincoming:%zu\toutgoing:%zu\n", i, graph[i].edges[1].size(), graph[i].edges[0].size());
        if(graph[i].edges[0].size() < 2 && graph[i].edges[1].size() < 2) // cannot be the source of a bubble
            continue;
        if(graph[i].edges[0].size() >= 2) // detect bubble at the outgoing edges
        {
            vector<uint32_t> best_path;
            set<pair<uint32_t, uint32_t>> bubble_edges;
            if(detect_super_bubble(graph, max_dist, i, 0, best_path, bubble_edges))
            {
                fprintf(fp, "bubble_src %u:%c\tbubble_sink %u:%c\n", i, '+', best_path.back() >> 1, "+-"[best_path.back() & 1]);
                fprintf(fp, "\tbest_path ");
                for(j = 0; j < best_path.size(); j++)
                {
                    fprintf(fp, "%u:%c ", best_path[j] >> 1, "+-"[best_path[j] & 1]);
                }
                fprintf(fp, "\n");
                for(j = 0; j < best_path.size() - 1; j++)
                {
                    bubble_edges.erase({best_path[j], best_path[j+1]});
                }
                fprintf(fp, "\tremoved_edges:\n");
                for(set<pair<uint32_t, uint32_t>>::iterator it = bubble_edges.begin(); it != bubble_edges.end(); it++)
                {
                    bbg_remove_edge(it->first >> 1, it->first & 1, it->second >> 1, it->second & 1, graph);
                    fprintf(fp, "\t\t%u:%c -> %u:%c\n", it->first >> 1, "+-"[it->first & 1], it->second >> 1, "+-"[it->second & 1]);
                }
                fprintf(fp, "\n");
                nb_removed++;
                // 
                i--;
                continue;
            }
            // else
            // {
            //     fprintf(stdout, "\tsink NOT FOUND\n");
            // }
        }
        if(graph[i].edges[1].size() >= 2) // detect bubble at the incoming edges
        {
            vector<uint32_t> best_path;
            set<pair<uint32_t, uint32_t>> bubble_edges;
            if(detect_super_bubble(graph, max_dist, i, 1, best_path, bubble_edges))
            {
                fprintf(fp, "bubble_src %u:%c\tbubble_sink %u:%c\n", i, '-', best_path.back() >> 1, "+-"[best_path.back() & 1]);
                fprintf(fp, "\tbest_path ");
                for(j = 0; j < best_path.size(); j++)
                {
                    fprintf(fp, "%u:%c ", best_path[j] >> 1, "+-"[best_path[j] & 1]);
                }
                fprintf(fp, "\n");
                for(j = 0; j < best_path.size() - 1; j++)
                {
                    bubble_edges.erase({best_path[j], best_path[j+1]});
                }
                fprintf(fp, "\tremoved_edges:\n");
                for(set<pair<uint32_t, uint32_t>>::iterator it = bubble_edges.begin(); it != bubble_edges.end(); it++)
                {
                    bbg_remove_edge(it->first >> 1, it->first & 1, it->second >> 1, it->second & 1, graph);
                    fprintf(fp, "\t\t%u:%c -> %u:%c\n", it->first >> 1, "+-"[it->first & 1], it->second >> 1, "+-"[it->second & 1]);
                }
                fprintf(fp, "\n");
                nb_removed++;
                // 
                i--;
                continue;
            }
            // else
            // {
            //     fprintf(stdout, "\tsink NOT FOUND\n");
            // }
        }
    }
    fclose(fp);
    return nb_removed;
}

// int remove_tips_old(vector<BBG_Node_t> &graph, int max_depth, string logpath)
// {
//     FILE *fp;
//     if(max_depth == 1)
//         fp = file_open_write(logpath.c_str());
//     else
//         fp = file_open_append(logpath.c_str());
//     // 
//     int nb_removed = 0;
//     uint32_t num = graph.size();
//     uint32_t i, j;
//     for(i = 0; i < num; i++)
//     {
//         uint32_t src_node, src_strand;
//         uint32_t curr_node, curr_strand;
//         if(graph[i].edges[1].size() == 0 && graph[i].edges[0].size() == 1) // no incoming one outgoing
//         {
//             src_node = i;
//             src_strand = 0;
//         }
//         else if(graph[i].edges[1].size() == 1 && graph[i].edges[0].size() == 0) // one incoming no outgoing
//         {
//             src_node = i;
//             src_strand = 1;
//         }
//         else continue;
//         // 
//         curr_node = src_node;
//         curr_strand = src_strand;
//         // fprintf(stdout, "tip: %u:%c\n", curr_node, "+-"[curr_strand]);
//         vector<tuple<uint32_t, uint32_t, uint32_t, uint32_t>> edges_to_remove;
//         int depth = 0;
//         do
//         {
//             // fprintf(stdout, "\tadd: %u:%c\n", curr_node, "+-"[curr_strand]);
//             if(++depth > max_depth) break;
//             if(graph[curr_node].edges[curr_strand].size() == 0) break;
//             map<uint32_t, BBG_Edge_t>::iterator it = graph[curr_node].edges[curr_strand].begin();
//             edges_to_remove.push_back(make_tuple(curr_node, curr_strand, it->first >> 1, it->first & 1));
//             curr_node = it->first >> 1;
//             curr_strand = it->first & 1;
//         } while (graph[curr_node].edges[curr_strand].size() <= 1 && graph[curr_node].edges[1 - curr_strand].size() <= 1);
//         // 
//         if(depth > max_depth) continue;
//         if(graph[curr_node].edges[curr_strand].size() == 0) continue;
//         // 
//         fprintf(fp, "tip_len:%zu\t%u:%c -> %u:%c\n", edges_to_remove.size(), get<0>(edges_to_remove.front()), "+-"[get<1>(edges_to_remove.front())], get<2>(edges_to_remove.back()), "+-"[get<3>(edges_to_remove.back())]);
//         for(j = 0; j < edges_to_remove.size(); j++)
//         {
//             bbg_remove_edge(get<0>(edges_to_remove[j]), get<1>(edges_to_remove[j]), get<2>(edges_to_remove[j]), get<3>(edges_to_remove[j]), graph);
//         }
//         nb_removed++;
//     }
//     fclose(fp);
//     return nb_removed;
// }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// uint32_t cl_estimate_path_length(vector<id_strand2_t> &path, vector<cns_info_t> &cns_list)
// {
//     int32_t num = cns_list.size();
//     uint32_t len = 0;
//     int32_t i;
//     // add estimated length of consensus between anchors
//     for(i = 0; i < num; i++)
//     {
//         len += cns_list[i].len_estimate;
//     }
//     // add the length of shared region of short read contigs
//     for(i = 0; i < num - 1; i++)
//     {
//         if(path[i+1].strand == 0)
//             len += cns_list[i+1].prev_end - cns_list[i].next_beg + 1;
//         else
//             len += cns_list[i].next_beg - cns_list[i+1].prev_end + 1;
//     }
//     return len;
// }
