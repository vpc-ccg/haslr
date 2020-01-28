/*****************************************************
 * Author: Ehsan Haghshenas (ehaghshe AT sfu DOT ca) *
 *****************************************************/

#include "Common.hpp"
#include "Commandline.hpp"
#include "Compressed_sequence.hpp"
#include "Contig.hpp"
#include "Longread.hpp"
#include "Backbone_graph.hpp"
// #include "Graph_repeat.hpp"
// #include "Align_LR2path.hpp"
#include "Cleaning.hpp"
#include "Assemble.hpp"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <inttypes.h>

using namespace std;

typedef struct
{
    uint32_t contig;
    uint8_t is_rev;
} contig_dir_pair_t;

int main(int argc, char *argv[])
{
    if(parse_command_line(argc, argv) == false)
        exit(EXIT_FAILURE);

    fprintf(stderr, "[NOTE] number of threads: %d\n\n", gopt.num_threads);
    double cputime_start = get_cpu_time();
    double realtime_start = get_real_time();
    srand(time(NULL));

    Contig_List_t contig_list;
    if(file_exists(gopt.out_dir + "/index.contig"))
    {
        // load contigs from index
        fprintf(stderr, "[NOTE] reading contig index: %s...\n", (gopt.out_dir + "/index.contig").c_str());
        read_contig_index(gopt.out_dir + "/index.contig", contig_list);
    }
    else
    {
        // load contigs from fasta file
        fprintf(stderr, "[NOTE] loading contig sequences...\n");
        initialize_contig(contig_list);
        load_contig_compressed(gopt.contig_path, contig_list);
        write_contig_index(gopt.out_dir + "/index.contig", contig_list);
    }
    fprintf(stderr, "       loaded %lu contigs\n", contig_list.contigs_size);
    fprintf(stderr, "       elapsed time %.2lf CPU seconds (%.2lf real seconds)\n\n", get_cpu_time() - cputime_start, get_real_time() - realtime_start);
    // check if contigs are stored correctly
    // print_loaded_contigs(contig_list);
    // return 0;
    // 
    fprintf(stderr, "[NOTE] calculating kmer frequency of unique contigs\n");
    calc_uniq_freq(contig_list);
    fprintf(stderr, "       mean: %.2lf\n", gopt.uniq_freq);
    fprintf(stderr, "       elapsed time %.2lf CPU seconds (%.2lf real seconds)\n\n", get_cpu_time() - cputime_start, get_real_time() - realtime_start);
    
    Longread_List_t lr_list;
    if(file_exists(gopt.out_dir + "/index.longread"))
    {
        // load long reads and alignments from index
        fprintf(stderr, "[NOTE] reading long read and alignment index: %s...\n", (gopt.out_dir + "/index.longread").c_str());
        read_longread_index(gopt.out_dir + "/index.longread", lr_list);
        fprintf(stderr, "       loaded %lu long reads\n", lr_list.reads_size);
        fprintf(stderr, "       loaded %lu alignments\n", lr_list.alignments_size);
    }
    else
    {
        // load long reads from fasta file
        fprintf(stderr, "[NOTE] loading long read sequences...\n");
        initialize_longread(lr_list);
        if(gopt.long_fofn)
        {
            load_longread_compressed_fofn(gopt.long_path, lr_list);
        }
        else
        {
            load_longread_compressed(gopt.long_path, lr_list);
            update_longreads(lr_list);
        }
        fprintf(stderr, "       loaded %lu long reads\n", lr_list.reads_size);
        fprintf(stderr, "       elapsed time %.2lf CPU seconds (%.2lf real seconds)\n\n", get_cpu_time() - cputime_start, get_real_time() - realtime_start);

        // load mapping of long reads onto contigs from PAF file
        fprintf(stderr, "[NOTE] loading alignment between contigs and long reads...\n");
        if(gopt.mapping_fofn)
        {
            load_alignment_fofn(gopt.mapping_path, contig_list, lr_list);
        }
        else
        {
            load_alignment(gopt.mapping_path, contig_list, lr_list);
            update_longreads(lr_list);
        }
        fprintf(stderr, "       loaded %lu alignments\n", lr_list.alignments_size);
        write_longread_index(gopt.out_dir + "/index.longread", lr_list);
    }
    fprintf(stderr, "       elapsed time %.2lf CPU seconds (%.2lf real seconds)\n\n", get_cpu_time() - cputime_start, get_real_time() - realtime_start);
    // check if long reads and alignments are stored correctly
    // print_loaded_lrs(lr_list);
    // print_loaded_alignments(lr_list, gopt.out_dir + "/alignments.original.paf");
    // return 0;

    // // loading long read comments
    // vector<string> comm_list; // list of long read comments
    // comm_list.resize(lr_list.reads_size);
    // load_longread_comments(gopt.long_path, comm_list);

    fprintf(stderr, "[NOTE] fixing overlapping alignments...\n");
    fix_alignments(lr_list);
    // fix_alignments_MT(lr_list);
    fprintf(stderr, "       elapsed time %.2lf CPU seconds (%.2lf real seconds)\n\n", get_cpu_time() - cputime_start, get_real_time() - realtime_start);
    // print_loaded_alignments(lr_list, gopt.out_dir + "/alignments.fixed.paf");

    // build compact long reads
    fprintf(stderr, "[NOTE] building compact long reads...\n");
    vector<vector<Align_Seq_t*>> compact_lr_list;
    build_compact_longreads(lr_list, compact_lr_list, contig_list, gopt.min_aln_block, 1);
    // build_compact_longreads(lr_list, compact_lr_list, contig_list, gopt.min_aln_block, 2);
    print_compact_longreads(compact_lr_list, gopt.out_dir + "/compact_uniq.txt");
    fprintf(stderr, "       elapsed time %.2lf CPU seconds (%.2lf real seconds)\n\n", get_cpu_time() - cputime_start, get_real_time() - realtime_start);

    fprintf(stderr, "[NOTE] building the backbone graph...\n");
    vector<BBG_Node_t> backbone_graph;
    bbg_build_graph(backbone_graph, contig_list, compact_lr_list);
    // bbg_build_graph_2(backbone_graph, contig_list.contigs_size, compact_lr_list);
    bbg_general_stats(backbone_graph, contig_list, gopt.out_dir + "/backbone.01.init.stat");
    bbg_print_graph_gfa(backbone_graph, contig_list, gopt.out_dir + "/backbone.01.init.gfa");
    fprintf(stderr, "       elapsed time %.2lf CPU seconds (%.2lf real seconds)\n\n", get_cpu_time() - cputime_start, get_real_time() - realtime_start);

    fprintf(stderr, "[NOTE] cleaning weak edges...\n");
    int nb_weak_edges = bbg_remove_weak_edges(backbone_graph);
    fprintf(stderr, "       removed %d edges\n", nb_weak_edges);
    bbg_general_stats(backbone_graph, contig_list, gopt.out_dir + "/backbone.02.weakEdge.stat");
    bbg_print_graph_gfa(backbone_graph, contig_list, gopt.out_dir + "/backbone.02.weakEdge.gfa");
    fprintf(stderr, "       elapsed time %.2lf CPU seconds (%.2lf real seconds)\n\n", get_cpu_time() - cputime_start, get_real_time() - realtime_start);

    // two different ways of adding repetitive SRCs to the backbone graph
    // bbg_build_graph_3(backbone_graph, contig_list, compact_lr_list);
    // bbg_add_repetitive_SRCs(backbone_graph, contig_list, compact_lr_list);
    // return 0;

    fprintf(stderr, "[NOTE] cleaning tips...\n");
    int nb_tips = clean_tips(backbone_graph, 1, gopt.out_dir + "/backbone.03.tip.log");
    nb_tips += clean_tips(backbone_graph, 2, gopt.out_dir + "/backbone.03.tip.log");
    nb_tips += clean_tips(backbone_graph, 3, gopt.out_dir + "/backbone.03.tip.log");
    fprintf(stderr, "       removed %d tips\n", nb_tips);
    bbg_general_stats(backbone_graph, contig_list, gopt.out_dir + "/backbone.03.tip.stat");
    bbg_print_graph_gfa(backbone_graph, contig_list, gopt.out_dir + "/backbone.03.tip.gfa");
    fprintf(stderr, "       elapsed time %.2lf CPU seconds (%.2lf real seconds)\n\n", get_cpu_time() - cputime_start, get_real_time() - realtime_start);

    // fprintf(stderr, "[NOTE] cleaning small bubbles...\n");
    // int nb_small_bubbles = clean_small_bubbles(backbone_graph, gopt.out_dir + "/backbone.03.smallbubble.log");
    // fprintf(stderr, "       removed %d small bubbles\n", nb_small_bubbles);
    // bbg_general_stats(backbone_graph, contig_list, gopt.out_dir + "/backbone.03.smallbubble.stat");
    // bbg_print_graph_gfa(backbone_graph, contig_list, gopt.out_dir + "/backbone.03.smallbubble.gfa");
    // fprintf(stderr, "       elapsed time %.2lf CPU seconds (%.2lf real seconds)\n\n", get_cpu_time() - cputime_start, get_real_time() - realtime_start);

    // fprintf(stderr, "[NOTE] cleaning tips...\n");
    // nb_tips = clean_tips(backbone_graph, 1, gopt.out_dir + "/backbone.044.tip.log");
    // nb_tips += clean_tips(backbone_graph, 2, gopt.out_dir + "/backbone.044.tip.log");
    // nb_tips += clean_tips(backbone_graph, 3, gopt.out_dir + "/backbone.044.tip.log");
    // fprintf(stderr, "       removed %d tips\n", nb_tips);
    // bbg_general_stats(backbone_graph, contig_list, gopt.out_dir + "/backbone.044.tip.stat");
    // bbg_print_graph_gfa(backbone_graph, contig_list, gopt.out_dir + "/backbone.044.tip.gfa");
    // fprintf(stderr, "       elapsed time %.2lf CPU seconds (%.2lf real seconds)\n\n", get_cpu_time() - cputime_start, get_real_time() - realtime_start);

    fprintf(stderr, "[NOTE] cleaning simple bubbles...\n");
    int nb_simple_bubbles = clean_simple_bubbles_old(backbone_graph, 4, gopt.out_dir + "/backbone.04.simplebubble.log"); // simple bubbles with at most 3 intermediate nodes
    // int nb_simple_bubbles = clean_simple_bubbles(backbone_graph, 4, gopt.out_dir + "/backbone.04.simplebubble.log"); // simple bubbles with at most 3 intermediate nodes
    fprintf(stderr, "       removed %d simple bubbles\n", nb_simple_bubbles);
    bbg_general_stats(backbone_graph, contig_list, gopt.out_dir + "/backbone.04.simplebubble.stat");
    bbg_print_graph_gfa(backbone_graph, contig_list, gopt.out_dir + "/backbone.04.simplebubble.gfa");
    fprintf(stderr, "       elapsed time %.2lf CPU seconds (%.2lf real seconds)\n\n", get_cpu_time() - cputime_start, get_real_time() - realtime_start);

    // print_untangling_graph_node_list(backbone_graph, contig_list, gopt.out_dir + "/nodes_uniq.txt");

    fprintf(stderr, "[NOTE] cleaning super bubbles...\n");
    int nb_super_bubbles = clean_super_bubbles(backbone_graph, 50000, gopt.out_dir + "/backbone.05.superbubble.log");
    fprintf(stderr, "       removed %d super bubbles\n", nb_super_bubbles);
    bbg_general_stats(backbone_graph, contig_list, gopt.out_dir + "/backbone.05.superbubble.stat");
    bbg_print_graph_gfa(backbone_graph, contig_list, gopt.out_dir + "/backbone.05.superbubble.gfa");
    fprintf(stderr, "       elapsed time %.2lf CPU seconds (%.2lf real seconds)\n\n", get_cpu_time() - cputime_start, get_real_time() - realtime_start);

    fprintf(stderr, "[NOTE] cleaning small bubbles...\n");
    int nb_small_bubbles = clean_small_bubbles(backbone_graph, gopt.out_dir + "/backbone.06.smallbubble.log");
    fprintf(stderr, "       removed %d small bubbles\n", nb_small_bubbles);
    bbg_general_stats(backbone_graph, contig_list, gopt.out_dir + "/backbone.06.smallbubble.stat");
    bbg_print_graph_gfa(backbone_graph, contig_list, gopt.out_dir + "/backbone.06.smallbubble.gfa");
    fprintf(stderr, "       elapsed time %.2lf CPU seconds (%.2lf real seconds)\n\n", get_cpu_time() - cputime_start, get_real_time() - realtime_start);

    // bbg_print_graph_gfa(backbone_graph, contig_list, gopt.out_dir + "/backbone.gfa");
    bbg_report_branching_nodes(backbone_graph, gopt.out_dir + "/backbone.branching.log");
    // return 0;

    fprintf(stderr, "[NOTE] calculating long read coordinates between anchors...\n");
    asm_calc_edge_coordinates_MT(&backbone_graph, &contig_list, &lr_list, &compact_lr_list, gopt.out_dir + "/log_coordinate.txt");
    fprintf(stderr, "       elapsed time %.2lf CPU seconds (%.2lf real seconds)\n\n", get_cpu_time() - cputime_start, get_real_time() - realtime_start);

    fprintf(stderr, "[NOTE] calling consensus sequence between anchors...\n");
    asm_cal_cns_seq_MT(&backbone_graph, &lr_list, gopt.out_dir + "/log_consensus.txt");
    fprintf(stderr, "       elapsed time %.2lf CPU seconds (%.2lf real seconds)\n\n", get_cpu_time() - cputime_start, get_real_time() - realtime_start);

    fprintf(stderr, "[NOTE] generating the assembly from the cleaned backbone graph...\n");
    asm_get_assembly(backbone_graph);
    fprintf(stderr, "       elapsed time %.2lf CPU seconds (%.2lf real seconds)\n\n", get_cpu_time() - cputime_start, get_real_time() - realtime_start);
    
    // free memory
    fprintf(stderr, "[NOTE] cleaning up the memory!\n");
    finalize_contig(contig_list);
    finalize_longread(lr_list);
    fprintf(stderr, "[NOTE] elapsed time %.2lf CPU seconds (%.2lf real seconds)\n\n", get_cpu_time() - cputime_start, get_real_time() - realtime_start);
    fprintf(stderr, "*** BYE ***\n\n");
    return(EXIT_SUCCESS);
}
