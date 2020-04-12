#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdint.h>
#include <zlib.h>
#include <string>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

using namespace std;

template <typename T>
T str2type(string str)
{
    T n;
    istringstream sin(str);
    sin >> n;
    return n;
}

template <typename T>
string type2str(T v)
{
    ostringstream sout;
    sout << v;
    return sout.str();
}

int main(int argc, char* argv[])
{
    if(argc == 2 && strcmp(argv[1], "-h") == 0)
    {
        fprintf(stderr, "usage: ./nooverlap unitigs.fa kmerSize\n");
        return EXIT_SUCCESS;
    }
    if(argc < 3)
    {
        fprintf(stderr, "usage: ./nooverlap unitigs.fa kmerSize\n");
        return EXIT_FAILURE;
    }
    //
    int overlapLen = str2type<int>(argv[2]) - 1;
    //
    gzFile fp = gzopen(argv[1], "r");
    if(fp == NULL)
    {
        fprintf(stderr, "[ERROR] could not open file: %s\n", argv[1]);
        return EXIT_FAILURE;
    }
    kseq_t *seq = kseq_init(fp);
    while(kseq_read(seq) >= 0)
    {
        fprintf(stdout, ">%s %s\n", seq->name.s, seq->comment.s);
        string dummy, link;
        istringstream sin(seq->comment.s);
        sin >> dummy;
        sin >> dummy;
        sin >> dummy;
        bool incoming = false;
        bool outgoing = false;
        while(sin >> link)
        {
            char from_type = link[2];
            // char to_type = link[link.size()-1];

            if(from_type == '+')
                outgoing = true;
            else if(from_type == '-')
                incoming = true;
        }
        //
        string node = seq->seq.s;
        // 
        if(incoming == true)
        {
            node = node.substr(overlapLen / 2);
        }
        if(outgoing == true)
        {
            node = node.substr(0, node.size() - (overlapLen / 2));
        }
        fprintf(stdout, "%s\n", node.c_str());
    }
    kseq_destroy(seq);
    gzclose(fp);
    // 
    return EXIT_SUCCESS;
}
