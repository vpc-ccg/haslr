/*****************************************************
 * Author: Ehsan Haghshenas (ehaghshe AT sfu DOT ca) *
 *****************************************************/

#include "Common.hpp"
#include <zlib.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

char tableRev[128] = {
    'N','N','N','N',   'N','N','N','N',   'N','N','N','N',   'N','N','N','N', 
    'N','N','N','N',   'N','N','N','N',   'N','N','N','N',   'N','N','N','N', 
    'N','N','N','N',   'N','N','N','N',   'N','N','N','N',   'N','N','N','N', 
    'N','N','N','N',   'N','N','N','N',   'N','N','N','N',   'N','N','N','N', 
    'N','T','N','G',   'N','N','N','C',   'N','N','N','N',   'N','N','N','N', 
    'N','N','N','N',   'A','N','N','N',   'N','N','N','N',   'N','N','N','N', 
    'N','T','N','G',   'N','N','N','C',   'N','N','N','N',   'N','N','N','N', 
    'N','N','N','N',   'A','N','N','N',   'N','N','N','N',   'N','N','N','N'
};

global_options_t gopt;

// void file_clean(string path)
// {
//     FILE *fp = fopen(path.c_str(), "w");
//     if(fp == NULL)
//     {
//         fprintf(stderr, "[ERROR] (Common::file_clean) could not open file: %s\n", path.c_str());
//         exit(EXIT_FAILURE);
//     }
//     fclose(fp);
// }

FILE* file_open_write(string path)
{
    FILE *fp = fopen(path.c_str(), "w");
    if(fp == NULL)
    {
        fprintf(stderr, "[ERROR] (Common::file_open_write) could not open file: %s\n", path.c_str());
        exit(EXIT_FAILURE);
    }
    return fp;
}

FILE* file_open_append(string path)
{
    FILE *fp = fopen(path.c_str(), "a");
    if(fp == NULL)
    {
        fprintf(stderr, "[ERROR] (Common::file_open_append) could not open file: %s\n", path.c_str());
        exit(EXIT_FAILURE);
    }
    return fp;
}

bool file_exists(string path)
{
    FILE *fp = fopen(path.c_str(), "rb");
    if(fp == NULL)
    {
        return false;
    }
    else
    {
        fclose(fp);
        return true;
    }
}

void str_split(string str, char delim, vector<string> &v)
{
    v.clear();
    size_t p1 = 0;
    size_t p2 = 0;
    while((p2 = str.find(delim, p1)) != string::npos)
    {
        v.push_back(str.substr(p1, p2-p1));
        p1 = p2+1;
    }
    v.push_back(str.substr(p1));
}

int32_t count_seqs(string path)
{
    int32_t cnt = 0;

    gzFile fp = gzopen(path.c_str(), "r");
    if(fp == NULL)
    {
        fprintf(stderr, "[ERROR] (Common::count_seqs) could not open file: %s\n", path.c_str());
        exit(EXIT_FAILURE);
    }

    kseq_t *seq = kseq_init(fp);
    while(kseq_read(seq) > 0)
    {
        cnt++;
    }
    kseq_destroy(seq);
    gzclose(fp);

    return cnt;
}

string expand_cigar(char *cigar)
{
    string res = "";
    int offset = 0;
    int consumed;
    uint32_t op_len;
    char op;
    while(sscanf(cigar + offset, "%u%c%n", &op_len, &op, &consumed) > 0)
    {
        offset += consumed;
        res += string(op_len, op);
    }
    return res;
}

string collapse_cigar(string cigar_exp)
{
    uint32_t i;
    char last = '?';
    int cnt = 0;
    string res = "";
    for(i = 0; i < cigar_exp.size(); i++)
    {
        if(cigar_exp[i] == last)
        {
            cnt++;
        }
        else
        {
            if(cnt > 0)
            {
                res += type2str<int>(cnt) + last;
            }
            cnt = 1;
            last = cigar_exp[i];
        }
    }
    if(cnt > 0)
    {
        res += type2str<int>(cnt) + last;
    }
    return res;
}

double get_cpu_time()
{
    struct rusage t;
    getrusage(RUSAGE_SELF, &t);
    return t.ru_utime.tv_sec + t.ru_utime.tv_usec / 1000000.0 + t.ru_stime.tv_sec + t.ru_stime.tv_usec / 1000000.0;
}

double get_real_time()
{
    struct timeval t;
    struct timezone tz;
    gettimeofday(&t, &tz);
    return t.tv_sec + t.tv_usec / 1000000.0;
}

// string str2Lower(string s)
// {
//     string sLow = "";
//     for(int i = 0; i < s.size(); i++)
//     {
//         sLow += tolower(s[i]);
//     }
//     return sLow;
// }

string reverseString(string str)
{
    int n = str.size();
    string revStr(n, 'N');
    for(int i = 0; i < n; i++)
        revStr[i] = str[n-i-1];
    return revStr;
}

string reverseComplement(string seq)
{
    int len = seq.size();
    string seq_rev(len, 'N');
    for(int i = 0; i < len; i++)
        seq_rev[i] = tableRev[(int)seq[len-i-1]];
    return seq_rev;
}
