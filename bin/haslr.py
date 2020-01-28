#!/usr/bin/env python3

haslr_version = '0.8a1'

import os
import sys
import argparse
import multiprocessing
import subprocess
import glob

def main():
    args = parse_options()
    
    path_exe = os.path.abspath(__file__) # should I use realpath instead?
    dir_exe  = os.path.dirname(path_exe)
    path_haslr_assemble = dir_exe + '/haslr_assemble'
    check_prog(path_haslr_assemble)
    path_nooverlap = dir_exe + '/minia_nooverlap'
    check_prog(path_nooverlap)
    path_fastutils = dir_exe + '/fastutils'
    check_prog(path_fastutils)
    path_minia = dir_exe + '/minia'
    check_prog(path_minia)
    path_minimap2 = dir_exe + '/minimap2'
    check_prog(path_minimap2)
    
    sys.stdout.write('number of threads: {}\n'.format(args.threads))
    sys.stdout.write('output directory: {}\n'.format(args.out))
    os.makedirs(args.out, exist_ok=True)
    # 
    lr_name = args.out + '/lr' + str(args.cov_lr) + 'x.fasta'
    sys.stdout.write('subsampling {0}x long reads to {1}... '.format(args.cov_lr, lr_name))
    sys.stdout.flush()
    if os.path.isfile(lr_name) == False:
        with open(lr_name, 'w') as fp:
            try:
                completed = subprocess.run([path_fastutils, 'subsample', '-i', args.long, '-d', str(args.cov_lr), '-g', args.genome, '-lnk'], stdout=fp)
            except subprocess.CalledProcessError as err:
                sys.stdout.write('failed\nERROR: {}\n'.format(err))
                sys.exit(os.EX_SOFTWARE)
        sys.stdout.write('done\n')
        sys.stdout.flush()
    else:
        sys.stdout.write('already exists\n')
        sys.stdout.flush()
    # 
    minimum_sr_contig = 250 # change if necessary
    sys.stdout.write('assembling short reads using Minia... '.format(args.minia_kmer, args.minia_solid))
    sys.stdout.flush()
    sr_name = args.out + '/sr.fofn'
    sr_asm_prefix = args.out + '/sr_k{}_a{}'.format(args.minia_kmer, args.minia_solid)
    sr_asm_name = sr_asm_prefix + '.' + args.minia_asm + '.fa'
    if os.path.isfile(sr_asm_name) == False:
        with open(sr_name, 'w') as fp:
            for fn in args.short:
                fp.write(fn + '\n')
        with open(sr_asm_prefix + '.log', 'w') as fp:
            try:
                completed = subprocess.run([path_minia, '-nb-cores', str(args.threads), '-out-dir', args.out, '-out-tmp', args.out, '-out', sr_asm_prefix, '-in', sr_name, '-kmer-size', str(args.minia_kmer), '-abundance-min', str(args.minia_solid), '-no-ec-removal'], stdout=fp, stderr=fp)
            except subprocess.CalledProcessError as err:
                sys.stdout.write('failed\nERROR: {}\n'.format(err))
                sys.exit(os.EX_SOFTWARE)
        sys.stdout.write('done\n')
        sys.stdout.flush()
    else:
        sys.stdout.write('already exists\n')
        sys.stdout.flush()
    # 
    fileList = glob.glob(sr_asm_prefix + '.unitigs.fa.glue*')
    for fn in fileList:
        try:
            os.remove(fn)
        except:
            sys.stdout.write('ERROR: cannot delete file: {}\n'.format(fn))
    # 
    sys.stdout.write('removing overlaps in short read assembly... ')
    sys.stdout.flush()
    sr_asm_name_noov = sr_asm_prefix + '.' + args.minia_asm + '.nooverlap.fa'
    if os.path.isfile(sr_asm_name_noov) == False:
        with open(sr_asm_name_noov, 'w') as fp:
            try:
                completed = subprocess.run([path_nooverlap, sr_asm_name, str(args.minia_kmer)], stdout=fp, stderr=subprocess.DEVNULL)
            except subprocess.CalledProcessError as err:
                sys.stdout.write('failed\nERROR: {}\n'.format(err))
                sys.exit(os.EX_SOFTWARE)
        sys.stdout.write('done\n')
    else:
        sys.stdout.write('already exists\n')
    sys.stdout.write('removing short sequences in short read assembly... ')
    sys.stdout.flush()
    sr_asm_name_good = '{0}.{1}.nooverlap.{2}.fa'.format(sr_asm_prefix, args.minia_asm, minimum_sr_contig)
    if os.path.isfile(sr_asm_name_good) == False:
        with open(sr_asm_name_good, 'w') as fp:
            try:
                completed = subprocess.run([path_fastutils, 'format', '-i', sr_asm_name_noov, '-m', str(minimum_sr_contig), '-c'], stdout=fp, stderr=subprocess.DEVNULL)
            except subprocess.CalledProcessError as err:
                sys.stdout.write('failed\nERROR: {}\n'.format(err))
                sys.exit(os.EX_SOFTWARE)
        sys.stdout.write('done\n')
        sys.stdout.flush()
    else:
        sys.stdout.write('already exists\n')
        sys.stdout.flush()
    # 
    sys.stdout.write('aligning long reads to short read assembly using minimap2... ')
    sys.stdout.flush()
    map_lr2contig = '{0}/map_{1}_k{2}_a{3}_lr{4}x'.format(args.out, args.minia_asm, args.minia_kmer, args.minia_solid, args.cov_lr)
    if args.type in ['corrected', 'ccs']:
        map_opt = '-k19'
    elif args.type == 'pacbio':
        map_opt = '-Hk17'
    elif args.type == 'nanopore':
        map_opt = '-k15'
    if os.path.isfile(map_lr2contig + '.paf') == False:
        with open(map_lr2contig + '.paf', 'w') as fpout, open(map_lr2contig + '.log', 'w') as fperr:
            try:
                completed = subprocess.run([path_minimap2, '-t', str(args.threads), '--secondary=no', '-c', map_opt, sr_asm_name_good, lr_name], stdout=fpout, stderr=fperr)
            except subprocess.CalledProcessError as err:
                sys.stdout.write('failed\nERROR: {}\n'.format(err))
                sys.exit(os.EX_SOFTWARE)
        sys.stdout.write('done\n')
        sys.stdout.flush()
    else:
        sys.stdout.write('already exists\n')
        sys.stdout.flush()
    # 
    sys.stdout.write('assembling long reads using HASLR... ')
    sys.stdout.flush()
    lr_asm_dir = '{0}/asm_{1}_k{2}_a{3}_lr{4}x_b{5}_s{6}_sim{7}'.format(args.out, args.minia_asm, args.minia_kmer, args.minia_solid, args.cov_lr, args.aln_block, args.edge_sup, args.aln_sim)
    if os.path.isfile(lr_asm_dir + '/asm.final.fa') == False:
        with open(lr_asm_dir + '.out', 'w') as fpout, open(lr_asm_dir + '.err', 'w') as fperr:
            try:
                completed = subprocess.run([path_haslr_assemble, '-t', str(args.threads), '-c', sr_asm_name_noov, '-l', lr_name, '-m', map_lr2contig + '.paf', '-d', lr_asm_dir, '--aln-block', str(args.aln_block), '--aln-sim', str(args.aln_sim), '--edge-sup', str(args.edge_sup)], stdout=fpout, stderr=fperr)
            except subprocess.CalledProcessError as err:
                sys.stdout.write('failed\nERROR: {}\n'.format(err))
                sys.exit(os.EX_SOFTWARE)
        sys.stdout.write('done\n')
        sys.stdout.flush()
    else:
        sys.stdout.write('already exists\n')
        sys.stdout.flush()
    # print(lrname)
    # print(os.path.isfile(lrname))

    # print(args.out)
    # print(args.genome)
    # print(args.long)
    # print(args.type)
    # print(args.short)
    # print(type(args.short))

    # print(args.threads)
    # print(args.cov_lr)
    # print(args.aln_block)
    # print(args.edge_sup)
    # print(args.minia_kmer)
    # print(args.minia_solid)

    sys.exit(os.EX_OK)

def check_prog(prog):
    sys.stdout.write('checking {}: '.format(prog))
    sys.stdout.flush()
    if(os.path.isfile(prog) == False):
        sys.stdout.write('not found\n')
        sys.exit(os.EX_SOFTWARE)
    try:
        completed = subprocess.run([prog, '-h'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if completed.returncode == 0:
            sys.stdout.write('ok\n')
        else:
            sys.stdout.write('failed\n')
            sys.exit(os.EX_SOFTWARE)
    except subprocess.CalledProcessError as err:
        sys.stdout.write('failed\nERROR: {}\n'.format(err))
        sys.exit(os.EX_SOFTWARE)

class CustomHelpFormatter(argparse.HelpFormatter):
    def _format_action_invocation(self, action):
        self._max_help_position = 50
        self._width = 1000
        if not action.option_strings or action.nargs == 0:
            return super()._format_action_invocation(action)
        default = self._get_default_metavar_for_optional(action)
        args_string = self._format_args(action, default)
        return ', '.join(action.option_strings) + ' ' + args_string

def parse_options():
    fmt = lambda prog: CustomHelpFormatter(prog)
    parser = argparse.ArgumentParser(add_help=False, formatter_class=fmt, usage='haslr.py [-t THREADS] -o OUT_DIR -g GENOME_SIZE -l LONG -x LONG_TYPE -s SHORT [SHORT ...]')
    # required
    parser_req = parser.add_argument_group(title='required arguments')
    parser_req.add_argument('-o', '--out', type=str, help='output directory', required=True, metavar='OUT_DIR')
    parser_req.add_argument('-g', '--genome', type=str, help='estimated genome size; accepted suffixes are k,m,g', required=True, metavar='GENOME_SIZE')
    parser_req.add_argument('-l', '--long', type=str, help='long read file', required=True)
    # parser_req.add_argument('-x', '--type', type=str, help='type of long reads chosen from {pacbio,nanopore, corrected, ccs}', required=True, choices=['pacbio', 'nanopore', 'corrected', 'ccs'], metavar='LONG_TYPE')
    parser_req.add_argument('-x', '--type', type=str, help='type of long reads chosen from {pacbio,nanopore}', required=True, choices=['pacbio', 'nanopore'], metavar='LONG_TYPE')
    parser_req.add_argument('-s', '--short', type=str, help='short read file', nargs='+')
    # parser_req.add_argument('-1', '--short1', type=str, help='forward file of a short read dataset', required=True, nargs='+')
    # parser_req.add_argument('-2', '--short2', type=str, help='backward file of a short read dataset', required=True, nargs='+')
    # optional
    parser_opt = parser.add_argument_group(title='optional arguments')
    parser_opt.add_argument('-t', '--threads', type=int, help='number of CPU threads to use [1]', default=1)
    parser_opt.add_argument('--cov-lr', type=int, help='amount of long read coverage to use for assembly [25]', default=25)
    parser_opt.add_argument('--aln-block', type=int, help='minimum length of alignment block [500]', default=500)
    parser_opt.add_argument('--aln-sim', type=float, help='minimum alignment similarity [0.85]', default=0.85)
    parser_opt.add_argument('--edge-sup', type=int, help='minimum number of long read supporting each edge [3]', default=3)
    parser_opt.add_argument('--minia-kmer', type=int, help='kmer size used by minia [49]', default=49)
    parser_opt.add_argument('--minia-solid', type=int, help='minimum kmer abundance used by minia [3]', default=3)
    parser_opt.add_argument('--minia-asm', type=str, help='type of minia assembly chosen from {contigs,unitigs} [contigs]', metavar='MINIA_ASM', default='contigs', choices=['contigs', 'unitigs'])
    parser_opt.add_argument('-v', '--version', action='version', help='print version', version=haslr_version)
    parser_opt.add_argument('-h', '--help', action='help', help='show this help message and exit')
    # 
    if len(sys.argv) == 1:
        parser.print_usage()
        sys.exit(os.EX_USAGE)
    args = parser.parse_args()
    # 
    if args.type != 'ccs' and args.short == None:
        sys.stdout.write('{0}: error: argument -s/--short is required when --type is "{1}"\n'.format(os.path.basename(__file__), args.type))
        sys.exit(os.EX_USAGE)
    # 
    args.threads = max(args.threads, 1)
    args.threads = min(args.threads, multiprocessing.cpu_count())
    # 
    if os.path.isfile(args.long) == False:
        sys.stdout.write('{0}: error: could not find file {1}\n'.format(os.path.basename(__file__), args.long))
        sys.exit(os.EX_USAGE)
    if args.short != None:
        for fn in args.short:
            if os.path.isfile(fn) == False:
                sys.stdout.write('{0}: error: could not find file {1}\n'.format(os.path.basename(__file__), fn))
                sys.exit(os.EX_USAGE)

    # print(args.short + [args.long])
    # 
    args.out = os.path.abspath(args.out)
    args.long = os.path.abspath(args.long)
    for i, f in enumerate(args.short):
        args.short[i] = os.path.abspath(f)

    return args

if __name__ == '__main__':
    main()
