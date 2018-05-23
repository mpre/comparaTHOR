#!/usr/bin/env python3

import argparse
import sys
import os
import csv
import glob

ANNOTATION_DIR='/home/prj_gralign/data/BMC-exp-sim/NewAnnos/'
RESULTS_DIR='/home/prj_gralign/data/BMC-exp-sim/Results/'

rmats_spladder_event_map = { 'RI'  : 'intron_retention',
                             'SE'  : 'exon_skip',
                             'A5SS': 'alt_5prime',
                             'A3SS': 'alt_3prime'}

def evaluate_asgal(args, gene_events, exp_events):
    asgal_introns = []
    with open(os.path.join(RESULTS_DIR, args.sample, args.chromosome, 'asgal',
                           args.gene,   args.exp,    args.exp + '.events.csv'),
              newline='') as asgal_csv:
        results = csv.reader(asgal_csv)
        next(results, None) # Drop header
        for row in results:
            asgal_introns += ((int(row[1]) - 1, int(row[2]) + 1),)

    tp=len(set(asgal_introns) & set(exp_events))
    fp=len(set(asgal_introns) - set(exp_events) - set(gene_events.values()))
    fn=len(set(exp_events) - set(asgal_introns))
    return tp, fp, fn

def evaluate_spladder(args, gene_events, exp_events):
    spladder_introns = []
    spladder_files_list = glob.glob(os.path.join(RESULTS_DIR, args.sample, args.chromosome,
                                                 'spladder',  args.gene,   args.exp, '*.txt'))
    for spladder_file in spladder_files_list:
        with open(spladder_file) as spladder_current_event_file:
            next(spladder_current_event_file, None) # Drop header
            for row in spladder_current_event_file:
                _, strand, event_id, _, *positions = row.strip('\n').split('\t')
                if event_id.startswith('mult'):
                    # here spladder outputs multiple positions in columns 2 and 3
                    # We can set them to 0 and avoid problems later on
                    positions[2], positions[3] = 0, 0
                positions = tuple(int(p) for p in positions[:6])
                spladder_introns += (spladder_parse_line(event_id,
                                                         strand,
                                                         positions))
    tp=len(set(spladder_introns) & set(exp_events))
    fp=len(set(spladder_introns) - set(exp_events) - set(gene_events.values()))
    fn=len(set(exp_events) - set(spladder_introns))
    return tp, fp, fn

def evaluate_rmats(args, gene_events, exp_events):
    # We do not consider MXE
    rmats_introns = []
    rmats_files_list = glob.glob(os.path.join(RESULTS_DIR, args.sample, args.chromosome,
                                              'rMATS',     args.gene,   args.exp,
                                              'fromGTF.[!(n,M)]*.txt'))
    for rmats_file in rmats_files_list:
        with open(rmats_file) as rmats_current_event_file:
            next(rmats_current_event_file) # Drop header
            event_type = rmats_file[rmats_file.rfind('.', 0, rmats_file.rfind('.')) + 1:
                                    rmats_file.rfind('.')]
            event_id = rmats_spladder_event_map[event_type]
            for row in rmats_current_event_file:
                _, _, _, _, strand, *positions = row.strip('\n').split('\t')
                positions = tuple(int(p) for p in positions)
                rmats_introns += (rmats_parse_line(event_id,
                                                   strand,
                                                   positions))
    tp=len(set(rmats_introns) & set(exp_events))
    fp=len(set(rmats_introns) - set(exp_events) - set(gene_events.values()))
    fn=len(set(exp_events) - set(rmats_introns))
    return tp, fp, fn

def evaluate_majiq(args, gene_events, exp_events):
    return

def rmats_parse_line(event_id, strand, positions):
    if not event_id.startswith('alt'):
        # Exon skipping, mult exon skipping, intron retention
        return ((positions[3], positions[4] + 1),)
    elif event_id.startswith('alt_5prime'):
        if strand is '+':
            return ((positions[1], positions[4] + 1),
                    (positions[3], positions[4] + 1))
        elif strand is '-':
            return ((positions[5], positions[0] + 1),
                    (positions[5], positions[2] + 1))
        else:
            return ((0, 0), (0, 0))
    elif event_id.startswith('alt_3prime'):
        if strand is '+':
            return ((positions[5], positions[0] + 1),
                    (positions[5], positions[2] + 1))
        elif strand is '-':
            return ((positions[1], positions[4] + 1),
                    (positions[3], positions[4] + 1))
        else:
            return ((0, 0), (0, 0))
    else:
        return (0, 0)

def spladder_parse_line(event_id, strand, positions):
    if not event_id.startswith('alt'):
        # Exon skipping, mult exon skipping, intron retention
        return ((positions[1], positions[4]),)
    elif event_id.startswith('alt_5prime'):
        if strand is '+':
            return ((positions[3], positions[0]),
                    (positions[5], positions[0]))
        elif strand is '-':
            return ((positions[1], positions[2]),
                    (positions[1], positions[4]))
        else:
            return ((0, 0), (0, 0))
    elif event_id.startswith('alt_3prime'):
        if strand is '+':
            return ((positions[1], positions[2]),
                    (positions[1], positions[4]))
        elif strand is '-':
            return ((positions[3], positions[0]),
                    (positions[5], positions[0]))
        else:
            return ((0, 0), (0, 0))
    else:
        return (0, 0)

def add_alternative_donor_acceptor(gene_events, ev_name, ev_num, positions):
    if positions[0] < positions[1] < positions[2]:
        gene_events[(ev_name, ev_num, 'i')] = (positions[0], positions[2])
        gene_events[(ev_name, ev_num, 'e')] = (positions[1], positions[2])
    else:
        gene_events[(ev_name, ev_num, 'e')] = (positions[0], positions[2])
        gene_events[(ev_name, ev_num, 'i')] = (positions[1], positions[2])
    return gene_events

def parse_gene_events(args):
    gene_events = {}
    with open(os.path.join(ANNOTATION_DIR,
                           args.chromosome,
                           args.gene + '.events')) as gene_events_file:
        for line in gene_events_file:
            _, _, *positions, ev_id = line.strip('\n').split(' ')
            positions = tuple(int(p) for p in positions)
            ev_type, ev_num = ev_id[:2], ev_id[2:]
            if not ev_type.startswith('A'):
                gene_events[(ev_type, ev_num)] = (positions[0], positions[1])
            else:
                gene_events = add_alternative_donor_acceptor(gene_events, ev_type,
                                                             ev_num, positions)
    return gene_events

def parse_exps_events(args, gene_events):
    exps_events = {}
    with open(os.path.join(ANNOTATION_DIR,
                           args.chromosome,
                           args.gene, 'equals.info')) as exps_mapping:
        for line in exps_mapping:
            exp_name, *exp_events_list = line.strip('\n').split(' ')
            exps_events[exp_name] = ()
            for alt_event_id in exp_events_list:
                ev_name, ev_num = alt_event_id[:2], alt_event_id[2:]
                if ev_name == "A5" or ev_name == "A3":
                    alt_ss_type = ev_num[-1:] # [i]nner or [e]xtern
                    ev_num = ev_num[:-1]
                    exps_events[exp_name] += (gene_events[(ev_name, ev_num, alt_ss_type)],)
                else:
                    exps_events[exp_name] += (gene_events[(ev_name, ev_num)],)
    return exps_events    

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-c', '--chr',         dest='chromosome', help='Chromosome name')
    parser.add_argument('-g', '--gene',        dest='gene',       help='Gene name')
    parser.add_argument('-e', '--exp',         dest='exp',        help='Exp name')
    parser.add_argument('-s', '--sample-size', dest='sample',     help='Sample Size')

    args = parser.parse_args()

    gene_events = parse_gene_events(args)
    exps_events = parse_exps_events(args, gene_events)

    print("TOOL,TP,FP,FN,PREC,REC")
    tp, fp, fn = evaluate_asgal(   args, gene_events, exps_events[args.exp])
    print("asgal,{},{},{},{:>.5},{:>.5}".format(tp, fp, fn, tp / max((tp + fn), 1),
                                                tp / max((tp + fp), 1)))
    tp, fp, fn = evaluate_spladder(args, gene_events, exps_events[args.exp])
    print("spladder,{},{},{},{:>.5},{:>.5}".format(tp, fp, fn, tp / max((tp + fn), 1),
                                                   tp / max((tp + fp), 1)))
    tp, fp, fn = evaluate_rmats(   args, gene_events, exps_events[args.exp])
    if tp == fp == 0:
        fp = 1
    print("rmats,{},{},{},{:>.5},{:>.5}".format(tp, fp, fn, tp / max((tp + fn), 1),
                                                   tp / max((tp + fp), 1)))
    evaluate_majiq(   args, gene_events, exps_events[args.exp])
    # evaluate_suppa()
    # evaluate_leafcutter()

if __name__ == "__main__":
    main()
