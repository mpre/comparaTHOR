#!/usr/bin/env python3

import argparse
import sys
import os
import csv
import glob

ANNOTATION_DIR='/home/prj_gralign/data/BMC-exp-sim/NewAnnos/'
RESULTS_DIR='/home/prj_gralign/data/BMC-exp-sim/Results/'

EVENT_TYPES = ('A5', 'A3', 'IR', 'ES')

spladder_ev_map = {'alt_5prime'       : 'A5',
                   'alt_3prime'       : 'A3',
                   'exon_skip'        : 'ES',
                   'intron_retention' : 'IR'}

rmats_ev_map = { 'A5SS' : 'A5',
                 'A3SS' : 'A3',
                 'SE'   : 'ES',
                 'RI'   : 'IR'}

def evaluate_asgal(args, gene_events, exp_events):
    asgal_introns = {ev_type : () for ev_type in EVENT_TYPES}
    with open(os.path.join(RESULTS_DIR, args.sample, args.chromosome, 'asgal',
                           args.gene,   args.exp,    args.exp + '.events.csv'),
              newline='') as asgal_csv:
        results = csv.reader(asgal_csv)
        next(results, None) # Drop header
        for row in results:
            asgal_introns[row[0]] += ((int(row[1]) - 1, int(row[2]) + 1),)

    results = {ev_type : {} for ev_type in EVENT_TYPES}
    for ev_type in EVENT_TYPES:
        nelems = len(set(exp_events[ev_type]))
        tp=len(set(asgal_introns[ev_type]) & set(exp_events[ev_type]))
        fp=len(set(asgal_introns[ev_type]) - set(exp_events[ev_type]) - set(gene_events[ev_type].values()))
        fn=len(set(exp_events[ev_type]) - set(asgal_introns[ev_type]))
        results[ev_type] = {'nelems' : nelems, 'tp' : tp, 'fp' : fp, 'fn' : fn}
    return results

def evaluate_spladder(args, gene_events, exp_events):
    spladder_introns = {ev_type : () for ev_type in EVENT_TYPES}
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
                spladder_introns = spladder_parse_line(spladder_introns,
                                                       event_id,
                                                       strand,
                                                       positions)
    results = {ev_type : {} for ev_type in EVENT_TYPES}
    for ev_type in EVENT_TYPES:
        nelems = len(set(exp_events[ev_type]))
        tp=len(set(spladder_introns[ev_type]) & set(exp_events[ev_type]))
        fp=len(set(spladder_introns[ev_type]) - set(exp_events[ev_type]) - set(gene_events[ev_type].values()))
        fn=len(set(exp_events[ev_type]) - set(spladder_introns[ev_type]))
        results[ev_type] = {'nelems' : nelems, 'tp' : tp, 'fp' : fp, 'fn' : fn}
    return results

def evaluate_rmats(args, gene_events, exp_events):
    # We do not consider MXE
    rmats_introns = {ev_type : () for ev_type in EVENT_TYPES}
    rmats_files_list = glob.glob(os.path.join(RESULTS_DIR, args.sample, args.chromosome,
                                              'rMATS',     args.gene,   args.exp,
                                              'fromGTF.[!(n,M)]*.txt'))
    for rmats_file in rmats_files_list:
        with open(rmats_file) as rmats_current_event_file:
            next(rmats_current_event_file) # Drop header
            rmats_ev_type = rmats_file[rmats_file.rfind('.', 0, rmats_file.rfind('.')) + 1:
                                       rmats_file.rfind('.')]
            ev_type = rmats_ev_map[rmats_ev_type]
            for row in rmats_current_event_file:
                _, _, _, _, strand, *positions = row.strip('\n').split('\t')
                positions = tuple(int(p) for p in positions)
                rmats_introns = rmats_parse_line(rmats_introns,
                                                 ev_type,
                                                 strand,
                                                 positions)
    results = {ev_type : {} for ev_type in EVENT_TYPES}
    for ev_type in EVENT_TYPES:
        nelems = len(set(exp_events[ev_type]))
        tp=len(set(rmats_introns[ev_type]) & set(exp_events[ev_type]))
        fp=len(set(rmats_introns[ev_type]) - set(exp_events[ev_type]) - set(gene_events[ev_type].values()))
        fn=len(set(exp_events[ev_type]) - set(rmats_introns[ev_type]))
        results[ev_type] = {'nelems' : nelems, 'tp' : tp, 'fp' : fp, 'fn' : fn}
    return results

def evaluate_majiq(args, gene_events, exp_events):
    return

def spladder_parse_line(spladder_introns, event_id, strand, positions):
    ev_type = spladder_ev_map[event_id[0:event_id.rfind('_')]]
    if ev_type in ['ES', 'IR']:
        # Exon skipping, mult exon skipping, intron retention
        spladder_introns[ev_type] += ((positions[1], positions[4]),)
    elif ev_type is 'A5':
        if strand is '+':
            spladder_introns[ev_type] += ((positions[3], positions[0]),
                                          (positions[5], positions[0]))
        elif strand is '-':
            spladder_introns[ev_type] += ((positions[1], positions[2]),
                                          (positions[1], positions[4]))
    elif ev_type is 'A3':
        if strand is '+':
            spladder_introns[ev_type] += ((positions[1], positions[2]),
                                          (positions[1], positions[4]))
        elif strand is '-':
            spladder_introns[ev_type] += ((positions[3], positions[0]),
                                          (positions[5], positions[0]))
    return spladder_introns

def rmats_parse_line(rmats_introns, ev_type, strand, positions):
    if ev_type in ['ES', 'IR']:
        # Exon skipping, mult exon skipping, intron retention
        rmats_introns[ev_type] += ((positions[3], positions[4] + 1),)
    elif ev_type is 'A5':
        if strand is '+':
            rmats_introns[ev_type] += ((positions[1], positions[4] + 1),
                                       (positions[3], positions[4] + 1))
        elif strand is '-':
            rmats_introns[ev_type] += ((positions[5], positions[0] + 1),
                                       (positions[5], positions[2] + 1))
    elif ev_type is 'A3':
        if strand is '+':
            rmats_introns[ev_type] += ((positions[5], positions[0] + 1),
                                       (positions[5], positions[2] + 1))
        elif strand is '-':
            rmats_introns[ev_type] += ((positions[1], positions[4] + 1),
                                       (positions[3], positions[4] + 1))
    return rmats_introns

def add_alternative_donor_acceptor(gene_events, ev_type, ev_num, positions):
    if positions[0] < positions[1] < positions[2]:
        gene_events[ev_type][(ev_type, ev_num, 'i')] = (positions[0], positions[2])
        gene_events[ev_type][(ev_type, ev_num, 'e')] = (positions[1], positions[2])
    else:
        gene_events[ev_type][(ev_type, ev_num, 'e')] = (positions[0], positions[2])
        gene_events[ev_type][(ev_type, ev_num, 'i')] = (positions[1], positions[2])
    return gene_events

def parse_gene_events(args):
    gene_events = {ev_type : {} for ev_type in EVENT_TYPES}
    with open(os.path.join(ANNOTATION_DIR,
                           args.chromosome,
                           args.gene + '.events')) as gene_events_file:
        for line in gene_events_file:
            _, ev_type, *positions, ev_id = line.strip('\n').split(' ')
            positions = tuple(int(p) for p in positions)
            _, ev_num = ev_id[:2], ev_id[2:]
            if not ev_type.startswith('A'):
                gene_events[ev_type][(ev_type, ev_num)] = (positions[0], positions[1])
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
            if exp_name not in exps_events:
                exps_events[exp_name] = {ev_type : () for ev_type in EVENT_TYPES}
            for alt_event_id in exp_events_list:
                ev_type, ev_num = alt_event_id[:2], alt_event_id[2:]
                if ev_type == "A5" or ev_type == "A3":
                    alt_ss_type = ev_num[-1:] # [i]nner or [e]xtern
                    ev_num = ev_num[:-1]
                    exps_events[exp_name][ev_type] += (gene_events[ev_type][(ev_type, ev_num, alt_ss_type)],)
                else:
                    exps_events[exp_name][ev_type] += (gene_events[ev_type][(ev_type, ev_num)],)
    return exps_events    

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-c', '--chr',         dest='chromosome', help='Chromosome name', required=True)
    parser.add_argument('-g', '--gene',        dest='gene',       help='Gene name',       required=True)
    parser.add_argument('-e', '--exp',         dest='exp',        help='Exp name',        required=True)
    parser.add_argument('-s', '--sample-size', dest='sample',     help='Sample Size',     required=True)

    args = parser.parse_args()

    gene_events = parse_gene_events(args)
    exps_events = parse_exps_events(args, gene_events)

    print("TOOL,EV_TYPE,NELEMS,TP,FP,FN")

    results = evaluate_asgal(args, gene_events, exps_events[args.exp])
    for ev_type in EVENT_TYPES:
        print("asgal,{},{},{},{},{}".format(ev_type, results[ev_type]['nelems'], results[ev_type]['tp'],
                                            results[ev_type]['fp'], results[ev_type]['fn']))
    
    results = evaluate_spladder(args, gene_events, exps_events[args.exp])
    for ev_type in EVENT_TYPES:
        print("spladder,{},{},{},{},{}".format(ev_type, results[ev_type]['nelems'], results[ev_type]['tp'],
                                               results[ev_type]['fp'], results[ev_type]['fn']))

    results =  evaluate_rmats(args, gene_events, exps_events[args.exp])
    for ev_type in EVENT_TYPES:
        print("spladder,{},{},{},{},{}".format(ev_type, results[ev_type]['nelems'], results[ev_type]['tp'],
                                               results[ev_type]['fp'], results[ev_type]['fn']))

#     # evaluate_majiq(   args, gene_events, exps_events[args.exp])
#     evaluate_suppa()
#     # evaluate_leafcutter()

if __name__ == "__main__":
    main()

