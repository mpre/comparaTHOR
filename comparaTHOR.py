#!/usr/bin/env python3

import argparse
import sys
import os
import csv

ANNOTATION_DIR='/home/prj_gralign/data/BMC-exp-sim/NewAnnos/'
RESULTS_DIR='/home/prj_gralign/data/BMC-exp-sim/old-Results/'

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
    print("TP,FP,FN,PREC,REC")
    print("{:},{},{},{:>.5},{:>.5}".format(tp, fp, fn, tp / (tp + fn), tp / (tp + fp)))

    return tp, fp, fn

def evaluate_spladder():
    return

def evaluate_rmats():
    return

def evaluate_majiq():
    return

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

    evaluate_asgal(args, gene_events, exps_events[args.exp])
    evaluate_spladder()
    evaluate_rmats()
    evaluate_majiq()
    # evaluate_suppa()
    # evaluate_leafcutter()

if __name__ == "__main__":
    main()
