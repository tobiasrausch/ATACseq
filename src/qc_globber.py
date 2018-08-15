#! /usr/bin/env python

from __future__ import print_function
import argparse
import sys
import gzip
import os
import csv

def linecount(fname):
    i = 0
    if os.path.isfile(fname):
        f = gzip.open(fname) if fname.endswith('.gz') else open(fname)
        for i, l in enumerate(f):
            pass
    return i + 1

# Parse command line
parser = argparse.ArgumentParser(description='Key Quality Control Metrics')
parser.add_argument('-p', '--prefix', metavar='prefix', required=True, dest='prefix', help='output prefix used in ATAC-Seq pipeline (required)')
args = parser.parse_args()

# Collect QC statistics
qc = dict()
if args.prefix:
    # Unfiltered BAM stats
    filename = args.prefix + ".bamStats.unfiltered.tsv.gz"
    if os.path.isfile(filename):
        with gzip.open(filename, mode="rt") as f:
            header = dict()
            headerCM = dict()
            f_reader = csv.reader(f, delimiter="\t")
            for row in f_reader:
                if row[0].startswith('ME'):
                    if not len(header.keys()):
                        for idx, key in enumerate(row):
                            header[key] = idx
                    else:
                        qc['UnmappedFraction'] = row[header['UnmappedFraction']]
                        qc['DuplicateFraction'] = row[header['DuplicateFraction']]
                        qc['MappedSameChrFraction'] = row[header['MappedSameChrFraction']]
                if row[0].startswith('CM'):
                    if not len(headerCM.keys()):
                        for idx, key in enumerate(row):
                            headerCM[key] = idx
                    else:
                        if (row[headerCM['Chrom']] == "M") or (row[headerCM['Chrom']] == "chrM") or (row[headerCM['Chrom']] == "MT") or (row[headerCM['Chrom']] == "chrMT"):
                            qc['FractionChrM'] = row[headerCM['MappedFraction']]
                        
    # Filtered BAM stats
    filename = args.prefix + ".bamStats.promoters.tsv.gz"
    if os.path.isfile(filename):
        with gzip.open(filename, mode="rt") as f:
            header = dict()
            f_reader = csv.reader(f, delimiter="\t")
            for row in f_reader:
                if row[0].startswith('ME'):
                    if not len(header.keys()):
                        for idx, key in enumerate(row):
                            header[key] = idx
                    else:
                        qc['Sample'] = row[header['Sample']].replace(".final","")
                        qc['MappedReads'] = row[header['#Mapped']]
                        qc['ErrorRate'] = row[header['ErrorRate']]
                        qc['SDCoverage'] = row[header['SDCoverage']]
                        qc['BpCov1ToCovNRatio'] = row[header['BpCov1ToCovNRatio']]
                        qc['BpCov1ToCov2Ratio'] = row[header['BpCov1ToCov2Ratio']]
                        qc['TssEnrichment'] = row[header['EnrichmentOverBed']]

    # Peak calling statistics
    qc['UnfilteredPeaks'] = str(linecount(args.prefix + ".unfiltered.peaks.gz"))
    qc['FilteredPeaks'] = str(linecount(args.prefix + ".peaks"))
    qc['FractionPeaksRetained'] = str(float(qc['FilteredPeaks'])/float(qc['UnfilteredPeaks']))
    filename = args.prefix + ".peaks.log"
    if os.path.isfile(filename):
        with open(filename, 'r') as f:
            f_reader = csv.DictReader(f, delimiter="\t")
            for row in f_reader:
                qc['FRiP'] = row['frip']
                qc['PeakSaturation'] = str(min(float(row['recallRep1']), float(row['recallRep2'])))

# Output summary QC information
cols = sorted(qc.keys())
print('\t'.join(cols))
first = True
for c in cols:
    if first:
        first = False
    else:
        print("\t", end="")
    print(qc[c], end="")
print()
