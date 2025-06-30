#!/usr/bin/env python
# -*- coding: utf-8 -*-

import edlib
import os
import sys, getopt

def add(a,b):
    both = (a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3], a[4] + b[4])
    return both

def mask(s, positions):  ##masks positions where similarity detected
    for pos in positions:
        length = pos[1]-pos[0]
        s = s[:pos[0]] + length*"N" + s[pos[1]:]
    return s

def resultc(result): ##counts non-overlapping loci detected by edlib with this editing distance
    prev=0
    count=0
    for a in result["locations"]:
        if a[0] > prev:
           count+=1
        prev=a[1]
    return count

def find_seq(seq, refone, dictofsteps, bed, id): #finds occurances of sequences in a chromosome that are closest to query_seq
    result = edlib.align(seq, refone.seq, mode = "HW", task = "locations")
    distance = result["editDistance"]
    if distance in dictofsteps:
       dictofsteps[distance].append({refone.id:resultc(result)})
    else:
       dictofsteps[distance] = [{refone.id:resultc(result)}]
    refone.seq = mask(refone.seq, result["locations"])
    for pos in result["locations"]:
        bed.write(refone.id + "\t" + str(pos[0]) + "\t" + str(pos[1]) + "\t" + id + "\t" + str(distance) + "\n")

def maxeditdistance(seq, ref, posbed, id): #returns the counts for edit distances 0-4
    refmasked = ref
    dictofsteps = {}
    for refone in refmasked:
        find_seq(seq, refone, dictofsteps, posbed, id)
    count=0
    stepsandcounts = {}
    for i in range(0,5):
        if i in dictofsteps:
            for a in dictofsteps[i]:
                for key in a:
                    count+=a[key]
            if i < 4:
                for b in dictofsteps[i]:
                    for k in b:
                        for record in refmasked:
                            if record.id == k:
                                find_seq(seq, record, dictofsteps, posbed, id)
        stepsandcounts[i] = count
    return stepsandcounts

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "hr:i:o:b:")
    except getopt.GetoptError:
        print("test.py [-S] -c <cram> -i <inputlist>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print("test.py -c <cram> -i <inputputlist>")
            sys.exit()
        elif opt in ("-i"):
            inputfile = arg
        elif opt in ("-o"):
            outputfile = arg
        elif opt in ("-r"):
            reference = arg
        elif opt in ("-b"):
            bedfile = arg

    records = list(SeqIO.parse(inputfile, "fasta"))
    ref=list(SeqIO.parse(reference, "fasta"))
    output = open(outputfile, 'w+')
    posbed = open(bedfile, 'w+')
    for record in records:
       print(record.id)
       dist = maxeditdistance(record.seq, ref, posbed, record.id)  #counts occurences for sequence with varying edit distances
       rev_dist = maxeditdistance(record.seq.reverse_complement(), ref, posbed, record.id) #counts occurences for reverse of the sequence with varying edit distances
       both = add(dist, rev_dist)
       output.write(str(record.id) + "\t" + str(record.seq) + "\t" + str(both[0]) + "\t" + str(both[1]) + "\t" + str(both[2]) + "\t" + str(both[3]) + "\t" + str(both[4]))

if __name__ == "__main__":
    main(sys.argv[1:])

