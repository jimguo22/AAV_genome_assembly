# This script assembles plus strand of AAV-MTM1 recombinant genome. The first searching tag is D sequence plus the first base that is different between begining and end(reverse-complmentary) of the recombinant genome. The reference sequence is built-into the script to reduce the complexity of operation. 
# python assmebler.py read1_fastq_file report.txt

from sys import argv
import time
import re
import itertools
#package for pyton2 and 3 compatibility
from __future__ import print_function
import future 
import builtins
import past 
import six  
## Define functions

### This function convert the ascII code to numerical quality score 
def phred33ToQ(qual):
    return ord(qual)-33

def reverse_complement(seq_tag):
    new_seq = ''
    reverse_compl = {"A":"T", "C":"G", "G":"C", "T":"A", "N":"N"}
    for i in range(len(seq_tag)):
        new_seq =reverse_compl[seq_tag[i]]+new_seq 
    return new_seq

def find_next_nt(seq_tag, read1_seq, read2_seq, read1_phred, read2_phred, read1_name, read2_name):
    next_nt,next_nt_final, count_sorted, read_names = ([] for i in range(4))
    next_nt_1,next_nt_2,len_next_nt = ('' for i in range(3))
    count = {'A':0, 'C':0, 'G':0, 'T':0, 'NA':0}
    ratio_1,ratio_2,ratio_3, n, i = (0 for i in range(5))
    seq_tag_rc = reverse_complement(seq_tag)
    for i, seq in enumerate(read1_seq):   
        for m in re.finditer(seq_tag, seq):
            if m.end()<len(seq):   
                next_nt.append([seq[m.end()],read1_phred[i][m.end()]])
                read_names.append(read1_name[i])
    for i, seq in enumerate(read2_seq):
        for n in re.finditer(seq_tag_rc, seq):
            if n.start()>0:
                next_nt.append([reverse_complement(seq[n.start()-1]),read2_phred[i][n.start()+1]])
                read_names.append(read2_name[i])
                
                # print next_nt
    #print "total # of next nt is:"
    #print len(next_nt)
    len_next_nt = len(next_nt)
    # convert phred score to numeric quality score
    for i, line in enumerate(next_nt):
        line[1]=phred33ToQ(line[1])
        # set quality score threshold
        if int(line[1]) > 30:
            next_nt_final.append(line)
        #print next_nt_final
    # generate key:value pairs equal to nt:frequency, case-insensitive
    for line in next_nt_final:
        if line[0] == 'A' or line[0] == 'a':
            count['A'] +=1
        elif line[0] == 'C' or line[0] == 'c':
            count['C'] +=1
        elif line[0] == 'G' or line[0] == 'g':
            count['G'] +=1
        elif line[0] == 'T' or line[0] == 't':
            count['T'] +=1
        else:
            count['NA'] +=1
    count_sorted = sorted(count, key=count.get, reverse=True)
    next_nt_1 = count_sorted[0]
 
    if (count['A']+count['C']+count['G']+count['T']+count['NA']) != 0:
        ratio_1 = float(count[count_sorted[0]])/(count['A']+count['C']+count['G']+count['T']+count['NA'])*100
        ratio_2 = float(count[count_sorted[1]]+count[count_sorted[2]]+count[count_sorted[3]]+count[count_sorted[4]])/(count['A']+count['C']+count['G']+count['T']+count['NA'])*100
        ratio_3 = float(count[count_sorted[1]])/(count['A']+count['C']+count['G']+count['T']+count['NA'])*100
        # print "next most abundant nt is %s and account for %.2f of all nts" %(next_nt_1, ratio_1)
        # print "2nd most abundant nt is %s and account for %.2f of all nts" %(count_sorted[1], ratio_3)
    return next_nt_1, ratio_1, ratio_2, read_names, len_next_nt            

def assembler(seq_tag, last_tag,r1_seq, r2_seq,r1_phred,r2_phred, r1_name, r2_name, reference_seq):
    m,per_reads_used = (0 for i in range(2))
    next_seq_tag = ''
    assembled_seq = seq_tag
    maximum_len = len(reference_seq)
    poly_T_tail = 'T'*21
    next_nt_1, ratio_1, ratio_2, read_names, coverage, read_names_list, coverage_list,error_rate,report = ([] for i in range(9))
    while 1:
        next_nt_1, ratio_1, ratio_2, read_names, coverage = find_next_nt(seq_tag,r1_seq, r2_seq,r1_phred,r2_phred, r1_name, r2_name)
        read_names_list = read_names+read_names_list
        # print "read_names_list is " 
        # print read_names_list
        # print len(read_names_list)
        coverage_list.append(coverage)
        if ratio_1 == 0:
            print "Error: no more nt could be added."
            report.append("Error: no more nt could be added.")
            break
        assembled_seq +=next_nt_1
        if len(assembled_seq) > len(reference_seq):
            print "Error: assembled sequence is longer than reference sequence."
            report.append("Error: assembled sequence is longer than reference sequence.")
            break
        next_seq_tag = seq_tag[1:]+next_nt_1        
        if next_seq_tag == last_tag:
            print "last tag reached! Sucessful assembly!"
            break
        elif next_seq_tag == poly_T_tail:
            print ("A 21nt polyT sequence was assembled. The assembly process was not successful, probably because the quality of sequencing reads is too low.")
            report.append("A 21nt polyT sequence was assembled. The assembly process was not successful, probably because the quality of sequencing reads is too low.")
        else:
            seq_tag = next_seq_tag
        # print "the assembled_seq is:"
        # print assembled_seq
        error_rate.append(ratio_2)
    if len(assembled_seq) != len(reference_seq):
        report.append("The length of the assembled sequence is %s, different from the length of reference sequence %s.\n" %(len(assembled_seq), len(reference_seq)))
        report.append("The final assembled sequence is %s.\n" %assembled_seq)
        m = -1
    for i in range(len(assembled_seq)):   
        if assembled_seq[i] != reference_seq[i]:
            m +=1
    if m == 0:
        report.append ("The assembled sequence matches the reference sequence 100%!\n")
    elif m>0:
        report.append("There are %d mismatches between assembled and reference sequence.\n" %m)
        report.append("The final assembled sequence is %s\n" %assembled_seq)
    read_names_list = list(set(read_names_list))
    # print len(read_names_list)
    # print len(r1_name)
    # calculate the percentage of reads used for mapping
    per_reads_used = float(len(read_names_list))/len(r1_name)*100
    report.append("%.2f percentage of total reads used for denovo assembly.\n" %per_reads_used) 
    return(error_rate, coverage_list, report)

## Main
t0 = time.time()

r1 = argv[1]
r2 = argv[2]
final_report = argv[3]
error_rate_fw, coverage_list_fw, report_fw, error_rate_rv, coverage_list_rv, report_rv, r1_name, r1_seq, r1_plus, r1_phred,r2_name, r2_seq, r2_plus, r2_phred, = ([] for i in range(14))

with open(r1, 'r') as f1:
    lines = f1.read().splitlines()
    # read the first 1 M reads for downstream analysis
    for i, line in enumerate(lines):
        if i%4 == 0:
            r1_name.append(line)
        elif i%4 == 1:
            r1_seq.append(line)
        elif i%4 == 2:
            r1_plus.append(line)
        elif i%4 ==3:
            r1_phred.append(line)

with open(r2, 'r') as f2:
    lines = f2.read().splitlines()
    # read the first 1 M reads for downstream analysis
    for i, line in enumerate(lines):
        if i%4 == 0:
            r2_name.append(line)
        elif i%4 == 1:
            r2_seq.append(line)
        elif i%4 == 2:
            r2_plus.append(line)
        elif i%4 ==3:
            r2_phred.append(line)

#seq_tag for actural run
seq_tag_fw ="TACCCCCTGCCCCCCACAGCT"
#seq_tag for testing
#seq_tag_fw = "ACTGTCCTGTGAGCCCTTCTT"
last_tag_fw = 'GCCTCCCCCACTCACAGTGAC'
reference_seq_fw = "TACCCCCTGCCCCCCACAGCTCCTCTCCTGTGCCTTGTTTCCCAGCCATGCGTTCTCCTCTATAAATACCCGCTCTGGTATTTGGGGTTGGCAGCTGTTGCTGCCAGGGAGATGGTTGGGTTGACATGCGGCTCCTGACAAAACACAAACCCCTGGTGTGTGTGGGCGTGGGTGGTGTGAGTAGGGGGATGAATCAGGGAGGGGGCGGGGGACCCAGGGGGCAGGAGCCACACAAAGTCTGTGCGGGGGTGGGAGCGCACATAGCAATTGGAAACTGAAAGCTTATCAGACCCTTTCTGGAAATCAGCCCACTGTTTATAAACTTGAGGCCCCACCCTCGACAGTACCGGGGAGGAAGAGGGCCTGCACTAGTCCAGAGGGAAACTGAGGCTCAGGGCCAGCTCGCCCATAGACATACATGGCAGGCAGGCTTTGGCCAGGATCCCTCCGCCTGCCAGGCGTCTCCCTGCCCTCCCTTCCTGCCTAGAGACCCCCACCCTCAAGCCTGGCTGGTCTTTGCCTGAGACCCAAACCTCTTCGACTTCAAGAGAATATTTAGGAACAAGGTGGTTTAGGGCCTTTCCTGGGAACAGGCCTTGACCCTTTAAGAAATGACCCAAAGTCTCTCCTTGACCAAAAAGGGGACCCTCAAACTAAAGGGAAGCCTCTCTTCTGCTGTCTCCCCTGACCCCACTCCCCCCCACCCCAGGACGAGGAGATAACCAGGGCTGAAAGAGGCCCGCCTGGGGGCTGCAGACATGCTTGCTGCCTGCCCTGGCGAAGGATTGGTAGGCTTGCCCGTCACAGGACCCCCGCTGGCTGACTCAGGGGCGCAGGCCTCTTGCGGGGGAGCTGGCCTCCCCGCCCCCACGGCCACGGGCCGCCCTTTCCTGGCAGGACAGCGGGATCTTGCAGCTGTCAGGGGAGGGGAGGCGGGGGCTGATGTCAGGAGGGATACAAATAGTGCCGACGGCTGGGGGCCCTGTCTCCCCTCGCCGCATCCACTCTCCGGCCGGCCGCCTGCCCGCCGCCTCCTCCGTGCGCCCGCCAGCCTCGCCCGGACTCTAGAGGATCCAGATCTAAGCTTCTCTGGTCACCGATCCTGAGAACTTCAGGGTGAGTCTATGGGACCCTTGATGTTTTCTTTCCCCTTCTTTTCTATGGTTAAGTTCATGTCATAGGAAGGGGAGAAGTAACAGGGTACACATATTGACCAAATCAGGGTAATTTTGCATTTGTAATTTTAAAAAATGCTTTCTTCTTTTAATATACTTTTTTGTTTATCTTATTTCTAATACTTTCCCTAATCTCTTTCTTTCAGGGCAATAATGATACAATGTATCATGCCTCTTTGCACCATTCTAAAGAATAACAGTGATAATTTCTGGGTTAAGGCAATAGCAATATTTCTGCATATAAATATTTCTGCATATAAATTGTAACTGATGTAAGAGGTTTCATATTGCTAATAGCAGCTACAATCCAGCTACCATTCTGCTTTTATTTTATGGTTGGGATAAGGCTGGATTATTCTGAGTCCAAGCTAGGCCCTTTTGCTAATCATGTTCATACCTCTTATCTTCCTCCCACAGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCCGCGGGCGGCCGCAAGTTTCCAGGATGGCTTCTGCATCAACTTCTAAATATAATTCACACTCCTTGGAGAATGAGTCTATTAAGAGGACGTCTCGAGATGGAGTCAATCGAGATCTCACTGAGGCTGTTCCTCGACTTCCAGGAGAAACACTAATCACTGACAAAGAAGTTATTTACATATGTCCTTTCAATGGCCCCATTAAGGGAAGAGTTTACATCACAAATTATCGTCTTTATTTAAGAAGTTTGGAAACGGATTCTTCTCTAATACTTGATGTTCCTCTGGGTGTGATCTCGAGAATTGAAAAAATGGGAGGCGCGACAAGTAGAGGAGAAAATTCCTATGGTCTAGATATTACTTGTAAAGACATGAGAAACCTGAGGTTCGCTTTGAAACAGGAAGGCCACAGCAGAAGAGATATGTTTGAGATCCTCACGAGATACGCGTTTCCCCTGGCTCACAGTCTGCCATTATTTGCATTTTTAAATGAAGAAAAGTTTAACGTGGATGGATGGACAGTTTACAATCCAGTGGAAGAATACAGGAGGCAGGGCTTGCCCAATCACCATTGGAGAATAACTTTTATTAATAAGTGCTATGAGCTCTGTGACACTTACCCTGCTCTTTTGGTGGTTCCGTATCGTGCCTCAGATGATGACCTCCGGAGAGTTGCAACTTTTAGGTCCCGAAATCGAATTCCAGTGCTGTCATGGATTCATCCAGAAAATAAGACGGTCATTGTGCGTTGCAGTCAGCCTCTTGTCGGTATGAGTGGGAAACGAAATAAAGATGATGAGAAATATCTCGATGTTATCAGGGAGACTAATAAACAAATTTCTAAACTCACCATTTATGATGCAAGACCCAGCGTAAATGCAGTGGCCAACAAGGCAACAGGAGGAGGATATGAAAGTGATGATGCATATCATAACGCCGAACTTTTCTTCTTAGACATTCATAATATTCATGTTATGCGGGAATCTTTAAAAAAAGTGAAGGACATTGTTTATCCTAATGTAGAAGAATCTCATTGGTTGTCCAGTTTGGAGTCTACTCATTGGTTAGAACATATCAAGCTCGTTTTGACAGGAGCCATTCAAGTAGCAGACAAAGTTTCTTCAGGGAAGAGTTCAGTGCTTGTGCATTGCAGTGACGGATGGGACAGGACTGCTCAGCTGACATCCTTGGCCATGCTGATGTTGGATAGCTTCTATAGGAGCATTGAAGGGTTCGAAATACTGGTACAAAAAGAATGGATAAGTTTTGGACATAAATTTGCATCTCGAATAGGTCATGGTGATAAAAACCACACCGATGCTGACCGTTCTCCTATTTTTCTCCAGTTTATTGATTGTGTGTGGCAAATGTCAAAACAGTTCCCTACAGCTTTTGAATTCAATGAACAATTTTTGATTATAATTTTGGATCATCTGTATAGTTGCCGATTTGGTACTTTCTTATTCAACTGTGAATCTGCTCGAGAAAGACAGAAGGTTACAGAAAGGACTGTTTCTTTATGGTCACTGATAAACAGTAATAAAGAAAAATTCAAAAACCCCTTCTATACTAAAGAAATCAATCGAGTTTTATATCCAGTTGCCAGTATGCGTCACTTGGAACTCTGGGTGAATTACTACATTAGATGGAACCCCAGGATCAAGCAACAACAGCCGAATCCAGTGGAGCAGCGTTACATGGAGCTCTTAGCCTTACGCGACGAATACATAAAGCGGCTTGAGGAACTGCAGCTCGCCAACTCTGCCAAGCTTTCTGATCCCCCAACTTCACCTTCCAGTCCTTCGCAAATGATGCCCCATGTGCAAACTCACTTCTGACCGGTCCGAGGGCCCAGATCTAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAAGCTCGCTTTCTTGCTGTCCAATTTCTATTAAAGGTTCCTTTGTTCCCTAAGTCCAACTACTAAACTGGGGGATATTATGAAGGGCCTTGAGCATCTGGATTCTGCCTAATAAAAAACATTTATTTTCATTGCAATGATGTATTTAAATTATTTCTGAATATTTTACTAAAAAGGGAATGTGGGAGGTCAGTGCATTTAAAACATAAAGAAATGAAGAGCTAGTTCAAACCTTGGGAAAATACACTATATCTTAAACTCCATGAAAGAAGGTGAGGCTGCAAACAGCTAATGCACATTGGCAACAGCCCCTGATGCCTATGCCTTATTCATCCCTCAGAAAAGGATTCAAGTAGAGGCTTGATTTGGAGGTTAAAGTTTTGCTATGCTGTATTTTACATTACTTATTGTTTTAGCTGTCCTCATGAATGTCTTTTCACTACCCATTTGCTTATCCTGCATCTCTCAGCCTTGACTCCACTCAGTTCTCTTGCTTAGAGATACCACCTTTCCCCTGAAGTGTTCCTTCCATGTTTTACGGCGAGATGGTTTCTCCTCGCCTGGCCACTCAGCCTTAGTTGTCTCTGTTGTCTTATAGAGGTCTACTTGAAGAAGGAAAAACAGGGGGCATGGTTTGACTGTCCTGTGAGCCCTTCTTCCCTGCCTCCCCCACTCACAGTGAC"

#seq_tag for actural run
seq_tag_rv ="GTCACTGTGAGTGGGGGAGGC"
#seq_tag for testing
#seq_tag_rv = "GAGAACGCATGGCTGGGAAAC"
last_tag_rv = 'AGCTGTGGGGGGCAGGGGGTA'
reference_seq_rv = "GTCACTGTGAGTGGGGGAGGCAGGGAAGAAGGGCTCACAGGACAGTCAAACCATGCCCCCTGTTTTTCCTTCTTCAAGTAGACCTCTATAAGACAACAGAGACAACTAAGGCTGAGTGGCCAGGCGAGGAGAAACCATCTCGCCGTAAAACATGGAAGGAACACTTCAGGGGAAAGGTGGTATCTCTAAGCAAGAGAACTGAGTGGAGTCAAGGCTGAGAGATGCAGGATAAGCAAATGGGTAGTGAAAAGACATTCATGAGGACAGCTAAAACAATAAGTAATGTAAAATACAGCATAGCAAAACTTTAACCTCCAAATCAAGCCTCTACTTGAATCCTTTTCTGAGGGATGAATAAGGCATAGGCATCAGGGGCTGTTGCCAATGTGCATTAGCTGTTTGCAGCCTCACCTTCTTTCATGGAGTTTAAGATATAGTGTATTTTCCCAAGGTTTGAACTAGCTCTTCATTTCTTTATGTTTTAAATGCACTGACCTCCCACATTCCCTTTTTAGTAAAATATTCAGAAATAATTTAAATACATCATTGCAATGAAAATAAATGTTTTTTATTAGGCAGAATCCAGATGCTCAAGGCCCTTCATAATATCCCCCAGTTTAGTAGTTGGACTTAGGGAACAAAGGAACCTTTAATAGAAATTGGACAGCAAGAAAGCGAGCTTAGTGATACTTGTGGGCCAGGGCATTAGCCACACCAGCCACCACTTTCTGATAGGCAGCCTGCACTGGTGGGGTGAATTAGATCTGGGCCCTCGGACCGGTCAGAAGTGAGTTTGCACATGGGGCATCATTTGCGAAGGACTGGAAGGTGAAGTTGGGGGATCAGAAAGCTTGGCAGAGTTGGCGAGCTGCAGTTCCTCAAGCCGCTTTATGTATTCGTCGCGTAAGGCTAAGAGCTCCATGTAACGCTGCTCCACTGGATTCGGCTGTTGTTGCTTGATCCTGGGGTTCCATCTAATGTAGTAATTCACCCAGAGTTCCAAGTGACGCATACTGGCAACTGGATATAAAACTCGATTGATTTCTTTAGTATAGAAGGGGTTTTTGAATTTTTCTTTATTACTGTTTATCAGTGACCATAAAGAAACAGTCCTTTCTGTAACCTTCTGTCTTTCTCGAGCAGATTCACAGTTGAATAAGAAAGTACCAAATCGGCAACTATACAGATGATCCAAAATTATAATCAAAAATTGTTCATTGAATTCAAAAGCTGTAGGGAACTGTTTTGACATTTGCCACACACAATCAATAAACTGGAGAAAAATAGGAGAACGGTCAGCATCGGTGTGGTTTTTATCACCATGACCTATTCGAGATGCAAATTTATGTCCAAAACTTATCCATTCTTTTTGTACCAGTATTTCGAACCCTTCAATGCTCCTATAGAAGCTATCCAACATCAGCATGGCCAAGGATGTCAGCTGAGCAGTCCTGTCCCATCCGTCACTGCAATGCACAAGCACTGAACTCTTCCCTGAAGAAACTTTGTCTGCTACTTGAATGGCTCCTGTCAAAACGAGCTTGATATGTTCTAACCAATGAGTAGACTCCAAACTGGACAACCAATGAGATTCTTCTACATTAGGATAAACAATGTCCTTCACTTTTTTTAAAGATTCCCGCATAACATGAATATTATGAATGTCTAAGAAGAAAAGTTCGGCGTTATGATATGCATCATCACTTTCATATCCTCCTCCTGTTGCCTTGTTGGCCACTGCATTTACGCTGGGTCTTGCATCATAAATGGTGAGTTTAGAAATTTGTTTATTAGTCTCCCTGATAACATCGAGATATTTCTCATCATCTTTATTTCGTTTCCCACTCATACCGACAAGAGGCTGACTGCAACGCACAATGACCGTCTTATTTTCTGGATGAATCCATGACAGCACTGGAATTCGATTTCGGGACCTAAAAGTTGCAACTCTCCGGAGGTCATCATCTGAGGCACGATACGGAACCACCAAAAGAGCAGGGTAAGTGTCACAGAGCTCATAGCACTTATTAATAAAAGTTATTCTCCAATGGTGATTGGGCAAGCCCTGCCTCCTGTATTCTTCCACTGGATTGTAAACTGTCCATCCATCCACGTTAAACTTTTCTTCATTTAAAAATGCAAATAATGGCAGACTGTGAGCCAGGGGAAACGCGTATCTCGTGAGGATCTCAAACATATCTCTTCTGCTGTGGCCTTCCTGTTTCAAAGCGAACCTCAGGTTTCTCATGTCTTTACAAGTAATATCTAGACCATAGGAATTTTCTCCTCTACTTGTCGCGCCTCCCATTTTTTCAATTCTCGAGATCACACCCAGAGGAACATCAAGTATTAGAGAAGAATCCGTTTCCAAACTTCTTAAATAAAGACGATAATTTGTGATGTAAACTCTTCCCTTAATGGGGCCATTGAAAGGACATATGTAAATAACTTCTTTGTCAGTGATTAGTGTTTCTCCTGGAAGTCGAGGAACAGCCTCAGTGAGATCTCGATTGACTCCATCTCGAGACGTCCTCTTAATAGACTCATTCTCCAAGGAGTGTGAATTATATTTAGAAGTTGATGCAGAAGCCATCCTGGAAACTTGCGGCCGCCCGCGGAATTCTTTGCCAAAGTGATGGGCCAGCACACAGACCAGCACGTTGCCCAGGAGCTGTGGGAGGAAGATAAGAGGTATGAACATGATTAGCAAAAGGGCCTAGCTTGGACTCAGAATAATCCAGCCTTATCCCAACCATAAAATAAAAGCAGAATGGTAGCTGGATTGTAGCTGCTATTAGCAATATGAAACCTCTTACATCAGTTACAATTTATATGCAGAAATATTTATATGCAGAAATATTGCTATTGCCTTAACCCAGAAATTATCACTGTTATTCTTTAGAATGGTGCAAAGAGGCATGATACATTGTATCATTATTGCCCTGAAAGAAAGAGATTAGGGAAAGTATTAGAAATAAGATAAACAAAAAAGTATATTAAAAGAAGAAAGCATTTTTTAAAATTACAAATGCAAAATTACCCTGATTTGGTCAATATGTGTACCCTGTTACTTCTCCCCTTCCTATGACATGAACTTAACCATAGAAAAGAAGGGGAAAGAAAACATCAAGGGTCCCATAGACTCACCCTGAAGTTCTCAGGATCGGTGACCAGAGAAGCTTAGATCTGGATCCTCTAGAGTCCGGGCGAGGCTGGCGGGCGCACGGAGGAGGCGGCGGGCAGGCGGCCGGCCGGAGAGTGGATGCGGCGAGGGGAGACAGGGCCCCCAGCCGTCGGCACTATTTGTATCCCTCCTGACATCAGCCCCCGCCTCCCCTCCCCTGACAGCTGCAAGATCCCGCTGTCCTGCCAGGAAAGGGCGGCCCGTGGCCGTGGGGGCGGGGAGGCCAGCTCCCCCGCAAGAGGCCTGCGCCCCTGAGTCAGCCAGCGGGGGTCCTGTGACGGGCAAGCCTACCAATCCTTCGCCAGGGCAGGCAGCAAGCATGTCTGCAGCCCCCAGGCGGGCCTCTTTCAGCCCTGGTTATCTCCTCGTCCTGGGGTGGGGGGGAGTGGGGTCAGGGGAGACAGCAGAAGAGAGGCTTCCCTTTAGTTTGAGGGTCCCCTTTTTGGTCAAGGAGAGACTTTGGGTCATTTCTTAAAGGGTCAAGGCCTGTTCCCAGGAAAGGCCCTAAACCACCTTGTTCCTAAATATTCTCTTGAAGTCGAAGAGGTTTGGGTCTCAGGCAAAGACCAGCCAGGCTTGAGGGTGGGGGTCTCTAGGCAGGAAGGGAGGGCAGGGAGACGCCTGGCAGGCGGAGGGATCCTGGCCAAAGCCTGCCTGCCATGTATGTCTATGGGCGAGCTGGCCCTGAGCCTCAGTTTCCCTCTGGACTAGTGCAGGCCCTCTTCCTCCCCGGTACTGTCGAGGGTGGGGCCTCAAGTTTATAAACAGTGGGCTGATTTCCAGAAAGGGTCTGATAAGCTTTCAGTTTCCAATTGCTATGTGCGCTCCCACCCCCGCACAGACTTTGTGTGGCTCCTGCCCCCTGGGTCCCCCGCCCCCTCCCTGATTCATCCCCCTACTCACACCACCCACGCCCACACACACCAGGGGTTTGTGTTTTGTCAGGAGCCGCATGTCAACCCAACCATCTCCCTGGCAGCAACAGCTGCCAACCCCAAATACCAGAGCGGGTATTTATAGAGGAGAACGCATGGCTGGGAAACAAGGCACAGGAGAGGAGCTGTGGGGGGCAGGGGGTA"

print "start assembling the plus strand"
error_rate_fw, coverage_list_fw, report_fw=assembler(seq_tag_fw, last_tag_fw,r1_seq, r2_seq,r1_phred,r2_phred, r1_name, r2_name, reference_seq_fw)
print "start assembling the minus strand"
error_rate_rv, coverage_list_rv, report_rv=assembler(seq_tag_rv, last_tag_rv,r1_seq, r2_seq,r1_phred,r2_phred, r1_name, r2_name, reference_seq_rv)

# reverse the order or error rate and coverage for minus strand so the same coordinate could be used in both plus and minus strands in plot
error_rate_rv.reverse()
coverage_list_rv.reverse()


import matplotlib
# set default backend
matplotlib.use('Agg')

import matplotlib.pyplot as plt
f=plt.figure(1)
plt.subplot(211)
plt.plot(error_rate_fw)
plt.show()
plt.title("Base error rate (plus strand)")
plt.ylabel("% of non-reference sequence")
plt.subplot(212)
plt.plot(error_rate_rv)
plt.show()
plt.title("Base error rate (minus strand)")
plt.ylabel("% of non-reference sequence")
plt.xlabel("recombinant AAV genome")
f.savefig("error_rate.pdf")

g=plt.figure(2)
plt.subplot(211)
plt.plot(coverage_list_fw)
plt.show()
plt.title("Base coverage (plus strand)")
plt.ylabel("coverage")
plt.subplot(212)
plt.plot(coverage_list_rv)
plt.show()
plt.title("Base coverage (minus strand)")
plt.xlabel("recombinant AAV genome")
plt.ylabel("coverage")
g.savefig("coverage.pdf")

t1 = time.time()
total = t1-t0
report_rv.append("\n")
report_rv.append("\n")
report_rv.append("total running time is %.2f" % total)
with open (final_report, 'w') as f:
    f.write("\tSummary of MTM1 recombinant genome assembly (plus strand)\n")
    f.write("\n")
    f.write('\n')
    for line in report_fw:
        f.write(line)
    f.write("\n")
    f.write('\n')
    f.write("\tSummary of MTM1 recombinant genome assembly (minus strand)\n")
    f.write("\n")
    f.write('\n')
    for line in report_rv:
        f.write(line)
        
        