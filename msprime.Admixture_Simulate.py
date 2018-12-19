from __future__ import print_function
import sys
import msprime
import math
import numpy as np
import pybedtools
import random
import scipy
import subprocess
import glob
import collections
import gzip
import os.path
from scipy import spatial
from optparse import OptionParser
import msprime_demo_models

usage = str("usage: python <script> -p nonAfr -o Tenn -s 1 -i 2 -n 0.02 -d 0.0 -t 350 -c F4Dstat -l 1e6 -e 1006 -a 1008\n"+
            "WARNING: options 'haplo' and 'vcf' now write to stdout\n"+
            "WARNING: one-pulse option is deprecated; use -i 2 -n X -d 0.0 for comparable results")
parser = OptionParser(usage=usage)
parser.add_option("-p", "--population", action = "store", type = "string", dest = "pop", default="nonAfr", help="call introg_haplotypes in AFR, EUR, EAS, nonAfr, or all modHum; default=nonAfr")
parser.add_option("-o", "--outdir", action = "store", type = "string", dest = "outdir", default="Tenn", help="specify which demographic model to use by name; default=Tenn")
parser.add_option("-s", "--seed", action = "store", type = "int", dest = "seed", default=1, help="Set random seed for replicate chromosomes; default=1")
parser.add_option("-i", "--introgress_pulses", action = "store", type = "int", dest = "pulses", default=2, help="Set number of introgression pulses; default=2")
parser.add_option("-n", "--neand1_admixture_proportion", action = "store", type = "float", dest = "n1_admix_prop", default=0.02, help="Set N1 admixture proportion; default=0.02")
parser.add_option("-d", "--deni_admixture_propotion", action = "store", type = "float", dest = "n2_admix_prop", default=0.0, help="Set N1 2nd admixture proportion; default=0.0")
parser.add_option("--Neand1_sample_size", action = "store", type = "int", dest = "s_n1", default = 2, help="Set Neand1 sample size (haplotypes); default=2")
parser.add_option("--Neand2_sample_size", action = "store", type = "int", dest = "s_n2", default = 2, help="Set Neand2 sample size (haplotypes); default=2")
parser.add_option("-t", "--time_N1_N2_split", action = "store", type = "float", dest = "t_n1_n2", default=145, help="Set N1 N2 split time in kya; default=145")
parser.add_option("--time_N1_sample", action = "store", type = "float", dest = "t_n1_sample", default=55, help="Set N1 (Vindija) sample time in kya; default=55")
parser.add_option("--time_N2_sample", action = "store", type = "float", dest = "t_n2_sample", default=125, help="Set N2 (Altai) sample time in kya; default=125")
parser.add_option("-c", "--calls_or_stats", action="store", type = "string", dest = "haplo", default="F4Dstat", help="Pick haplotype calls 'haplo' or 'F4Dstat' or 'vcf' output, or 'debug'; default=F4Dstat")
parser.add_option("-l", "--length_chrom", action="store", type= "float", dest = "length", default=1e6, help="Define length of simulated chromosome; default=1e6")
parser.add_option("-e", "--european_sample_size", action="store", type="int", dest = "EU_sample_size", default=1006, help="Set EU haploid sample size ; default=1006")
parser.add_option("-a", "--asian_sample_size", action="store", type="int", dest="AS_sample_size", default=1008, help="Set AS haploid sample size ; default=1008")
parser.add_option("-r", "--reference_sample_size", action="store", type="int", dest="AF_sample_size", default=2, help="Set AF haploid sample size ; default=2")

parser.add_option("--migration_AF_B", action="store", type="float", dest="m_AF_B", default=15e-5, help="M[dest,soure]: M[AF,B] (backward time) ; Set African -> ancestral Eurasian migration rate (forward time) ; default=15e-5")
parser.add_option("--migration_B_AF", action="store", type="float", dest="m_B_AF", default=15e-5, help="M[dest,soure]: M[B,AF] (backward time) ; Set African <- ancestral Eurasian migration rate ; default=15e-5")

parser.add_option("--migration_AF_EU", action="store", type="float", dest="m_AF_EU", default=2.5e-5, help="M[dest,soure]: M[AF,EU] (backward time) ; Set African -> European migration rate ; default=2.5e-5")
parser.add_option("--migration_EU_AF", action="store", type="float", dest="m_EU_AF", default=2.5e-5, help="M[dest,soure]: M[EU,AF] (backward time) ; Set African <- European migration rate ; default=2.5e-5")

parser.add_option("--migration_AF_AS", action="store", type="float", dest="m_AF_AS", default=0.78e-5, help="M[dest,soure]: M[AF,AS] (backward time) ; Set African -> Asian migration rate ; default=0.78e-5")
parser.add_option("--migration_AS_AF", action="store", type="float", dest="m_AS_AF", default=0.78e-5, help="M[dest,soure]: M[AS,AF] (backward time) ; Set African <- Asian migration rate ; default=0.78e-5")

parser.add_option("--migration_EU_AS", action="store", type="float", dest="m_EU_AS", default=3.11e-5, help="M[dest,soure]: M[EU,AS] (backward time) ; Set European -> Asian migration rate ; default=3.11e-5")
parser.add_option("--migration_AS_EU", action="store", type="float", dest="m_AS_EU", default=3.11e-5, help="M[dest,soure]: M[AS,EU] (backward time) ; Set European <- Asian migration rate ; default=3.11e-5")


(options, args) = parser.parse_args()



class introgression_object(object):
    def __init__(self, mutations, haplotypes, intervals):
        self.mutations = mutations
        self.haplotypes = haplotypes
        self.intervals = intervals

def introgressed_samples_fn(ts, neanderthal_mrca, neanderthal_samples, segments):
    # Define the samples that carry introgressed segments?
    trees = ts.trees()
    tree = next(trees)
    for left, right in segments:
##        print("\tINTERVAL:", left, right)
        # Skip ahead to the tree that intersects with this segment.
        # NOTE: Regading tree data; The "Records" group consists of four pieces of information:
        # 1) the left and 2) right coordinates of the coalescing interval, 3) the list of child nodes (modern human samples) and 4) the parent node.
        # Each record returned can be accessed via the attributes 'left', 'right', 'node', 'children', 'time' and 'population'
        # A record represents the assignment of a pair of children 'c' to a parent 'u'. This assignment happens at 't' generations in the past within the population with ID 'd'
        while tree.get_interval()[0] < left:
            ## Pass through the potential trees until you get to the one where the start_position matches the start of the defined segment
            tree = next(trees)
            #print(tree.get_interval()[0], left)
        assert tree.get_interval()[0] == left
        start = None
        last_human_leaves = None
##        print(start, last_human_leaves)
        while tree is not None and tree.get_interval()[1] <= right:
            #print( set(tree.leaves(neanderthal_mrca)) )
            #print( neanderthal_samples )
            human_leaves = set(tree.leaves(neanderthal_mrca)) - set(neanderthal_samples)
            #print( human_leaves )
            ## Create a set() of all the human leaves represented in this tree that share an mrca in the Neandertal lineage
            # (i.e. mrca node is in Neand population)
            ## neanderthal_mrca is supplied from the node_map dictionary as the parent node for the tree present in the Neandertal population
            ## In python set() is an unordered collection of unique values;
            # x = [1, 1, 2, 2, 3, 3] --> set(x) --> set([1, 2, 3])
            ## tree.leaves(n) returns an iterator over all the leaves in this tree underneath the specified node (in this case the neanderthal_mrca)
            #print("\t\tTREE:", tree.get_interval(), human_leaves)
            #print("\t\tTREE:", tree.get_interval(), set(tree.leaves(neanderthal_mrca)))
            if start is None:
                last_human_leaves = human_leaves
                start = tree.get_interval()[0]
            elif human_leaves != last_human_leaves:
                end = tree.get_interval()[0]
                yield start, end, last_human_leaves
                start = end
                last_human_leaves = human_leaves
            tree = next(trees, None)
        yield start, right, last_human_leaves


def ILS_samples_fn(ts, neanderthal_mrca, neanderthal_samples, segments):
    # Define the samples that carry introgressed segments?
    trees = ts.trees()
    tree = next(trees)
    for left, right in segments:
##        print("\tINTERVAL:", left, right)
        # Skip ahead to the tree that intersects with this segment.
        # NOTE: Regading tree data; The "Records" group consists of four pieces of information:
        # 1) the left and 2) right coordinates of the coalescing interval, 3) the list of child nodes (modern human samples) and 4) the parent node.
        # Each record returned can be accessed via the attributes 'left', 'right', 'node', 'children', 'time' and 'population'
        # A record represents the assignment of a pair of children 'c' to a parent 'u'. This assignment happens at 't' generations in the past within the population with ID 'd'
        while tree.get_interval()[0] < left:
            ## Pass through the potential trees until you get to the one where the start_position matches the start of the defined segment
            tree = next(trees)
            #print(tree.get_interval()[0], left)
        if set(neanderthal_samples).issubset(set(tree.leaves(neanderthal_mrca))):
            ## Check that we are looking at a coalescent tree that has leaves in both the Neand population and the human population
            # It is likely many trees will have leaves just in the MH lineage, we dont want to pick these
            assert tree.get_interval()[0] == left
            start = None
            last_human_leaves = None
    ##        print(start, last_human_leaves)
            while tree is not None and tree.get_interval()[1] <= right:
                #print( set(tree.leaves(neanderthal_mrca)) )
                #print( neanderthal_samples )
                human_leaves = set(tree.leaves(neanderthal_mrca)) - set(neanderthal_samples)
                #print( human_leaves )
                ## Create a set() of all the human leaves represented in this tree that share an mrca in the Neandertal lineage
                # (i.e. mrca node is in Neand population)
                ## neanderthal_mrca is supplied from the node_map dictionary as the parent node for the tree present in the Neandertal population
                ## In python set() is an unordered collection of unique values;
                # x = [1, 1, 2, 2, 3, 3] --> set(x) --> set([1, 2, 3])
                ## tree.leaves(n) returns an iterator over all the leaves in this tree underneath the specified node (in this case the neanderthal_mrca)
                #print("\t\tTREE:", tree.get_interval(), human_leaves)
                #print("\t\tTREE:", tree.get_interval(), set(tree.leaves(neanderthal_mrca)))
                if start is None:
                    last_human_leaves = human_leaves
                    start = tree.get_interval()[0]
                elif human_leaves != last_human_leaves:
                    end = tree.get_interval()[0]
                    yield start, end, last_human_leaves
                    start = end
                    last_human_leaves = human_leaves
                tree = next(trees, None)
            yield start, right, last_human_leaves



####### SIMULATION RUNS ############

print("Model: ", options.outdir, "Seed: ", options.seed, "Neand_Pulse1: ", options.n1_admix_prop, "Neand_Pulse2: ", options.n2_admix_prop, "Length: ", options.length, file=sys.stderr, sep ='\t')

S_N1 = options.s_n1
S_N2 = options.s_n2

if (options.outdir == "Tenn"):
        #Tenn_demography(S_N1, S_N2, S_AF, S_EU, S_AS, pulses, seed, n1_admix_prop, n2_admix_prop, outdir, t_n1_n2, haplo, length):
        simulation = msprime_demo_models.Tenn_demography(
                        options.s_n1, options.s_n2,
                        options.AF_sample_size, options.EU_sample_size, options.AS_sample_size,
                        options.pulses,
                        options.seed,
                        options.n1_admix_prop, options.n2_admix_prop,
                        options.outdir,
                        options.t_n1_n2,
                        options.t_n1_sample,
                        options.t_n2_sample,
                        options.haplo,
                        options.length,
                        options.m_AF_B, options.m_B_AF,
                        options.m_AF_AS, options.m_AS_AF,
                        options.m_AF_EU, options.m_EU_AF,
                        options.m_EU_AS, options.m_AS_EU)

elif (options.outdir == "Sriram"):
        simulation = msprime_demo_models.Sriram_demography(
                        options.s_n1, options.s_n2,
                        options.AF_sample_size, options.EU_sample_size, options.AS_sample_size,
                        options.pulses,
                        options.seed,
                        options.n1_admix_prop, options.n2_admix_prop,
                        options.outdir,
                        options.t_n1_n2,
                        options.t_n1_sample,
                        options.t_n2_sample,
                        options.haplo,
                        options.length,
                        options.m_AF_B, options.m_B_AF,
                        options.m_AF_AS, options.m_AS_AF,
                        options.m_AF_EU, options.m_EU_AF,
                        options.m_EU_AS, options.m_AS_EU)

elif (options.outdir == "SplitPop"):
        simulation = msprime_demo_models.SplitPop_demography(
                        options.s_n1, options.s_n2,
                        options.AF_sample_size, options.EU_sample_size, options.AS_sample_size,
                        options.pulses,
                        options.seed,
                        options.n1_admix_prop, options.n2_admix_prop,
                        options.outdir,
                        options.t_n1_n2,
                        options.t_n1_sample,
                        options.t_n2_sample,
                        options.haplo,
                        options.length,
                        options.m_AF_B, options.m_B_AF,
                        options.m_AF_AS, options.m_AS_AF,
                        options.m_AF_EU, options.m_EU_AF,
                        options.m_EU_AS, options.m_AS_EU)

elif (options.outdir == "test"):
        simulation = msprime_demo_models.out_of_africa(
                        options.AF_sample_size, options.EU_sample_size, options.AS_sample_size,
                        options.pulses,
                        options.seed,
                        options.n1_admix_prop, options.n2_admix_prop,
                        options.outdir,
                        options.t_n1_n2,
                        options.haplo,
                        options.length,
                        options.m_AF_B, options.m_B_AF,
                        options.m_AF_AS, options.m_AS_AF,
                        options.m_AF_EU, options.m_EU_AF,
                        options.m_EU_AS, options.m_AS_EU)

elif (options.outdir == "Tenn_nomod"):
        simulation = msprime_demo_models.Tenn_demography(
                        options.s_n1, options.s_n2,
                        options.AF_sample_size, options.EU_sample_size, options.AS_sample_size,
                        3,
                        options.seed,
                        options.n1_admix_prop, options.n2_admix_prop,
                        options.outdir,
                        options.t_n1_n2,
                        options.t_n1_sample,
                        options.t_n2_sample,
                        options.haplo,
                        options.length,
                        options.m_AF_B, options.m_B_AF,
                        options.m_AF_AS, options.m_AS_AF,
                        options.m_AF_EU, options.m_EU_AF,
                        options.m_EU_AS, options.m_AS_EU)

elif (options.outdir == "Tenn_pulsed"):
        simulation = msprime_demo_models.Tenn_demography(
                        options.s_n1, options.s_n2,
                        options.AF_sample_size, options.EU_sample_size, options.AS_sample_size,
                        4,
                        options.seed,
                        options.n1_admix_prop, options.n2_admix_prop,
                        options.outdir,
                        options.t_n1_n2,
                        options.t_n1_sample,
                        options.t_n2_sample,
                        options.haplo,
                        options.length,
                        options.m_AF_B, options.m_B_AF,
                        options.m_AF_AS, options.m_AS_AF,
                        options.m_AF_EU, options.m_EU_AF,
                        options.m_EU_AS, options.m_AS_EU)

## Define the sample indices
EU_count = options.EU_sample_size
AS_count = options.AS_sample_size
AF_count = options.AF_sample_size
N_samples  = range(0, S_N1+S_N2)                # range(0,4) --> 0,1,2,3
AF_samples = range(S_N1+S_N2, S_N1+S_N2+AF_count)        # range(4,6) --> 4,5
non_trgt = len(N_samples + AF_samples)                # 6
nonAfr_samples = range(non_trgt, (EU_count + AS_count + non_trgt))   # range(6, 3046 + 6) ; 3046 <-- 1006 + 2040
modHum_samples = range( (S_N1+S_N2), (AF_count + EU_count + AS_count + (S_N1+S_N2)) )
EU_samples = range(non_trgt, (EU_count + non_trgt))
AS_samples = range((EU_count + non_trgt), (EU_count+AS_count+non_trgt))
Chimp_samples = range((EU_count+AS_count+non_trgt), (EU_count+AS_count+non_trgt)+2)
Deni_samples = range((EU_count+AS_count+non_trgt) + 2, (EU_count+AS_count+non_trgt) + 2 + 2)

if (options.pop == "EAS"):
        human_samples = AS_samples
elif (options.pop == "EUR"):
        human_samples = EU_samples
elif (options.pop == "nonAfr"):
        human_samples = nonAfr_samples
elif (options.pop == "AFR"):
        human_samples = AF_samples
elif (options.pop == "modHum"):
        human_samples = modHum_samples


###### GET HAPLOTYPES FROM SIMULATED TREES #######

if (options.haplo == "haplo"):
        # Create a .bed file to write to for the simulation
        haplo_outfile = sys.stdout
        #haplo_outfile = gzip.open( options.outdir+'_'+options.pop+'_'+str(options.seed)+'_n1_'+str(options.n1_admix_prop)+'_n2_'+str(options.n2_admix_prop)+'.bed.merged.gz' , 'wb' )
        haplo_entry_list=[]
##  FOR EACH SIMULATED CHROMOSOME, PRINT ALL THE INTROGRESSED HAPLOTYPES BELONGING TO THE SPECIFIED NON-AFR POPULATION IN BED FORMAT  ###
        for t, ts in enumerate(simulation):
            ## t is the tree ID, ts is a given tree from the simulation
            node_map = collections.defaultdict(list)
            ## 'collections' = pythons high-performance container types, 'defaultdict()' is one such container type
            # defaultdict() = dict subclass that calls a factory function to supply missing values
            # Using 'list' as the default_factory, it is easy to group a sequence of key-value pairs into a dictionary of lists:
            ## Defines "node_map" as an empty dictionary of lists ; [('n1', [0, 10]), ('n2', [11, 14]), ('n3', [15, 20])]
            for record in ts.records():
                ## for a single record from all the records of a given tree ts
                # We are interested in coalescence events that occured in the Neandertal population
                if record.population <= 1:
                    ## If the population of the tree is Neand (0 or 1)
                    # Since our Neanderthal sample is only two, we can easily exclude events that just concern Neanderthals
                    if set(record.children) != N_samples:
                        ## If the child/leaf nodes for the record are not Neandertal samples
                        # add the node ID and the position of the segments to node_map
#                        print("record:", record)
                        node_map[record.node].append((record.left, record.right))
            for neanderthal_mrca, segments in node_map.items():
                ## >>> node_map.items()
                # [('n1', [0, 10]), ('n2', [11, 14]), ('n3', [15, 20])]
                # where the node is the neanderthal_mrca (e.g. 'n1'), and the segments is the tree interval (e.g. [0,10])
#                print(neanderthal_mrca, "->", segments)
                iterator = introgressed_samples_fn(ts, neanderthal_mrca, N_samples, segments)
                ## Run the introgressed_samples function,
                # using the given tree, the defined Neandertal_mrca/node, the defined tree interval/'segments', and the defined Neand_samples
                for left, right, samples in iterator:
#                        print("\t", left, right, samples, file=sys.stderr)
                        for s in samples:
                                if s in human_samples:
                                        if int(math.ceil(left)) < int(math.ceil(right)):
#                                                print(s, int(math.ceil(left)), int(math.ceil(right)), s)
                                                haplo_entry=str(s)+'\t'+str(int(math.ceil(left)))+'\t'+str(int(math.ceil(right)))+'\t'+str(s)
                                                haplo_entry_list.append(haplo_entry+'\n')
                                                # Only print nonAfr samples
                                                ## NOTE: After printing to a bed file, still need to sort and merge the bedfile.
        ## Join together the list of introgressed haplotypes into a string, and convert this to a BED file using pybedtools
        ## Then perform the sort() and merge() functions in python on the BEDfile object
        haplo_entry_string = ''.join(haplo_entry_list)
        pybedtools.set_tempdir('/scratch/tmp/abwolf/msprime/')
        #print(pybedtools.get_tempdir(), file=sys.stderr)
        BEDFILE = pybedtools.BedTool(haplo_entry_string, from_string=True)
        BEDFILE_SORTED_MERGED = pybedtools.BedTool.sort(BEDFILE).merge()

#        print(BEDFILE_SORTED_MERGED.head(),file=sys.stdout)
#        print(BEDFILE_SORTED_MERGED)

        ## Read the BEDfile line by line, add in the chr#, reorder so that the columns are: chr, strt, end, ind
        ## write to haplo_outfile in gzip
        for bed_line in BEDFILE_SORTED_MERGED:
                new_bed_line = str(bed_line).strip() + '\t' + str(options.seed)
                new_bed_line = new_bed_line.split('\t')
                new_bed_line[0], new_bed_line[3] = new_bed_line[3], new_bed_line[0]
#                print(new_bed_line, file=sys.stdout)
#                print('\t'.join(new_bed_line))
                haplo_outfile.write('\t'.join(new_bed_line)+'\n')

        pybedtools.cleanup(verbose=False)
        haplo_outfile.close()



elif (options.haplo == "ILS"):
        # Create a .bed file to write to for the simulation
        ILS_outfile = sys.stdout
        #haplo_outfile = gzip.open( options.outdir+'_'+options.pop+'_'+str(options.seed)+'_n1_'+str(options.n1_admix_prop)+'_n2_'+str(options.n2_admix_prop)+'.bed.merged.gz' , 'wb' )
        haplo_entry_list=[]
##  FOR EACH SIMULATED CHROMOSOME, PRINT ALL THE INTROGRESSED HAPLOTYPES BELONGING TO THE SPECIFIED NON-AFR POPULATION IN BED FORMAT  ###
        for t, tree_sequence in enumerate(simulation):
            ## t is the tree ID, ts is a given tree from the simulation
            node_map = collections.defaultdict(list)
            ## 'collections' = pythons high-performance container types, 'defaultdict()' is one such container type
            # defaultdict() = dict subclass that calls a factory function to supply missing values
            # Using 'list' as the default_factory, it is easy to group a sequence of key-value pairs into a dictionary of lists:
            ## Defines "node_map" as an empty dictionary of lists ; [('n1', [0, 10]), ('n2', [11, 14]), ('n3', [15, 20])]
            for record in tree_sequence.records():
                ## for a single record from all the records of a given tree ts
                # We are interested in coalescence events that occured in the Neandertal population
#                if record.population <= 1:
                if record.population == 2:
                    ## Modified to look for ILS segments, where colescent event occurs in the ancestral Neand-EUR_ASN population
                    if record.time >= 28000 and record.time < 50000:
                        ## To ensure we pick up on only some ILS segements
                        # limit the time range for coalescent events...
                        # If you allow coalescent events to occur at any point from T_N_MH to T_MH_CH, the whole genome will called as ILS
                        node_map[record.node].append((record.left, record.right))
            for neanderthal_mrca, segments in node_map.items():
                ## >>> node_map.items()
                # [('n1', [0, 10]), ('n2', [11, 14]), ('n3', [15, 20])]
                # where the node is the neanderthal_mrca (e.g. 'n1'), and the segments is the tree interval (e.g. [0,10])
#                print(neanderthal_mrca, "->", segments)
                iterator = ILS_samples_fn(tree_sequence, neanderthal_mrca, N_samples, segments)
                ## Run the introgressed_samples function,
                # using the given tree, the defined Neandertal_mrca/node, the defined tree interval/'segments', and the defined Neand_samples
                for left, right, samples in iterator:
#                        print("\t", left, right, samples, file=sys.stderr)
                        for s in samples:
                                if s in human_samples:
                                        if int(math.ceil(left)) < int(math.ceil(right)):
#                                                print(s, int(math.ceil(left)), int(math.ceil(right)), s)
                                                haplo_entry=str(s)+'\t'+str(int(math.ceil(left)))+'\t'+str(int(math.ceil(right)))+'\t'+str(s)
                                                haplo_entry_list.append(haplo_entry+'\n')
                                                # Only print nonAfr samples
                                                ## NOTE: After printing to a bed file, still need to sort and merge the bedfile.
        ## Join together the list of introgressed haplotypes into a string, and convert this to a BED file using pybedtools
        ## Then perform the sort() and merge() functions in python on the BEDfile object
        haplo_entry_string = ''.join(haplo_entry_list)
        pybedtools.set_tempdir('/scratch/tmp/abwolf/msprime/')
        #print(pybedtools.get_tempdir(), file=sys.stderr)
        BEDFILE = pybedtools.BedTool(haplo_entry_string, from_string=True)
        BEDFILE_SORTED_MERGED = pybedtools.BedTool.sort(BEDFILE).merge()

#        print(BEDFILE_SORTED_MERGED.head(),file=sys.stdout)
#        print(BEDFILE_SORTED_MERGED)

        ## Read the BEDfile line by line, add in the chr#, reorder so that the columns are: chr, strt, end, ind
        ## write to haplo_outfile in gzip
        for bed_line in BEDFILE_SORTED_MERGED:
                new_bed_line = str(bed_line).strip() + '\t' + str(options.seed)
                new_bed_line = new_bed_line.split('\t')
                new_bed_line[0], new_bed_line[3] = new_bed_line[3], new_bed_line[0]
#                print(new_bed_line, file=sys.stdout)
#                print('\t'.join(new_bed_line))
                ILS_outfile.write('\t'.join(new_bed_line)+'\n')

        pybedtools.cleanup(verbose=False)
        ILS_outfile.close()



##########    GET EIGENSTRATGENO FILES AND SNP FILES    #########

elif (options.haplo == "F4Dstat"):
            outdir = options.outdir
            n1_admix_prop = options.n1_admix_prop
            n2_admix_prop = options.n2_admix_prop
            t_n1_n2 = options.t_n1_n2
            length = options.length
            seed = options.seed
            m_AF_B = options.m_AF_B
            m_B_AF = options.m_B_AF
            m_AF_EU = options.m_AF_EU
            m_EU_AF = options.m_EU_AF

            tag = 'n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_mAfB_'+str(m_AF_B)+'_mBAf_'+str(m_B_AF)+'_mAfEu_'+str(m_AF_EU)+'_mEuAf_'+str(m_EU_AF)+'_'+str(seed)

            # geno_outfile = gzip.open(outdir+'.eigenstratgeno.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed)+'.gz', 'wb')  ##  sim.eigenstratgeno.n1_0.01_n2_0.05_200
            # snp_outfile = gzip.open(outdir+'.snp.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed)+'.gz', 'wb')
            # ind_outfile = gzip.open(outdir+'.ind.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed)+'.gz', 'wb')

            geno_outfile = gzip.open(outdir+'.eigenstratgeno.'+tag+'.gz', 'wb')  ##  sim.eigenstratgeno.n1_0.01_n2_0.05_200
            snp_outfile = gzip.open(outdir+'.snp.'+tag+'.gz', 'wb')
            ind_outfile = gzip.open(outdir+'.ind.'+tag+'.gz', 'wb')

            rs_num = 0


            for t, tree_sequence in enumerate(simulation): # enumerate creates a list of tuples, where each chrom is a tuple with an index and the tree e.g. [(1,Tree1), (2,Tree2)...]
                ###  WRITE THE .IND FILE  ###
                if t<1:                                                        ## Write .ind file based on output of Tree1
                        for i in range(tree_sequence.get_sample_size()):        ## Get the sample size from the tree
                                if tree_sequence.get_population(i)==0:          ## For each individual in Tree1, get the corresponding population
                                        pop='Neand1'
                                elif tree_sequence.get_population(i)==1:
                                        pop='Neand2'
                                elif tree_sequence.get_population(i)==2:
                                        pop='AFR'
                                elif tree_sequence.get_population(i)==3:
                                        pop='EUR'
                                elif tree_sequence.get_population(i)==4:
                                        pop='ASN'
                                elif tree_sequence.get_population(i)==5:
                                        pop='Chimp'
                                elif tree_sequence.get_population(i)==6:
                                        pop='Deni'
                                elif tree_sequence.get_population(i)==7:
                                        pop='Split'
                                ind_entry = 'Sample_'+str(i)+'\t'+'U'+'\t'+pop  ## Write, as a string, the sample entry using the population read from Tree1
                                ind_outfile.write(ind_entry+'\n')
                ###  WRITE THE .EIGENSTRATGENO AND .SNP FILES  ###
                chr_num = t+1
                for variant in tree_sequence.variants(as_bytes=True):
                    rs_num+=1
                    geno_outfile.write(variant.genotypes+'\n') # write genotypes to .eigenstratgeno file
                    line = str('rs'+str(rs_num)+'\t'+str(chr_num)+'\t'+str(variant.position/length)+'\t'+str(int(variant.position))+'\t'+'A'+'\t'+'T'+'\n')
                    snp_outfile.write(line) # write snp_allele info to .snp file


            # parF4_outfile = gzip.open('parfile.F4stat.'+outdir+'.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed)+'.gz', 'wb')
            # parF4_outfile.write('genotypename: '+outdir+'.eigenstratgeno.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed)+'\n')
            # parF4_outfile.write('snpname: '+outdir+'.snp.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed)+'\n')
            # parF4_outfile.write('indivname: '+outdir+'.ind.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed)+'\n')
            # parF4_outfile.write('popfilename: sim.popfile_F4stat'+'\n')

            parF4_outfile = gzip.open('parfile.F4stat.'+outdir+'.'+tag+'.gz', 'wb')
            parF4_outfile.write('genotypename: '+outdir+'.eigenstratgeno.'+tag+'\n')
            parF4_outfile.write('snpname: '+outdir+'.snp.'+tag+'\n')
            parF4_outfile.write('indivname: '+outdir+'.ind.'+tag+'\n')
            parF4_outfile.write('popfilename: sim.popfile_F4stat'+'\n')
#
#            parD_outfile = open('parfile.Dstat.'+outdir+'.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed), 'w')
#            parD_outfile.write('genotypename: '+outdir+'.eigenstratgeno.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed)+'\n')
#            parD_outfile.write('snpname: '+outdir+'.snp.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed)+'\n')
#            parD_outfile.write('indivname: '+outdir+'.ind.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed)+'\n')
#            parD_outfile.write('popfilename: sim.popfile_Dstat'+'\n')


            geno_outfile.close()
            snp_outfile.close()
            ind_outfile.close()
            parF4_outfile.close()
#            parD_outfile.close()


##########    WRITE VCF FILE FORMAT FOR S* CALCULATIONS , ALONG WITH POPULATION FILE   #########

elif (options.haplo == "vcf"):
            outdir = options.outdir
            n1_admix_prop = options.n1_admix_prop
            n2_admix_prop = options.n2_admix_prop
            pop = options.pop
            t_n1_n2 = options.t_n1_n2
            length = options.length
            seed = options.seed

            for t, tree_sequence in enumerate(simulation):

                    if os.path.isfile(outdir+'.popfile') == True:        # If a population file does not exist yet, write one
                        print('popfile already exists', file=sys.stderr)
                    elif os.path.isfile(outdir+'.popfile') != True:
                        pop_outfile = open(outdir+'.popfile', 'w') ## Tenn.popfile.gz
                        pop_outfile.write('samp'+'\t'+'pop'+'\t'+'super_pop'+'\n')                #write header to pop_outfile

                        for i in range(tree_sequence.get_sample_size()):
                            if tree_sequence.get_population(i)==0:          ## For each individual in Tree1, get the corresponding population
                                    s_pop='Neand1'
                            elif tree_sequence.get_population(i)==1:
                                    s_pop='Neand2'
                            elif tree_sequence.get_population(i)==2:
                                    s_pop='AFR'
                            elif tree_sequence.get_population(i)==3:
                                    s_pop='EUR'
                            elif tree_sequence.get_population(i)==4:
                                    s_pop='ASN'
                            elif tree_sequence.get_population(i)==5:
                                    s_pop='Chimp'
                            elif tree_sequence.get_population(i)==6:
                                    s_pop='Deni'
                            elif tree_sequence.get_population(i)==7:
                                    s_pop='Split'
                            if i%2==0:
                                    ind_entry = 'msp_'+str(i/2)+'\t'+s_pop+'\t'+s_pop  ## Write, as a string, the sample entry using the population read from Tree1
                                    pop_outfile.write(ind_entry+'\n')
                        pop_outfile.close()

                    vcf_outfile = sys.stdout
                    #vcf_outfile = gzip.open(outdir+'_'+pop+'_'+str(seed)+'_n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'.vcf.gz', 'wb')  ## Tenn_nonAfr_100_n1_0.01_n2_0.01.vcf.gz
                    tree_sequence.write_vcf(vcf_outfile,2)
                    #tree_sequence.write_vcf(sys.stdout,2)

            vcf_outfile.close()



print('fin',file=sys.stderr)
