from __future__ import print_function
import msprime
import math
import numpy as np
import pybedtools
import random
import scipy
import subprocess
import glob
import gzip
from scipy import spatial
from optparse import OptionParser

# usage = str("usage: python <script> -p nonAfr -o Tenn -s 1 -i 2 -n 0.02 -d 0.0 -t 350 -c F4Dstat -l 1e6 -e 1006 -a 1008\n"+
#             "WARNING: options 'haplo' and 'vcf' now write to stdout")
# parser = OptionParser(usage=usage)
# parser.add_option("-p", "--population", action = "store", type = "string", dest = "pop", default="nonAfr", help="call introg_haplotypes in EUR, EAS, or all nonAfr; default=nonAfr")
# parser.add_option("-o", "--outdir", action = "store", type = "string", dest = "outdir", default="Tenn", help="output directory name, best to include model name; default=Tenn")
# parser.add_option("-s", "--seed", action = "store", type = "int", dest = "seed", default=1, help="Set random seed for replicate chromosomes; default=1")
# parser.add_option("-i", "--introgress_pulses", action = "store", type = "int", dest = "pulses", default=2, help="Set number of introgression pulses; default=2")
# parser.add_option("-n", "--neand1_admixture_proportion", action = "store", type = "float", dest = "n1_admix_prop", default=0.02, help="Set N1 admixture proportion; default=0.02")
# parser.add_option("-d", "--deni_admixture_propotion", action = "store", type = "float", dest = "n2_admix_prop", default=0.0, help="Set N2 admixture proportion; default=0.0")
# parser.add_option("-t", "--time_N1_N2_split", action = "store", type = "int", dest = "t_n1_n2", default=350, help="Set N1 N2 split time in kya; default=350")
# parser.add_option("-c", "--calls_or_stats", action="store", type = "string", dest = "haplo", default="F4Dstat", help="Pick haplotype calls ('haplo') or 'F4Dstat' or 'vcf' output, or 'debug'; default=F4Dstat")
# parser.add_option("-l", "--length_chrom", action="store", type= "float", dest = "length", default=1e6, help="Define length of simulated chromosome; default=1e6")
# parser.add_option("-e", "--european_sample_size", action="store", type="int", dest = "EU_sample_size", default=1006, help="Set EU haploid sample size ; default=1006")
# parser.add_option("-a", "--asian_sample_size", action="store", type="int", dest="AS_sample_size", default=1008, help="Set AS haploid sample size ; default=1008")
# parser.add_option("-r", "--reference_sample_size", action="store", type="int", dest="AF_sample_size", default=2, help="Set AF haploid sample size ; default=2")
#
# parser.add_option("--migration_AF_B", action="store", type="float", dest="m_AF_B", default=15e-5, help="Set African -> ancestral Eurasian migration rate ; default=15e-5")
# parser.add_option("--migration_B_AF", action="store", type="float", dest="m_B_AF", default=15e-5, help="Set African <- ancestral Eurasian migration rate ; default=15e-5")
#
# parser.add_option("--migration_AF_EU", action="store", type="float", dest="m_AF_EU", default=2.5e-5, help="Set African -> European migration rate ; default=2.5e-5")
# parser.add_option("--migration_EU_AF", action="store", type="float", dest="m_EU_AF", default=2.5e-5, help="Set African <- European migration rate ; default=2.5e-5")
#
# parser.add_option("--migration_AF_AS", action="store", type="float", dest="m_AF_AS", default=0.78e-5, help="Set African -> Asian migration rate ; default=0.78e-5")
# parser.add_option("--migration_AS_AF", action="store", type="float", dest="m_AS_AF", default=0.78e-5, help="Set African <- Asian migration rate ; default=0.78e-5")
#
# parser.add_option("--migration_EU_AS", action="store", type="float", dest="m_EU_AS", default=3.11e-5, help="Set European -> Asian migration rate ; default=3.11e-5")
# parser.add_option("--migration_AS_EU", action="store", type="float", dest="m_AS_EU", default=3.11e-5, help="Set European <- Asian migration rate ; default=3.11e-5")


#length = 1e6        # simulate 1Mb chromosomes for faster F4 and D calculation
recombination_rate = 1e-8          # 2e-8 original from Rajiv, not sure where it is from.
#mutation_rate = 1e-8                 # see Shendure and Akey, 2015
#mutation_rate = 2.5e-8         # used in Vernot 2014 and 2016 ; matches 0migration and CWZ num_s_star_snps, matches 0migration s_star scores
mutation_rate = 1.2e-8         # value best matches CWZ data n_region_ind_snps and s_star scores

def Tenn_demography(S_N1, S_N2, S_AF, S_EU, S_AS,
                pulses, seed,
                n1_admix_prop, n2_admix_prop,
                outdir, t_n1_n2, t_n1_sample, t_n2_sample,
                haplo, length,
                m_AF_B, m_B_AF,
                m_AF_AS, m_AS_AF,
                m_AF_EU, m_EU_AF,
                m_EU_AS, m_AS_EU):
    N_A = 7310 # effective population size used for scaling
    N_N1 = 1000 # Altai/Vindija lineage effective population size
    N_N2 = 1000 # Eastern Neandertal effective population size
    N_B = 1861 # Population size after EUR + ASN join
    N_B_BN = 100 # Population size of EUR_ASN bottleneck (following initial admixture event)
    N_AF0 = 14474
    N_EU0 = 1032
    N_AS0 = 554
    N_CH = 1000
    N_DE = 1000
    # Times are provided in years, so we convert into generations.
    generation_time = 25
    T_MH_CH = 7e6 / generation_time # Chimp_Human lineage split
    T_MH_N = 700e3 / generation_time # Neanderthal / modern human split
    T_DE_N = 500e3 / generation_time # Denisovan / Neandertal split time
    T_N1_N2 = ((t_n1_n2) * 1e3) / generation_time # Altai/Vindija and Eastern Neandertal popualation split, default set to 145kya
    T_N1_SAMPLE = ((t_n1_sample * 1e3) / generation_time)
    T_N2_SAMPLE = ((t_n2_sample * 1e3) / generation_time)
    T_AF = 148e3 / generation_time # African population expansion
    T_B = 100e3 / generation_time # out of Africa
    T_PULSE1 = 55e3 / generation_time # Neandertal introgression pulse 1
    T_B_BN = 45e3 / generation_time # A 3rd Bottleneck following shortly after initial admixture event
    T_EU_AS = 23e3 / generation_time # European / East Asian split
    T_PULSE2 = 18e3 / generation_time # Neandertal introgression pulse 2
    T_ACL_GRW = 5.1e3 / generation_time # time of rapid population growth
    # We need to work out the starting (diploid) population sizes based on
    # the growth rates provided for these two populations
    r_EU1 = 0.00307
    r_AS1 = 0.0048
    #N_EU1 = N_EU0 / math.exp(-r_EU1 * (T_EU_AS-T_ACL_GRW))
    #N_AS1 = N_AS0 / math.exp(-r_AS1 * (T_EU_AS-T_ACL_GRW))
    r_EU2 = 0.0195        # rapid population growth at 5.1kya
    r_AS2 = 0.025        # rapid population growth at 5.1kya
    r_AF2 = 0.0166        # rapid popualtion growth at 5.1kya
    #N_EU2 = N_EU1 / math.exp(-r_EU2 * T_ACL_GRW)        # should reach ~512k
    #N_AS2 = N_AS1 / math.exp(-r_AS2 * T_ACL_GRW)        # should reach ~640k
    #N_AF2 = N_AF / math.exp(-r_AF2 * T_ACL_GRW)        # should reach ~424k
    N_EU = 512e3
    N_AS = 640e3
    N_AF = 424e3
    # Migration rates during the various epochs.
    #######
    #m_AF_B = 15e-5
    #m_AF_EU = 2.5e-5
    #m_AF_AS = 0.78e-5
    #m_EU_AS = 3.11e-5
    #######
    m_AF_B = m_AF_B
    m_AF_EU = m_AF_EU
    m_AF_AS = m_AF_AS
    m_EU_AS = m_EU_AS

    m_B_AF = m_B_AF
    m_EU_AF = m_EU_AF
    m_AS_AF = m_AS_AF
    m_AS_EU = m_AS_EU


#     m_AF_B_s = (float(T_B_BN - T_EU_AS) * m_AF_B)
#     m_AF_EU_s = m_AF_EU
#     m_AF_AS_s = m_AF_AS
#     m_EU_AS_s = m_EU_AS
#
#     m_B_AF_s = (float(T_B_BN - T_EU_AS) * m_B_AF)
#     m_EU_AF_s = m_EU_AF
#     m_AS_AF_s = m_AS_AF
#     m_AS_EU_s = m_AS_EU

    ####
    m_N1_N2 = 0
    m_N1_AF = 0
    m_N1_EU = 0
    m_N1_AS = 0
    m_N2_AF = 0
    m_N2_EU = 0
    m_N2_AS = 0
    m_N1_CH = 0
    m_N2_CH = 0
    m_AF_CH = 0
    m_EU_CH = 0
    m_AS_CH = 0
    m_DE_N1 = 0
    m_DE_N2 = 0
    m_DE_AF = 0
    m_DE_EU = 0
    m_DE_AS = 0
    m_DE_CH = 0
    # Set Percent introgression level
    m_PULSE1 = n1_admix_prop # an instantaneous movement of population from EU_AS to Neand1 (backward in time)
    m_PULSE2 = n2_admix_prop
    #print(m_PULSE1)

    # Population IDs correspond to their indexes in the population
    # configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB
    # initially.
    population_configurations = [
        msprime.PopulationConfiguration(
        initial_size = N_N1), # population 0
        msprime.PopulationConfiguration(
        initial_size = N_N2), # population 1
        msprime.PopulationConfiguration(
        initial_size = N_AF, growth_rate = r_AF2), # population 2
        msprime.PopulationConfiguration(
        initial_size = N_EU, growth_rate = r_EU2), # population 3
        msprime.PopulationConfiguration(
        initial_size = N_AS, growth_rate = r_AS2), # population 4
        msprime.PopulationConfiguration(
        initial_size = N_CH), #popultaion 5
        msprime.PopulationConfiguration(
        initial_size = N_DE)  #population 6
    ]
    # specify desired sample
    N1_sample = [msprime.Sample(population = 0, time = T_N1_SAMPLE)] * S_N1
    N2_sample = [msprime.Sample(population = 1, time = T_N2_SAMPLE)] * S_N2
    AF_sample = [msprime.Sample(population = 2, time = 0)] * S_AF
    EU_sample = [msprime.Sample(population = 3, time = 0)] * S_EU
    AS_sample = [msprime.Sample(population = 4, time = 0)] * S_AS
    CH_sample = [msprime.Sample(population = 5, time = 0)] * 2
    DE_sample = [msprime.Sample(population = 6, time = 50e3/ generation_time)] * 2
    samples = N1_sample + N2_sample + AF_sample + EU_sample + AS_sample + CH_sample + DE_sample
    # set up the initial migration matrix
    ## This migration matrix is the transposition of the migration matrix printed during debugger
    ## however, the format here is the format you should reference and the matrix indecese
    ## you should use when changing the migration parameters
    ## i.e. though m_AF_EU is at position (2,3) in the debugger matix, here and in the code it should be refered to at (3,2)
    migration_matrix = [
        [      0, m_N1_N2, m_N1_AF, m_N1_EU, m_N1_AS, m_N1_CH, m_DE_N1],
        [m_N1_N2,       0, m_N2_AF, m_N2_EU, m_N2_AS, m_N2_CH, m_DE_N2],
        [m_N1_AF, m_N2_AF,       0, m_EU_AF, m_AS_AF, m_AF_CH, m_DE_AF],
        [m_N1_EU, m_N2_EU, m_AF_EU,       0, m_AS_EU, m_EU_CH, m_DE_EU],
        [m_N1_AS, m_N2_AS, m_AF_AS, m_EU_AS,       0, m_AS_CH, m_DE_AS],
        [m_N1_CH, m_N2_CH, m_AF_CH, m_EU_CH, m_AS_CH,       0, m_DE_CH],
        [m_DE_N1, m_DE_N2, m_DE_AF, m_DE_EU, m_DE_AS, m_DE_CH,       0]
    ]

    one_pulse = [
        # msprime.PopulationParametersChange(time = T_ACL_GRW, initial_size = N_AF0, growth_rate = 0, population_id = 2),        # stop rapid population growth in AF
        # msprime.PopulationParametersChange(time = T_ACL_GRW, growth_rate = r_EU1, population_id = 3),        # stop rapid population growth in EU
        # msprime.PopulationParametersChange(time = T_ACL_GRW, growth_rate = r_AS1, population_id = 4),        # stop rapid population growth in AS
        # msprime.PopulationParametersChange(time = T_PULSE2-1, initial_size = N_AS0, growth_rate = 0, population_id = 4), # set AS popsize to AS0
        # msprime.PopulationParametersChange(time = T_EU_AS-1, initial_size = N_EU0, growth_rate = 0, population_id = 3),        # set EU popsize at EU0
        # msprime.PopulationParametersChange(time = T_EU_AS-1, initial_size = N_AS0, growth_rate = 0, population_id = 4),         # set AS popsize to AS0
        # msprime.MassMigration(time = T_EU_AS, source = 4, destination = 3, proportion = 1.0), # AS merges into EU, now termed "B"
        # msprime.MigrationRateChange(time = T_EU_AS, rate = 0), # set all migration rates to zero
        # msprime.MigrationRateChange(time = T_EU_AS, rate = m_AF_B, matrix_index = (2, 3)), # migration between "B" and Africa begins
        # msprime.MigrationRateChange(time = T_EU_AS, rate = m_AF_B, matrix_index = (3, 2)),
        # msprime.PopulationParametersChange(time = T_EU_AS, initial_size = N_B, growth_rate = 0, population_id = 3), # set parameters of population "B"
        # #msprime.PopulationParametersChange(time = T_B_BN, initial_size = N_B_BN, growth_rate = 0, population_id = 3),  # population bottleneck of B begins
        # #msprime.PopulationParametersChange(time = T_B_BN + 20, initial_size = N_B, growth_rate = 0, population_id = 3), # population bottleneck ends shortly before initial admixture
        # msprime.MassMigration(time = T_PULSE1, source = 3, destination = 0, proportion = m_PULSE1 ), # Neand1 to EUR_EAS pulse of introgression
        # msprime.MassMigration(time = T_B, source = 3, destination = 2, proportion = 1.0), # Population B merges into Africa at T_B
        # msprime.MigrationRateChange(time = T_B, rate = 0), # set all migration rates to zero
        # msprime.PopulationParametersChange(time = T_AF, initial_size = N_A, population_id = 2), # set parameters of ancestral modern human population
        # msprime.MassMigration(time = T_N1_N2, source = 1, destination = 0, proportion = 1.0), # N_2 merges with N_1 at T_N1_N2
        # msprime.PopulationParametersChange(time = T_N1_N2, initial_size = N_N1, population_id = 0), # set parameters of ancestral Neandertal population
        # msprime.MassMigration(time = T_DE_N, source = 6, destination = 0, proportion = 1.0), # DE merges with N1
        # msprime.PopulationParametersChange(time = T_DE_N, initial_size = N_N1, population_id = 0),
        # msprime.MassMigration(time = T_MH_N, source = 0, destination = 2, proportion = 1.0), # Neandertals merge into modern human lineage at time T_MH_N
        # msprime.PopulationParametersChange(time = T_MH_N, initial_size = N_A, population_id = 2), # set parameters of ancetral hominin population
        # msprime.MassMigration(time = T_MH_CH, source = 5, destination = 2, proportion = 1.0), # Chimp lineage merges into ancestral hominin population at time T_MH_Ch
        # msprime.PopulationParametersChange(time = T_MH_CH, initial_size = N_A, population_id = 2) # set parameters of ancestral hominin population
    ]

    two_pulse = [
        msprime.PopulationParametersChange(time = T_ACL_GRW, initial_size = N_AF0, growth_rate = 0, population_id = 2),        # stop rapid population growth in AF
        msprime.PopulationParametersChange(time = T_ACL_GRW, growth_rate = r_EU1, population_id = 3),        # stop rapid population growth in EU
        msprime.PopulationParametersChange(time = T_ACL_GRW, growth_rate = r_AS1, population_id = 4),        # stop rapid population growth in AS

        msprime.PopulationParametersChange(time = T_PULSE2-1, initial_size = N_AS0, growth_rate = 0, population_id = 4), # set AS popsize to AS0
        msprime.MassMigration(time = T_PULSE2, source = 4, destination = 0, proportion = m_PULSE2), # Neand1 to EAS pulse of introgression

        msprime.PopulationParametersChange(time = T_EU_AS-1, initial_size = N_EU0, growth_rate = 0, population_id = 3),        # set EU popsize at EU0
        msprime.PopulationParametersChange(time = T_EU_AS-1, initial_size = N_AS0, growth_rate = 0, population_id = 4),         # set AS popsize to AS0
        msprime.MassMigration(time = T_EU_AS, source = 4, destination = 3, proportion = 1.0), # AS merges into EU, now termed "B"
        msprime.MigrationRateChange(time = T_EU_AS, rate = 0), # set all migration rates to zero
        msprime.MigrationRateChange(time = T_EU_AS, rate = m_AF_B, matrix_index = (3, 2)), # Change rate m_AF_EU to m_AF_B, m(row, column) ; m(source, dest) backward in time
        msprime.MigrationRateChange(time = T_EU_AS, rate = m_B_AF, matrix_index = (2, 3)),
        msprime.PopulationParametersChange(time = T_EU_AS, initial_size = N_B, growth_rate = 0, population_id = 3), # set parameters of population "B"
        #msprime.PopulationParametersChange(time = T_B_BN, initial_size = N_B_BN, growth_rate = 0, population_id = 3), # population bottleneck of B begins
        #msprime.PopulationParametersChange(time = T_B_BN + 20, initial_size = N_B, growth_rate = 0, population_id = 3), # population bottleneck of B ends shortly before initial admixture
        msprime.MassMigration(time = T_PULSE1, source = 3, destination = 0, proportion = m_PULSE1), # Neand1 to EUR_EAS pulse of introgression
        msprime.MassMigration(time = T_B, source = 3, destination = 2, proportion = 1.0), # Population B merges into Africa at T_B
        msprime.MigrationRateChange(time = T_B, rate = 0), # set all migration rates to zero
        msprime.MassMigration(time = T_N1_N2, source = 1, destination = 0, proportion = 1.0), # N_2 merges with N_1 at T_N1_N2
        msprime.PopulationParametersChange(time = T_N1_N2, initial_size = N_N1, population_id = 0), # set parameters of ancestral Neandertal population
        msprime.PopulationParametersChange(time = T_AF, initial_size = N_A, population_id = 2), # set parameters of ancestral modern human population
        msprime.MassMigration(time = T_DE_N, source = 6, destination = 0, proportion = 1.0), # DE merges with N1
        msprime.PopulationParametersChange(time = T_DE_N, initial_size = N_N1, population_id = 0),
        msprime.MassMigration(time = T_MH_N, source = 0, destination = 2, proportion = 1.0), # Neandertals merge into modern human lineage at time T_MH_N
        msprime.PopulationParametersChange(time = T_MH_N, initial_size = N_A, population_id = 2), # set parameters of ancetral hominin population
        msprime.MassMigration(time = T_MH_CH, source = 5, destination = 2, proportion = 1.0), # Chimp lineage merges into ancestral hominin population at time T_MH_CH
        msprime.PopulationParametersChange(time = T_MH_CH, initial_size = N_A, population_id = 2) # set parameters of ancestral hominin population
    ]

    no_modern_migration = [
            msprime.MigrationRateChange(time = 0, rate = 0),    # turn off all migration from t=0 until T_ACL_GRW (5kya)

            msprime.PopulationParametersChange(time = T_ACL_GRW, initial_size = N_AF0, growth_rate = 0, population_id = 2),        # stop rapid population growth in AF
            msprime.PopulationParametersChange(time = T_ACL_GRW, growth_rate = r_EU1, population_id = 3),        # stop rapid population growth in EU
            msprime.PopulationParametersChange(time = T_ACL_GRW, growth_rate = r_AS1, population_id = 4),        # stop rapid population growth in AS

            msprime.MigrationRateChange(time = T_ACL_GRW, rate = m_AF_EU, matrix_index = (3,2)),    #start modern migration parameters at 5kya
            msprime.MigrationRateChange(time = T_ACL_GRW, rate = m_EU_AF, matrix_index = (2,3)),
            msprime.MigrationRateChange(time = T_ACL_GRW, rate = m_AF_AS, matrix_index = (4,2)),
            msprime.MigrationRateChange(time = T_ACL_GRW, rate = m_AS_AF, matrix_index = (2,4)),
            msprime.MigrationRateChange(time = T_ACL_GRW, rate = m_EU_AS, matrix_index = (4,3)),
            msprime.MigrationRateChange(time = T_ACL_GRW, rate = m_AS_EU, matrix_index = (3,4)),

            msprime.PopulationParametersChange(time = T_PULSE2-1, initial_size = N_AS0, growth_rate = 0, population_id = 4), # set AS popsize to AS0
            msprime.MassMigration(time = T_PULSE2, source = 4, destination = 0, proportion = m_PULSE2), # Neand1 to EAS pulse of introgression

            msprime.PopulationParametersChange(time = T_EU_AS-1, initial_size = N_EU0, growth_rate = 0, population_id = 3),        # set EU popsize at EU0
            msprime.PopulationParametersChange(time = T_EU_AS-1, initial_size = N_AS0, growth_rate = 0, population_id = 4),         # set AS popsize to AS0
            msprime.MassMigration(time = T_EU_AS, source = 4, destination = 3, proportion = 1.0), # AS merges into EU, now termed "B"
            msprime.MigrationRateChange(time = T_EU_AS, rate = 0), # set all migration rates to zero
            msprime.MigrationRateChange(time = T_EU_AS, rate = m_AF_B, matrix_index = (3, 2)), # Change rate m_AF_EU to m_AF_B, m(row, column) ; m(source, dest) backward in time
            msprime.MigrationRateChange(time = T_EU_AS, rate = m_B_AF, matrix_index = (2, 3)),
            msprime.PopulationParametersChange(time = T_EU_AS, initial_size = N_B, growth_rate = 0, population_id = 3), # set parameters of population "B"
            msprime.MassMigration(time = T_PULSE1, source = 3, destination = 0, proportion = m_PULSE1), # Neand1 to EUR_EAS pulse of introgression
            msprime.MassMigration(time = T_B, source = 3, destination = 2, proportion = 1.0), # Population B merges into Africa at T_B
            msprime.MigrationRateChange(time = T_B, rate = 0), # set all migration rates to zero
            msprime.MassMigration(time = T_N1_N2, source = 1, destination = 0, proportion = 1.0), # N_2 merges with N_1 at T_N1_N2
            msprime.PopulationParametersChange(time = T_N1_N2, initial_size = N_N1, population_id = 0), # set parameters of ancestral Neandertal population
            msprime.PopulationParametersChange(time = T_AF, initial_size = N_A, population_id = 2), # set parameters of ancestral modern human population
            msprime.MassMigration(time = T_DE_N, source = 6, destination = 0, proportion = 1.0), # DE merges with N1
            msprime.PopulationParametersChange(time = T_DE_N, initial_size = N_N1, population_id = 0),
            msprime.MassMigration(time = T_MH_N, source = 0, destination = 2, proportion = 1.0), # Neandertals merge into modern human lineage at time T_MH_N
            msprime.PopulationParametersChange(time = T_MH_N, initial_size = N_A, population_id = 2), # set parameters of ancetral hominin population
            msprime.MassMigration(time = T_MH_CH, source = 5, destination = 2, proportion = 1.0), # Chimp lineage merges into ancestral hominin population at time T_MH_CH
            msprime.PopulationParametersChange(time = T_MH_CH, initial_size = N_A, population_id = 2) # set parameters of ancestral hominin population
    ]

    pulsed_migration = [
            msprime.MigrationRateChange(time = 0, rate = 0),    # turn off all migration from t=0 until T_ACL_GRW (5kya)

            msprime.PopulationParametersChange(time = T_ACL_GRW, initial_size = N_AF0, growth_rate = 0, population_id = 2),        # stop rapid population growth in AF
            msprime.PopulationParametersChange(time = T_ACL_GRW, growth_rate = r_EU1, population_id = 3),        # stop rapid population growth in EU
            msprime.PopulationParametersChange(time = T_ACL_GRW, growth_rate = r_AS1, population_id = 4),        # stop rapid population growth in AS

            msprime.PopulationParametersChange(time = T_PULSE2-1, initial_size = N_AS0, growth_rate = 0, population_id = 4), # set AS popsize to AS0
            msprime.MassMigration(time = T_PULSE2, source = 4, destination = 0, proportion = m_PULSE2), # Neand1 to EAS pulse of introgression

            msprime.MigrationRateChange(time = T_PULSE2, rate = ( (T_EU_AS - T_ACL_GRW) * m_AF_EU), matrix_index = (3,2)),    #start modern migration parameters at 5kya
            msprime.MigrationRateChange(time = T_PULSE2, rate = ( (T_EU_AS - T_ACL_GRW) * m_EU_AF), matrix_index = (2,3)),
            msprime.MigrationRateChange(time = T_PULSE2, rate = ( (T_EU_AS - T_ACL_GRW) * m_AF_AS), matrix_index = (4,2)),
            msprime.MigrationRateChange(time = T_PULSE2, rate = ( (T_EU_AS - T_ACL_GRW) * m_AS_AF), matrix_index = (2,4)),
            msprime.MigrationRateChange(time = T_PULSE2, rate = ( (T_EU_AS - T_ACL_GRW) * m_EU_AS), matrix_index = (4,3)),
            msprime.MigrationRateChange(time = T_PULSE2, rate = ( (T_EU_AS - T_ACL_GRW) * m_AS_EU), matrix_index = (3,4)),
            msprime.MigrationRateChange(time = T_PULSE2+10, rate = 0),                                                         #stop migration parameters

            msprime.PopulationParametersChange(time = T_EU_AS-1, initial_size = N_EU0, growth_rate = 0, population_id = 3),        # set EU popsize at EU0
            msprime.PopulationParametersChange(time = T_EU_AS-1, initial_size = N_AS0, growth_rate = 0, population_id = 4),         # set AS popsize to AS0
            msprime.MassMigration(time = T_EU_AS, source = 4, destination = 3, proportion = 1.0), # AS merges into EU, now termed "B"
            msprime.MigrationRateChange(time = T_EU_AS, rate = 0), # set all migration rates to zero
            msprime.PopulationParametersChange(time = T_EU_AS, initial_size = N_B, growth_rate = 0, population_id = 3), # set parameters of population "B"

            msprime.MigrationRateChange(time = T_B_BN, rate = ( (T_B - T_EU_AS) * m_AF_B ), matrix_index = (3, 2)),         #start ancient migration parameters at 45kya
            msprime.MigrationRateChange(time = T_B_BN, rate = ( (T_B - T_EU_AS) * m_B_AF ), matrix_index = (2, 3)),
            msprime.MigrationRateChange(time = T_B_BN+10, rate = 0),                                                        #stop migration parameters

            msprime.MassMigration(time = T_PULSE1, source = 3, destination = 0, proportion = m_PULSE1), # Neand1 to EUR_EAS pulse of introgression
            msprime.MassMigration(time = T_B, source = 3, destination = 2, proportion = 1.0), # Population B merges into Africa at T_B
            msprime.MigrationRateChange(time = T_B, rate = 0), # set all migration rates to zero
            msprime.MassMigration(time = T_N1_N2, source = 1, destination = 0, proportion = 1.0), # N_2 merges with N_1 at T_N1_N2
            msprime.PopulationParametersChange(time = T_N1_N2, initial_size = N_N1, population_id = 0), # set parameters of ancestral Neandertal population
            msprime.PopulationParametersChange(time = T_AF, initial_size = N_A, population_id = 2), # set parameters of ancestral modern human population
            msprime.MassMigration(time = T_DE_N, source = 6, destination = 0, proportion = 1.0), # DE merges with N1
            msprime.PopulationParametersChange(time = T_DE_N, initial_size = N_N1, population_id = 0),
            msprime.MassMigration(time = T_MH_N, source = 0, destination = 2, proportion = 1.0), # Neandertals merge into modern human lineage at time T_MH_N
            msprime.PopulationParametersChange(time = T_MH_N, initial_size = N_A, population_id = 2), # set parameters of ancetral hominin population
            msprime.MassMigration(time = T_MH_CH, source = 5, destination = 2, proportion = 1.0), # Chimp lineage merges into ancestral hominin population at time T_MH_CH
            msprime.PopulationParametersChange(time = T_MH_CH, initial_size = N_A, population_id = 2) # set parameters of ancestral hominin population
    ]

    if (pulses == 1):
        demographic_events = one_pulse
    elif (pulses == 2):
        demographic_events = two_pulse
    elif (pulses == 3):
        demographic_events = no_modern_migration
    elif (pulses == 4):
        demographic_events = pulsed_migration

    if (haplo == "debug"):
######### DEBUG OPTION #############
            db = msprime.DemographyDebugger(
                        Ne = N_A,
                        population_configurations =  [
                        msprime.PopulationConfiguration(
                                initial_size = N_N1, sample_size = 2), # population 0
                        msprime.PopulationConfiguration(
                                initial_size = N_N2, sample_size = 2), # population 1
                        msprime.PopulationConfiguration(
                                initial_size = N_AF, growth_rate = r_AF2, sample_size = 2), # population 2
                        msprime.PopulationConfiguration(
                                initial_size = N_EU, growth_rate = r_EU2, sample_size = 2), # population 3
                        msprime.PopulationConfiguration(
                                initial_size = N_AS, growth_rate = r_AS2, sample_size = 2), # population 4
                        msprime.PopulationConfiguration(
                                initial_size = N_CH, sample_size = 2), #popultaion 5
                        msprime.PopulationConfiguration(
                                initial_size = N_DE, sample_size = 2)  #population 6
                            ],
                        migration_matrix = migration_matrix,
                        demographic_events = demographic_events
                        )


            print(migration_matrix)
            db.print_history()


    elif (haplo == "F4Dstat"):
####### HAPLOTYPE SIMULATION FOR F4 and DSTAT CALCULATIONS ############
        return msprime.simulate(
                            Ne = N_A,
                            length = length,
                            recombination_rate = recombination_rate,
                            mutation_rate = mutation_rate,
                            samples = samples,
                            population_configurations = population_configurations,
                            migration_matrix = migration_matrix,
                            demographic_events = demographic_events,
                            num_replicates = 20,
                            random_seed = seed
                            )


    else:
####### HAPLOTYPE SIMULATION FOR DESERT DISTRIBUTION ############
            return msprime.simulate(
                            Ne = N_A,
                            length = length,        # Simulate 10Mb chromosomes for looking at desert distribution
                            recombination_rate = recombination_rate,
                            mutation_rate = mutation_rate,
                            samples = samples,
                            population_configurations = population_configurations,
                            migration_matrix = migration_matrix,
                            demographic_events = demographic_events,
                            num_replicates = 1,
                            random_seed = seed
                            )

################################
################################

def Sriram_demography(S_N1, S_N2, S_AF, S_EU, S_AS,
                pulses, seed,
                n1_admix_prop, n2_admix_prop,
                outdir, t_n1_n2, t_n1_sample, t_n2_sample,
                haplo, length,
                m_AF_B, m_B_AF,
                m_AF_AS, m_AS_AF,
                m_AF_EU, m_EU_AF,
                m_EU_AS, m_AS_EU):
    N_A = 10000 # effective population size used for scaling
    N_N1 = 10000 # Altai/Vindija lineage effective population size
    N_N2 = 10000 # Eastern Neandertal effective population size
    N_N_BN = 100 # Neandertal population bottleneck size
    N_B = 10000 # Population size after EUR + ASN join
    N_B_BN = 100 # Population size of EUR_ASN bottleneck (following initial admixture event)
    N_AF = 10000
    #N_EU0 = 1000
    #N_AS0 = 510
    N_CH = 10000
    N_DE = 10000
    # Times are provided in years, so we convert into generations.
    generation_time = 25
    T_MH_CH = 7e6 / generation_time # Chimp_Human lineage split (7 MYA)
    T_MH_N = 700e3 / generation_time # Neanderthal / modern human split (700 KYA)
    T_DE_N = 500e3 / generation_time # Denisovan / Neandertal split time
    T_N1_N2 = ((t_n1_n2) * 1e3) / generation_time # Altai/Vindija and Eastern Neandertal popualation split, default set to 145kya
    T_N1_SAMPLE = ((t_n1_sample * 1e3) / generation_time)
    T_N2_SAMPLE = ((t_n2_sample * 1e3) / generation_time)
    #T_N1_N2 = 175e3 / generation_time # Altai/Vindija and Eastern Neandertal popualation split        ### DIFFERENT FROM TENNESSEN
    #T_AF = 220e3 / generation_time # African population expansion
    T_N_BN = 150e3 / generation_time # Neandertal bottleneck                        ### DIFFERENT FROM TENNESSEN
    T_B = 75e3 / generation_time # out of Africa                                ### DIFFERENT FROM TENNESSEN
    T_PULSE1 = 50e3 / generation_time # Neandertal introgression pulse 1        ### DIFFERENT FROM TENNESSEN
    T_B_BN = 45e3 / generation_time # A 3rd Bottleneck following shortly after initial admixture event        ### SET AT 1800 GEN
    T_EU_AS = 23e3 / generation_time # European / East Asian split
    T_PULSE2 = 18e3 / generation_time # Neandertal introgression pulse 2
    ## We need to work out the starting (diploid) population sizes based on
    ## the growth rates provided for these two populations
    #r_EU = 0.004
    #r_AS = 0.0055
    #N_EU = N_EU0 / math.exp(-r_EU * T_EU_AS)
    #N_AS = N_AS0 / math.exp(-r_AS * T_EU_AS)
    N_EU = 10000
    N_AS = 10000

    ## Migration rates during the various epochs.
    #m_AF_B = 25e-5
    #m_AF_EU = 3e-5
    #m_AF_AS = 1.9e-5
    #m_EU_AS = 9.6e-5
    ####
    m_AF_B = m_AF_B
    m_AF_EU = m_AF_EU
    m_AF_AS = m_AF_AS
    m_EU_AS = m_EU_AS

    m_B_AF = m_B_AF
    m_EU_AF = m_EU_AF
    m_AS_AF = m_AS_AF
    m_AS_EU = m_AS_EU

    ####
    m_N1_N2 = 0
    m_N1_AF = 0
    m_N1_EU = 0
    m_N1_AS = 0
    m_N2_AF = 0
    m_N2_EU = 0
    m_N2_AS = 0
    m_N1_CH = 0
    m_N2_CH = 0
    m_AF_CH = 0
    m_EU_CH = 0
    m_AS_CH = 0
    m_DE_N1 = 0
    m_DE_N2 = 0
    m_DE_AF = 0
    m_DE_EU = 0
    m_DE_AS = 0
    m_DE_CH = 0
    # Set Percent introgression level
    m_PULSE1 = n1_admix_prop # an instantaneous movement of population from EU_AS to Neand1 (backward in time)
    m_PULSE2 = n2_admix_prop
    #print(m_PULSE1)

    # Population IDs correspond to their indexes in the population
    # configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB
    # initially.
    population_configurations = [
        msprime.PopulationConfiguration(
        initial_size = N_N1), # population 0
        msprime.PopulationConfiguration(
        initial_size = N_N2), # population 1
        msprime.PopulationConfiguration(
        initial_size = N_AF), # population 2
        msprime.PopulationConfiguration(
        initial_size = N_EU), # population 3
        msprime.PopulationConfiguration(
        initial_size = N_AS), # population 4
        msprime.PopulationConfiguration(
        initial_size = N_CH), #popultaion 5
        msprime.PopulationConfiguration(
        initial_size = N_DE)  #population 6
    ]
    # specify desired sample
    N1_sample = [msprime.Sample(population = 0, time = T_N1_SAMPLE)] * S_N1
    N2_sample = [msprime.Sample(population = 1, time = T_N2_SAMPLE)] * S_N2
    AF_sample = [msprime.Sample(population = 2, time = 0)] * S_AF
    EU_sample = [msprime.Sample(population = 3, time = 0)] * S_EU
    AS_sample = [msprime.Sample(population = 4, time = 0)] * S_AS
    CH_sample = [msprime.Sample(population = 5, time = 0)] * 2
    DE_sample = [msprime.Sample(population = 6, time = 50e3/ generation_time)] * 2
    samples = N1_sample + N2_sample + AF_sample + EU_sample + AS_sample + CH_sample + DE_sample
    # set up the initial migration matrix
    migration_matrix = [
        [      0, m_N1_N2, m_N1_AF, m_N1_EU, m_N1_AS, m_N1_CH, m_DE_N1],
        [m_N1_N2,       0, m_N2_AF, m_N2_EU, m_N2_AS, m_N2_CH, m_DE_N2],
        [m_N1_AF, m_N2_AF,       0, m_EU_AF, m_AS_AF, m_AF_CH, m_DE_AF],
        [m_N1_EU, m_N2_EU, m_AF_EU,       0, m_AS_EU, m_EU_CH, m_DE_EU],
        [m_N1_AS, m_N2_AS, m_AF_AS, m_EU_AS,       0, m_AS_CH, m_DE_AS],
        [m_N1_CH, m_N2_CH, m_AF_CH, m_EU_CH, m_AS_CH,       0, m_DE_CH],
        [m_DE_N1, m_DE_N2, m_DE_AF, m_DE_EU, m_DE_AS, m_DE_CH,       0]
    ]
    one_pulse = [
        # msprime.MassMigration(time = T_EU_AS, source = 4, destination = 3, proportion = 1.0), # AS merges into EU, now termed "B"
        # msprime.MigrationRateChange(time = T_EU_AS, rate = 0), # set all migration rates to zero
        # msprime.PopulationParametersChange(time = T_EU_AS, initial_size = N_B, growth_rate = 0, population_id = 3), # set parameters of population "B"
        # msprime.PopulationParametersChange(time = T_B_BN, initial_size = N_B_BN, growth_rate = 0, population_id = 3),  # population bottleneck of B begins
        # msprime.PopulationParametersChange(time = T_B_BN + 20, initial_size = N_B, growth_rate = 0, population_id = 3), # population bottleneck ends shortly before initial admixture
        # msprime.MassMigration(time = T_PULSE1, source = 3, destination = 0, proportion = m_PULSE1 ), # Neand1 to EUR_EAS pulse of introgression
        # msprime.MassMigration(time = T_B, source = 3, destination = 2, proportion = 1.0), # Population B merges into Africa at T_B
        # msprime.MigrationRateChange(time = T_B, rate = 0), # set all migration rates to zero
        # msprime.PopulationParametersChange(time = T_N_BN, initial_size = N_N_BN, growth_rate = 0, population_id = 0), # Neand1 has short bottleneck
        # msprime.PopulationParametersChange(time = T_N_BN + 20, initial_size = N_N1, growth_rate = 0, population_id = 0), # Neand1 bottleneck ends
        # msprime.MassMigration(time = T_N1_N2, source = 1, destination = 0, proportion = 1.0), # N_2 merges with N_1 at T_N1_N2
        # msprime.PopulationParametersChange(time = T_N1_N2, initial_size = N_N1, population_id = 0), # set parameters of ancestral Neandertal population
        # msprime.MassMigration(time = T_DE_N, source = 6, destination = 0, proportion = 1.0), # DE merges with N1
        # msprime.PopulationParametersChange(time = T_DE_N, initial_size = N_N1, population_id = 0),
        # msprime.MassMigration(time = T_MH_N, source = 0, destination = 2, proportion = 1.0), # Neandertals merge into modern human lineage at time T_MH_N
        # msprime.PopulationParametersChange(time = T_MH_N, initial_size = N_A, population_id = 2), # set parameters of ancetral hominin population
        # msprime.MassMigration(time = T_MH_CH, source = 5, destination = 2, proportion = 1.0), # Chimp lineage merges into ancestral hominin population at time T_MH_Ch
        # msprime.PopulationParametersChange(time = T_MH_CH, initial_size = N_A, population_id = 2) # set parameters of ancestral hominin population
    ]
    two_pulse = [
        msprime.MassMigration(time = T_PULSE2, source = 4, destination = 0, proportion = m_PULSE2), # Neand1 to EAS pulse of introgression
        msprime.MassMigration(time = T_EU_AS, source = 4, destination = 3, proportion = 1.0), # AS merges into EU, now termed "B"
        msprime.MigrationRateChange(time = T_EU_AS, rate = 0), # set all migration rates to zero
        msprime.MigrationRateChange(time = T_EU_AS, rate = m_AF_B, matrix_index = (3, 2)), # migration between "B" and Africa begins
        msprime.MigrationRateChange(time = T_EU_AS, rate = m_B_AF, matrix_index = (2, 3)),
        msprime.PopulationParametersChange(time = T_EU_AS, initial_size = N_B, growth_rate = 0, population_id = 3), # set parameters of population "B"
        msprime.PopulationParametersChange(time = T_B_BN, initial_size = N_B_BN, growth_rate = 0, population_id = 3), # population bottleneck of B begins
        msprime.PopulationParametersChange(time = T_B_BN + 20, initial_size = N_B, growth_rate = 0, population_id = 3), # population bottleneck of B ends shortly before initial admixture
        msprime.MassMigration(time = T_PULSE1, source = 3, destination = 0, proportion = m_PULSE1), # Neand1 to EUR_EAS pulse of introgression
        msprime.MassMigration(time = T_B, source = 3, destination = 2, proportion = 1.0), # Population B merges into Africa at T_B
        msprime.MigrationRateChange(time = T_B, rate = 0), # set all migration rates to zero
        msprime.MassMigration(time = T_N1_N2, source = 1, destination = 0, proportion = 1.0), # N_2 merges with N_1 at T_N1_N2
        msprime.PopulationParametersChange(time = T_N1_N2, initial_size = N_A, population_id = 0), # set parameters of ancestral Neandertal population
        msprime.PopulationParametersChange(time = T_N_BN, initial_size = N_N_BN, growth_rate = 0, population_id = 0), # Neand1 has short bottleneck
        msprime.PopulationParametersChange(time = T_N_BN + 20, initial_size = N_A, growth_rate = 0, population_id = 0), # Neand1 bottleneck ends
        msprime.MassMigration(time = T_DE_N, source = 6, destination = 0, proportion = 1.0), # DE merges with N1
        msprime.PopulationParametersChange(time = T_DE_N, initial_size = N_A, population_id = 0),
        msprime.MassMigration(time = T_MH_N, source = 0, destination = 2, proportion = 1.0), # Neandertals merge into modern human lineage at time T_MH_N
        msprime.PopulationParametersChange(time = T_MH_N, initial_size = N_A, population_id = 2), # set parameters of ancetral hominin population
        msprime.MassMigration(time = T_MH_CH, source = 5, destination = 2, proportion = 1.0), # Chimp lineage merges into ancestral hominin population at time T_MH_CH
        msprime.PopulationParametersChange(time = T_MH_CH, initial_size = N_A, population_id = 2) # set parameters of ancestral hominin population
    ]


    if (pulses == 1):
        demographic_events = one_pulse
    elif (pulses == 2):
        demographic_events = two_pulse

    if (haplo == "debug"):
######### DEBUG OPTION #############
            db = msprime.DemographyDebugger(
                        Ne = N_A,
                        population_configurations =  [
                        msprime.PopulationConfiguration(
                                initial_size = N_N1, sample_size = 2), # population 0
                        msprime.PopulationConfiguration(
                                initial_size = N_N2, sample_size = 2), # population 1
                        msprime.PopulationConfiguration(
                                initial_size = N_AF, sample_size = 2), # population 2
                        msprime.PopulationConfiguration(
                                initial_size = N_EU, sample_size = 2), # population 3
                        msprime.PopulationConfiguration(
                                initial_size = N_AS, sample_size = 2), # population 4
                        msprime.PopulationConfiguration(
                                initial_size = N_CH, sample_size = 2), #popultaion 5
                        msprime.PopulationConfiguration(
                                initial_size = N_DE, sample_size = 2)  #population 6
                            ],
                        migration_matrix = migration_matrix,
                        demographic_events = demographic_events
                        )
            print(migration_matrix)
            db.print_history()



    elif (haplo == "F4Dstat"):
####### HAPLOTYPE SIMULATION FOR F4 and DSTAT CALCULATIONS ############
            return msprime.simulate(
                            Ne = N_A,
                            length = length,
                            recombination_rate = recombination_rate,
                            mutation_rate = mutation_rate,
                            samples = samples,
                            population_configurations = population_configurations,
                            migration_matrix = migration_matrix,
                            demographic_events = demographic_events,
                            num_replicates = 20,
                            random_seed = seed
                            )

    else:
####### HAPLOTYPE SIMULATION FOR DESERT DISTRIBUTION ############
            return msprime.simulate(
                            Ne = N_A,
                            length = length,        # Simulate 10Mb chromosomes for looking at desert distribution
                            recombination_rate = recombination_rate,
                            mutation_rate = mutation_rate,
                            samples = samples,
                            population_configurations = population_configurations,
                            migration_matrix = migration_matrix,
                            demographic_events = demographic_events,
                            num_replicates = 1,
                            random_seed = seed
                            )

################################
################################

def SplitPop_demography(S_N1, S_N2, S_AF, S_EU, S_AS,
                pulses, seed,
                n1_admix_prop, n2_admix_prop,
                outdir, t_n1_n2, t_n1_sample, t_n2_sample,
                haplo, length,
                m_AF_B, m_B_AF,
                m_AF_AS, m_AS_AF,
                m_AF_EU, m_EU_AF,
                m_EU_AS, m_AS_EU):
    N_A = 7310 # effective population size used for scaling
    N_N1 = 1000 # Altai/Vindija lineage effective population size
    N_N2 = 1000 # Eastern Neandertal effective population size
    N_B = 1861 # Population size after EUR + ASN join
    #N_B_BN = 100 # Population size of EUR_ASN bottleneck (following initial admixture event)
    N_AF0 = 14474
    N_EU0 = 1032
    N_AS0 = 554
    N_CH = 1000
    N_DE = 1000
    N_SP = 100
    # Times are provided in years, so we convert into generations.
    generation_time = 25
    T_MH_CH = 7e6 / generation_time # Chimp_Human lineage split
    T_MH_N = 700e3 / generation_time # Neanderthal / modern human split
    T_DE_N = 500e3 / generation_time # Denisovan / Neandertal split time
    T_N1_N2 = ((t_n1_n2) * 1e3) / generation_time # Altai/Vindija and Eastern Neandertal popualation split, default set to 145kya
    T_N1_SAMPLE = ((t_n1_sample * 1e3) / generation_time)
    T_N2_SAMPLE = ((t_n2_sample * 1e3) / generation_time)
    T_AF = 148e3 / generation_time # African population expansion
    T_B = 100e3 / generation_time # out of Africa
    T_PULSE1 = 55e3 / generation_time # Neandertal introgression pulse 1
    T_SP = 54e3 /generation_time # split population off of B
    #T_B_BN = 45e3 / generation_time # A 3rd Bottleneck following shortly after initial admixture event
    T_EU_AS = 23e3 / generation_time # European / East Asian split
    T_PULSE2 = 18e3 / generation_time # Neandertal introgression pulse 2
    T_ACL_GRW = 5.1e3 / generation_time # time of rapid population growth
    # We need to work out the starting (diploid) population sizes based on
    # the growth rates provided for these two populations
    r_EU1 = 0.00307
    r_AS1 = 0.0048
    #N_EU1 = N_EU0 / math.exp(-r_EU1 * (T_EU_AS-T_ACL_GRW))
    #N_AS1 = N_AS0 / math.exp(-r_AS1 * (T_EU_AS-T_ACL_GRW))
    r_EU2 = 0.0195        # rapid population growth at 5.1kya
    r_AS2 = 0.025        # rapid population growth at 5.1kya
    r_AF2 = 0.0166        # rapid popualtion growth at 5.1kya
    #N_EU2 = N_EU1 / math.exp(-r_EU2 * T_ACL_GRW)        # should reach ~512k
    #N_AS2 = N_AS1 / math.exp(-r_AS2 * T_ACL_GRW)        # should reach ~640k
    #N_AF2 = N_AF / math.exp(-r_AF2 * T_ACL_GRW)        # should reach ~424k
    N_EU = 512e3
    N_AS = 640e3
    N_AF = 424e3

    ## Migration rates during the various epochs.
    #m_AF_B = 15e-5
    #m_AF_EU = 2.5e-5
    #m_AF_AS = 0.78e-5
    #m_EU_AS = 3.11e-5
    ####
    m_AF_B = m_AF_B
    m_AF_EU = m_AF_EU
    m_AF_AS = m_AF_AS
    m_EU_AS = m_EU_AS

    m_B_AF = m_B_AF
    m_EU_AF = m_EU_AF
    m_AS_AF = m_AS_AF
    m_AS_EU = m_AS_EU

    ####
    m_N1_N2 = 0
    m_N1_AF = 0
    m_N1_EU = 0
    m_N1_AS = 0
    m_N2_AF = 0
    m_N2_EU = 0
    m_N2_AS = 0
    m_N1_CH = 0
    m_N2_CH = 0
    m_AF_CH = 0
    m_EU_CH = 0
    m_AS_CH = 0
    m_DE_N1 = 0
    m_DE_N2 = 0
    m_DE_AF = 0
    m_DE_EU = 0
    m_DE_AS = 0
    m_DE_CH = 0
    m_SP_N1 = 0
    m_SP_N2 = 0
    m_SP_AF = 0
    m_SP_EU = 0
    m_SP_AS = 0
    m_SP_CH = 0
    m_SP_DE = 0

    # Set Percent introgression level
    m_PULSE1 = n1_admix_prop # an instantaneous movement of population from EU_AS to Neand1 (backward in time)
    m_PULSE2 = n2_admix_prop
    #print(m_PULSE1)

    # Population IDs correspond to their indexes in the population
    # configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB
    # initially.
    population_configurations = [
        msprime.PopulationConfiguration(
        initial_size = N_N1), # population 0
        msprime.PopulationConfiguration(
        initial_size = N_N2), # population 1
        msprime.PopulationConfiguration(
        initial_size = N_AF, growth_rate = r_AF2), # population 2
        msprime.PopulationConfiguration(
        initial_size = N_EU, growth_rate = r_EU2), # population 3
        msprime.PopulationConfiguration(
        initial_size = N_AS, growth_rate = r_AS2), # population 4
        msprime.PopulationConfiguration(
        initial_size = N_CH), # popultaion 5
        msprime.PopulationConfiguration(
        initial_size = N_DE),  # population 6
        msprime.PopulationConfiguration(
        initial_size = N_SP)  # population 7
    ]

    # specify desired sample
    N1_sample = [msprime.Sample(population = 0, time = T_N1_SAMPLE)] * S_N1
    N2_sample = [msprime.Sample(population = 1, time = T_N2_SAMPLE)] * S_N2
    AF_sample = [msprime.Sample(population = 2, time = 0)] * S_AF
    EU_sample = [msprime.Sample(population = 3, time = 0)] * S_EU
    AS_sample = [msprime.Sample(population = 4, time = 0)] * S_AS
    CH_sample = [msprime.Sample(population = 5, time = 0)] * 2
    DE_sample = [msprime.Sample(population = 6, time = 50e3/ generation_time)] * 2
    samples = N1_sample + N2_sample + AF_sample + EU_sample + AS_sample + CH_sample + DE_sample
    # set up the initial migration matrix
    migration_matrix = [
        [      0, m_N1_N2, m_N1_AF, m_N1_EU, m_N1_AS, m_N1_CH, m_DE_N1, m_SP_N1],
        [m_N1_N2,       0, m_N2_AF, m_N2_EU, m_N2_AS, m_N2_CH, m_DE_N2, m_SP_N2],
        [m_N1_AF, m_N2_AF,       0, m_EU_AF, m_AS_AF, m_AF_CH, m_DE_AF, m_SP_AF],
        [m_N1_EU, m_N2_EU, m_AF_EU,       0, m_AS_EU, m_EU_CH, m_DE_EU, m_SP_EU],
        [m_N1_AS, m_N2_AS, m_AF_AS, m_EU_AS,       0, m_AS_CH, m_DE_AS, m_SP_AS],
        [m_N1_CH, m_N2_CH, m_AF_CH, m_EU_CH, m_AS_CH,       0, m_DE_CH, m_SP_CH],
        [m_DE_N1, m_DE_N2, m_DE_AF, m_DE_EU, m_DE_AS, m_DE_CH,       0, m_SP_DE],
        [m_SP_N1, m_SP_N2, m_SP_AF, m_SP_EU, m_SP_AS, m_SP_CH, m_SP_DE,       0]
    ]
    one_pulse = [
        # msprime.PopulationParametersChange(time = T_ACL_GRW, initial_size = N_AF0, growth_rate = 0, population_id = 2),        # stop rapid population growth in AF
        # msprime.PopulationParametersChange(time = T_ACL_GRW, growth_rate = r_EU1, population_id = 3),        # stop rapid population growth in EU
        # msprime.PopulationParametersChange(time = T_ACL_GRW, growth_rate = r_AS1, population_id = 4),        # stop rapid population growth in AS
        # msprime.PopulationParametersChange(time = T_PULSE2-1, initial_size = N_AS0, growth_rate = 0, population_id = 4), # set AS popsize to AS0
        # msprime.PopulationParametersChange(time = T_EU_AS-1, initial_size = N_EU0, growth_rate = 0, population_id = 3),        # set EU popsize at EU0
        # msprime.PopulationParametersChange(time = T_EU_AS-1, initial_size = N_AS0, growth_rate = 0, population_id = 4),         # set AS popsize to AS0
        # msprime.MassMigration(time = T_EU_AS, source = 4, destination = 3, proportion = 1.0), # AS merges into EU, now termed "B"
        # msprime.MigrationRateChange(time = T_EU_AS, rate = 0), # set all migration rates to zero
        # msprime.MigrationRateChange(time = T_EU_AS, rate = m_AF_B, matrix_index = (2, 3)), # migration between "B" and Africa begins
        # msprime.MigrationRateChange(time = T_EU_AS, rate = m_AF_B, matrix_index = (3, 2)),
        # msprime.PopulationParametersChange(time = T_EU_AS, initial_size = N_B, growth_rate = 0, population_id = 3), # set parameters of population "B"
        # msprime.MassMigration(time = T_SP, source = 3, destination = 7, proportion = 0.1),  # Population splits off of EUR_EAS, 10% migration rate to SplitPop
        # msprime.PopulationParametersChange(time = T_SP, initial_size = N_SP, growth_rate = 0, population_id = 7), # SplitPop is set at size Ne=100
        # msprime.MassMigration(time = T_PULSE1, source = 7, destination = 0, proportion = m_PULSE1 ), # Neand1 to SplitPop pulse of introgression
        # msprime.MassMigration(time = T_B, source = 3, destination = 2, proportion = 1.0), # Population B merges into African at T_B
        # msprime.MassMigration(time = T_B, source = 7, destination = 2, proportion = 1.0), # SplitPop merges into African at T_B
        # msprime.MigrationRateChange(time = T_B, rate = 0), # set all migration rates to zero
        # msprime.PopulationParametersChange(time = T_AF, initial_size = N_A, population_id = 2), # set parameters of ancestral modern human population
        # msprime.MassMigration(time = T_N1_N2, source = 1, destination = 0, proportion = 1.0), # N_2 merges with N_1 at T_N1_N2
        # msprime.PopulationParametersChange(time = T_N1_N2, initial_size = N_N1, population_id = 0), # set parameters of ancestral Neandertal population
        # msprime.MassMigration(time = T_DE_N, source = 6, destination = 0, proportion = 1.0), # DE merges with N1
        # msprime.PopulationParametersChange(time = T_DE_N, initial_size = N_N1, population_id = 0),
        # msprime.MassMigration(time = T_MH_N, source = 0, destination = 2, proportion = 1.0), # Neandertals merge into modern human lineage at time T_MH_N
        # msprime.PopulationParametersChange(time = T_MH_N, initial_size = N_A, population_id = 2), # set parameters of ancetral hominin population
        # msprime.MassMigration(time = T_MH_CH, source = 5, destination = 2, proportion = 1.0), # Chimp lineage merges into ancestral hominin population at time T_MH_Ch
        # msprime.PopulationParametersChange(time = T_MH_CH, initial_size = N_A, population_id = 2) # set parameters of ancestral hominin population
    ]

    two_pulse = [
        msprime.PopulationParametersChange(time = T_ACL_GRW, initial_size = N_AF0, growth_rate = 0, population_id = 2),        # stop rapid population growth in AF
        msprime.PopulationParametersChange(time = T_ACL_GRW, growth_rate = r_EU1, population_id = 3),        # stop rapid population growth in EU
        msprime.PopulationParametersChange(time = T_ACL_GRW, growth_rate = r_AS1, population_id = 4),        # stop rapid population growth in AS
        msprime.PopulationParametersChange(time = T_PULSE2-1, initial_size = N_AS0, growth_rate = 0, population_id = 4), # set AS popsize to AS0
        msprime.MassMigration(time = T_PULSE2, source = 4, destination = 0, proportion = m_PULSE2), # Neand1 to EAS pulse of introgression
        msprime.PopulationParametersChange(time = T_EU_AS-1, initial_size = N_EU0, growth_rate = 0, population_id = 3),        # set EU popsize at EU0
        msprime.PopulationParametersChange(time = T_EU_AS-1, initial_size = N_AS0, growth_rate = 0, population_id = 4),         # set AS popsize to AS0
        msprime.MassMigration(time = T_EU_AS, source = 4, destination = 3, proportion = 1.0), # AS merges into EU, now termed "B"
        msprime.MigrationRateChange(time = T_EU_AS, rate = 0), # set all migration rates to zero
        msprime.MigrationRateChange(time = T_EU_AS, rate = m_AF_B, matrix_index = (3, 2)), # migration between "B" and Africa begins
        msprime.MigrationRateChange(time = T_EU_AS, rate = m_B_AF, matrix_index = (2, 3)),
        msprime.PopulationParametersChange(time = T_EU_AS, initial_size = N_B, growth_rate = 0, population_id = 3), # set parameters of population "B"
        msprime.MassMigration(time = T_SP, source = 3, destination = 7, proportion = 0.1),  # Population splits off of EUR_EAS, 10% migration rate to SplitPop
        msprime.PopulationParametersChange(time = T_SP, initial_size = N_SP, growth_rate = 0, population_id = 7), # SplitPop is set at size Ne=100
        msprime.MassMigration(time = T_PULSE1, source = 7, destination = 0, proportion = m_PULSE1 ), # Neand1 to SplitPop pulse of introgression
        msprime.MassMigration(time = T_B, source = 3, destination = 2, proportion = 1.0), # Population B merges into African at T_B
        msprime.MassMigration(time = T_B, source = 7, destination = 2, proportion = 1.0), # SplitPop merges into African at T_B
        msprime.MigrationRateChange(time = T_B, rate = 0), # set all migration rates to zero
        msprime.MassMigration(time = T_N1_N2, source = 1, destination = 0, proportion = 1.0), # N_2 merges with N_1 at T_N1_N2
        msprime.PopulationParametersChange(time = T_N1_N2, initial_size = N_N1, population_id = 0), # set parameters of ancestral Neandertal population
        msprime.PopulationParametersChange(time = T_AF, initial_size = N_A, population_id = 2), # set parameters of ancestral modern human population
        msprime.MassMigration(time = T_DE_N, source = 6, destination = 0, proportion = 1.0), # DE merges with N1
        msprime.PopulationParametersChange(time = T_DE_N, initial_size = N_N1, population_id = 0),
        msprime.MassMigration(time = T_MH_N, source = 0, destination = 2, proportion = 1.0), # Neandertals merge into modern human lineage at time T_MH_N
        msprime.PopulationParametersChange(time = T_MH_N, initial_size = N_A, population_id = 2), # set parameters of ancetral hominin population
        msprime.MassMigration(time = T_MH_CH, source = 5, destination = 2, proportion = 1.0), # Chimp lineage merges into ancestral hominin population at time T_MH_CH
        msprime.PopulationParametersChange(time = T_MH_CH, initial_size = N_A, population_id = 2) # set parameters of ancestral hominin population
    ]

    if (pulses == 1):
        demographic_events = one_pulse
    elif (pulses == 2):
        demographic_events = two_pulse

    if (haplo == "debug"):
######### DEBUG OPTION #############
            db = msprime.DemographyDebugger(
                        Ne = N_A,
                        population_configurations =  [
                        msprime.PopulationConfiguration(
                                initial_size = N_N1, sample_size = 2), # population 0
                        msprime.PopulationConfiguration(
                                initial_size = N_N2, sample_size = 2), # population 1
                        msprime.PopulationConfiguration(
                                initial_size = N_AF, growth_rate = r_AF2, sample_size = 2), # population 2
                        msprime.PopulationConfiguration(
                                initial_size = N_EU, growth_rate = r_EU2, sample_size = 2), # population 3
                        msprime.PopulationConfiguration(
                                initial_size = N_AS, growth_rate = r_AS2, sample_size = 2), # population 4
                        msprime.PopulationConfiguration(
                                initial_size = N_CH, sample_size = 2), #popultaion 5
                        msprime.PopulationConfiguration(
                                initial_size = N_DE, sample_size = 2),  #population 6
                        msprime.PopulationConfiguration(
                                initial_size = N_SP, sample_size = 2)        #population 7
                            ],
                        migration_matrix = migration_matrix,
                        demographic_events = demographic_events
                        )
            print(migration_matrix)
            db.print_history()


    elif (haplo == "F4Dstat"):
####### HAPLOTYPE SIMULATION FOR F4 and DSTAT CALCULATIONS ############
            return msprime.simulate(
                            Ne = N_A,
                            length = length,
                            recombination_rate = recombination_rate,
                            mutation_rate = mutation_rate,
                            samples = samples,
                            population_configurations = population_configurations,
                            migration_matrix = migration_matrix,
                            demographic_events = demographic_events,
                            num_replicates = 20,
                            random_seed = seed
                            )

    else:
####### HAPLOTYPE SIMULATION FOR DESERT DISTRIBUTION ############
            return msprime.simulate(
                            Ne = N_A,
                            length = length,        # Simulate 10Mb chromosomes for looking at desert distribution
                            recombination_rate = recombination_rate,
                            mutation_rate = mutation_rate,
                            samples = samples,
                            population_configurations = population_configurations,
                            migration_matrix = migration_matrix,
                            demographic_events = demographic_events,
                            num_replicates = 1,
                            random_seed = seed
                            )



#######################################
#######################################

def out_of_africa(S_AF, S_EU, S_AS,
                pulses, seed,
                n1_admix_prop, n2_admix_prop,
                outdir, t_n1_n2,
                haplo, length,
                m_AF_B, m_B_AF,
                m_AF_AS, m_AS_AF,
                m_AF_EU, m_EU_AF,
                m_EU_AS, m_AS_EU):
    # First we set out the maximum likelihood values of the various parameters
    # given in Table 1.
    N_A = 7300
    N_B = 2100
    N_AF = 12300
    N_oAF = 1000
    N_EU0 = 1000
    N_AS0 = 510
    N_CH = 1000
    # Times are provided in years, so we convert into generations.
    generation_time = 25
    T_AF = 220e3 / generation_time
    T_oAF_AF = 200e3/generation_time    # creat a split between and African outgroup and main African group
    T_B = 140e3 / generation_time
    T_EU_AS = 21.2e3 / generation_time
    # We need to work out the starting (diploid) population sizes based on
    # the growth rates provided for these two populations
    r_EU = 0.004
    r_AS = 0.0055
    N_EU = N_EU0 / math.exp(-r_EU * T_EU_AS)
    N_AS = N_AS0 / math.exp(-r_AS * T_EU_AS)
    # Migration rates during the various epochs.
    # m_AF_B = 25e-5
    # m_AF_EU = 3e-5
    # m_AF_AS = 1.9e-5
    # m_EU_AS = 9.6e-5

    m_AF_B = m_AF_B
    m_AF_EU = m_AF_EU
    m_AF_AS = m_AF_AS
    m_EU_AS = m_EU_AS

    m_B_AF = m_B_AF
    m_EU_AF = m_EU_AF
    m_AS_AF = m_AS_AF
    m_AS_EU = m_AS_EU


    # Population IDs correspond to their indexes in the population
    # configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB
    # initially.
    population_configurations = [
        msprime.PopulationConfiguration(
            initial_size=N_oAF),
        msprime.PopulationConfiguration(
            initial_size=N_AF),
        msprime.PopulationConfiguration(
            initial_size=N_EU, growth_rate=r_EU),
        msprime.PopulationConfiguration(
            initial_size=N_AS, growth_rate=r_AS),
        msprime.PopulationConfiguration(
            initial_size=N_CH)
    ]


    # specify desired sample
    # N1_sample = [msprime.Sample(population = 0, time = 50e3 / generation_time)] * S_N1
    # N2_sample = [msprime.Sample(population = 1, time = 50e3 / generation_time)] * S_N2
    oAF_sample = [msprime.Sample(population = 0, time = 0)] * 2
    AF_sample = [msprime.Sample(population = 1, time = 0)] * S_AF
    EU_sample = [msprime.Sample(population = 2, time = 0)] * S_EU
    AS_sample = [msprime.Sample(population = 3, time = 0)] * S_AS
    CH_sample = [msprime.Sample(population = 4, time = 0)] * 2
    # DE_sample = [msprime.Sample(population = 6, time = 50e3/ generation_time)] * 2
    # samples = N1_sample + N2_sample + AF_sample + EU_sample + AS_sample + CH_sample + DE_sample
    samples = oAF_sample + AF_sample + EU_sample + AS_sample + CH_sample

    migration_matrix = [
        [0,     0,      0,      0,      0],
        [0,     0, m_EU_AF, m_AS_AF,    0],
        [0, m_AF_EU,      0, m_AS_EU,   0],
        [0, m_AF_AS, m_EU_AS,       0,  0],
        [0,     0,      0,      0,      0]
    ]


    demographic_events = [
        # CEU and CHB merge into B with rate changes at T_EU_AS
        msprime.MassMigration(
            time=T_EU_AS, source=3, destination=2, proportion=1.0),
        msprime.MigrationRateChange(
            time=T_EU_AS, rate=0),
        msprime.MigrationRateChange(
            time=T_EU_AS, rate=m_AF_B, matrix_index=(2, 1)),
        msprime.MigrationRateChange(
            time=T_EU_AS, rate=m_AF_B, matrix_index=(1, 2)),
        msprime.PopulationParametersChange(
            time=T_EU_AS, initial_size=N_B, growth_rate=0, population_id=2),
        # Population B merges into YRI at T_B
        msprime.MassMigration(
            time=T_B, source=2, destination=1, proportion=1.0),
        msprime.MassMigration(
            time=T_oAF_AF, source = 0, destination = 1, proportion=1.0),
        # Size changes to N_A at T_AF
        msprime.PopulationParametersChange(
            time=T_AF, initial_size=N_A, population_id=1)
    ]
    # Use the demography debugger to print out the demographic history
    # that we have just described.
    if (haplo == "debug"):
######### DEBUG OPTION #############
            db = msprime.DemographyDebugger(
                        Ne = N_A,
                        population_configurations =  [
                        msprime.PopulationConfiguration(
                            initial_size=N_oAF, sample_size = 2),
                        msprime.PopulationConfiguration(
                            initial_size=N_AF, sample_size = 2),
                        msprime.PopulationConfiguration(
                            initial_size=N_EU, growth_rate=r_EU, sample_size = 2),
                        msprime.PopulationConfiguration(
                            initial_size=N_AS, growth_rate=r_AS, sample_size = 2),
                        msprime.PopulationConfiguration(
                            initial_size=N_CH, sample_size = 2)
                        ],
                        migration_matrix = migration_matrix,
                        demographic_events = demographic_events
                        )
            print(migration_matrix)
            db.print_history()

    elif (haplo == "F4Dstat"):
####### HAPLOTYPE SIMULATION FOR F4 and DSTAT CALCULATIONS ############
            return msprime.simulate(
                            Ne = N_A,
                            length = length,
                            recombination_rate = recombination_rate,
                            mutation_rate = mutation_rate,
                            samples = samples,
                            population_configurations = population_configurations,
                            migration_matrix = migration_matrix,
                            demographic_events = demographic_events,
                            num_replicates = 20,
                            random_seed = seed
                            )



#######################################
#######################################
