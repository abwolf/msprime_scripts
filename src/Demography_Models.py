import msprime
import math
import numpy as np


class Population(object):
    def __init__(self,
                 abbreviation,
                 population_size,
                 growth_rate=0.0,
                 samples=2,
                 generations=0,
                 long_name=''):
        self.abbreviation = abbreviation
        self.size = population_size
        self.rate = growth_rate
        self.samples = samples
        self.generations = generations
        self.long_name = long_name

    def get_indices(self, index):
        return [index] * self.samples

    def get_sample(self, index):
        return [msprime.Sample(
            population=index,
            time=self.generations)] * self.samples

    def get_configuration(self):
        """get msprime configuration for the population"""
        return msprime.PopulationConfiguration(
            initial_size=self.size,
            growth_rate=self.rate)

    def get_debug_configuration(self, includeRate=False):
        if includeRate:
            return msprime.PopulationConfiguration(
                initial_size=self.size,
                sample_size=2,
                growth_rate=self.rate)

        else:
            return msprime.PopulationConfiguration(
                initial_size=self.size,
                sample_size=2)


class Base_demography(object):
    def __init__(self, options):
        self.S_N1 = options.s_n1
        self.S_N2 = options.s_n2
        self.options = options

        # Set Percent introgression level
        # an instantaneous movement of population from
        # EU_AS to Neand1 (backward in time)
        self.m_PULSE1 = self.options.n1_admix_prop
        self.m_PULSE2 = self.options.n2_admix_prop
        self.set_constants()
        self.set_populations()
        self.set_demographic_events()

    def simulate(self, replicates):
        return msprime.simulate(
            Ne=self.N_A,
            length=self.length,
            recombination_rate=self.recombination_rate,
            mutation_rate=self.mutation_rate,
            samples=self.get_samples(),
            population_configurations=self.get_population_configuration(),
            migration_matrix=self.get_migration_matrix(),
            demographic_events=self.get_demographic_events(),
            num_replicates=replicates,
            random_seed=self.options.seed
        )

    def print_debug(self, output):
        print(np.array(self.get_migration_matrix()), file=output)
        msprime.DemographyDebugger(
            Ne=self.N_A,
            population_configurations=self.get_debug_configuration(),
            migration_matrix=self.get_migration_matrix(),
            demographic_events=self.get_demographic_events()
        ).print_history(output)

    def set_constants(self):
        self.generation_time = 25  # years/generation
        self.length = self.options.length
        # 2e-8 original from Rajiv, not sure where it is from.
        self.recombination_rate = 1e-8

        # used in Vernot 2014 and 2016
        # matches 0migration and CWZ num_s_star_snps,
        # matches 0migration s_star scores
        # self.mutation_rate = 2.5e-8
        # value best matches CWZ data n_region_ind_snps and s_star scores
        self.mutation_rate = 1.2e-8

        # effective population size used for scaling
        self.N_A = 7310
        # Altai/Vindija lineage effective population size
        self.N_N1 = 1000
        # Eastern Neandertal effective population size
        self.N_N2 = 1000
        # Population size after EUR + ASN join
        self.N_B = 1861
        # Population size of EUR_ASN bottleneck
        # (following initial admixture event)
        self.N_AF0 = 14474
        self.N_EU0 = 1032
        self.N_AS0 = 554
        self.N_CH = 1000
        self.N_DE = 1000

        # Chimp_Human lineage split
        self.T_MH_CH = 7e6 / self.generation_time
        # Neanderthal / modern human split
        self.T_MH_N = 700e3 / self.generation_time
        # Denisovan / Neandertal split time
        self.T_DE_N = 500e3 / self.generation_time
        # Altai/Vindija and Eastern Neandertal popualation split,
        self.T_N1_N2 = ((self.options.t_n1_n2) * 1e3) / self.generation_time
        self.T_N1_SAMPLE = ((self.options.t_n1_sample) * 1e3) \
            / self.generation_time
        self.T_N2_SAMPLE = ((self.options.t_n2_sample) * 1e3) \
            / self.generation_time
        # African population expansion
        self.T_AF = 148e3 / self.generation_time
        # out of Africa
        self.T_B = 100e3 / self.generation_time
        # Neandertal introgression pulse 1
        self.T_PULSE1 = 55e3 / self.generation_time
        # A 3rd bottleneck following shortly after initial admixture event
        self.T_B_BN = 45e3 / self.generation_time
        # European / East Asian split
        self.T_EU_AS = 23e3 / self.generation_time
        # Neandertal introgression pulse 2
        self.T_PULSE2 = 18e3 / self.generation_time
        # time of rapid population growth
        self.T_ACL_GRW = 5.1e3 / self.generation_time
        # We need to work out the starting (diploid) population sizes based on
        # the growth rates provided for these two populations
        self.r_EU1 = 0.00307
        self.r_AS1 = 0.0048
        # rapid popualtion growth at 5.1kya
        self.N_EU = 512e3
        self.N_AS = 640e3
        self.N_AF = 424e3

    def set_populations(self):
        """set the self.populations variable
        a list of Populations for the simulation"""
        self.populations = [
            # Altai/Vindija lineage effective population size
            Population(abbreviation='N1',
                       population_size=self.N_N1,
                       samples=self.S_N1,
                       generations=self.T_N1_SAMPLE,
                       long_name='Neand1'),

            # Eastern Neandertal effective population size
            Population(abbreviation='N2',
                       population_size=self.N_N2,
                       samples=self.S_N2,
                       generations=self.T_N2_SAMPLE,
                       long_name='Neand2'),

            Population(abbreviation='AF',
                       population_size=self.N_AF,
                       growth_rate=0.0166,
                       samples=self.options.AF_sample_size,
                       long_name='AFR'),

            Population(abbreviation='EU',
                       population_size=self.N_EU,
                       growth_rate=0.0195,
                       samples=self.options.EU_sample_size,
                       long_name='EUR'),

            Population(abbreviation='AS',
                       population_size=self.N_AS,
                       growth_rate=0.025,
                       samples=self.options.AS_sample_size,
                       long_name='ASN'),

            Population(abbreviation='CH',
                       population_size=self.N_CH,
                       long_name='Chimp'),

            Population(abbreviation='DE',
                       population_size=self.N_DE,
                       generations=50e3 / self.generation_time,
                       long_name='Deni'),

        ]

    def set_demographic_events(self):
        """set the list of demographic events in self.events"""
        ids = self.get_population_map()
        migrations = self.get_later_migrations()

        self.events = []

        self.events += [
            # stop rapid population growth in AF
            msprime.PopulationParametersChange(
                time=self.T_ACL_GRW,
                initial_size=self.N_AF0,
                growth_rate=0,
                population_id=ids['AF']),

            # stop rapid population growth in EU
            msprime.PopulationParametersChange(
                time=self.T_ACL_GRW,
                growth_rate=self.r_EU1,
                population_id=ids['EU']),

            # stop rapid population growth in AS
            msprime.PopulationParametersChange(
                time=self.T_ACL_GRW,
                growth_rate=self.r_AS1,
                population_id=ids['AS']),

            # set AS popsize to AS0
            msprime.PopulationParametersChange(
                time=self.T_PULSE2-1,
                initial_size=self.N_AS0,
                growth_rate=0,
                population_id=ids['AS']),

            # Neand1 to EAS pulse of introgression
            msprime.MassMigration(
                time=self.T_PULSE2,
                source=ids['AS'],
                destination=ids['N1'],
                proportion=self.m_PULSE2),

            # set EU popsize at EU0
            msprime.PopulationParametersChange(
                time=self.T_EU_AS-1,
                initial_size=self.N_EU0,
                growth_rate=0,
                population_id=ids['EU']),

            # set AS popsize to AS0
            msprime.PopulationParametersChange(
                time=self.T_EU_AS-1,
                initial_size=self.N_AS0,
                growth_rate=0,
                population_id=ids['AS']),

            # AS merges into EU, now termed "B"
            msprime.MassMigration(
                time=self.T_EU_AS,
                source=ids['AS'],
                destination=ids['EU'],
                proportion=1.0)]

        ids['B'] = ids['EU']

        self.events += [
            # set all migration rates to zero
            msprime.MigrationRateChange(
                time=self.T_EU_AS,
                rate=0),

            # migration between "B" and Africa begins
            # note transposed order of indices
            msprime.MigrationRateChange(
                time=self.T_EU_AS,
                rate=migrations.get('AF_B', 0.0),
                matrix_index=(ids['B'], ids['AF'])),
            msprime.MigrationRateChange(
                time=self.T_EU_AS,
                rate=migrations.get('B_AF', 0.0),
                matrix_index=(ids['AF'], ids['B'])),

            # set parameters of population "B"
            msprime.PopulationParametersChange(
                time=self.T_EU_AS,
                initial_size=self.N_B,
                growth_rate=0,
                population_id=ids['B']),

            # Neand1 to EUR_EAS pulse of introgression
            msprime.MassMigration(
                time=self.T_PULSE1,
                source=ids['B'],
                destination=ids['N1'],
                proportion=self.m_PULSE1),

            # Population B merges into Africa at T_B
            msprime.MassMigration(
                time=self.T_B,
                source=ids['B'],
                destination=ids['AF'],
                proportion=1.0),

            # set all migration rates to zero
            msprime.MigrationRateChange(
                time=self.T_B,
                rate=0),

            # N_2 merges with N_1 at T_N1_N2
            msprime.MassMigration(
                time=self.T_N1_N2,
                source=ids['N2'],
                destination=ids['N1'],
                proportion=1.0),

            # set parameters of ancestral Neandertal population
            msprime.PopulationParametersChange(
                time=self.T_N1_N2,
                initial_size=self.N_N1,
                population_id=ids['N1']),

            # set parameters of ancestral modern human population
            msprime.PopulationParametersChange(
                time=self.T_AF,
                initial_size=self.N_A,
                population_id=ids['AF']),

            # DE merges with N1
            msprime.MassMigration(
                time=self.T_DE_N,
                source=ids['DE'],
                destination=ids['N1'],
                proportion=1.0),
            msprime.PopulationParametersChange(
                time=self.T_DE_N,
                initial_size=self.N_N1,
                population_id=ids['N1']),

            # Neandertals merge into modern human lineage at time T_MH_N
            msprime.MassMigration(
                time=self.T_MH_N,
                source=ids['N1'],
                destination=ids['AF'],
                proportion=1.0),

            # set parameters of ancetral hominin population
            msprime.PopulationParametersChange(
                time=self.T_MH_N,
                initial_size=self.N_A,
                population_id=ids['AF']),

            # Chimp lineage merges into ancestral hominin
            # population at time T_MH_CH
            msprime.MassMigration(
                time=self.T_MH_CH,
                source=ids['CH'],
                destination=ids['AF'],
                proportion=1.0),

            # set parameters of ancestral hominin population
            msprime.PopulationParametersChange(
                time=self.T_MH_CH,
                initial_size=self.N_A,
                population_id=ids['AF'])
        ]

    def get_population_map(self):
        result = {}
        for i, pop in enumerate(self.populations):
            result[pop.abbreviation] = i
        return result

    def get_long_name_map(self):
        return [pop.long_name for pop in self.populations]

    def get_population_configuration(self):
        """build array of population configurations from self.populations"""
        return [pop.get_configuration()
                for pop in self.populations]

    def get_debug_configuration(self):
        """build array of population configurations for debugging"""
        return [pop.get_debug_configuration()
                for pop in self.populations]

    def get_samples(self):
        result = []
        for i, p in enumerate(self.populations):
            result += p.get_sample(i)
        return result

    def get_sample_indices(self):
        result = []
        for i, p in enumerate(self.populations):
            result += p.get_indices(i)
        return result

    def get_initial_migrations(self):
        '''convert the list of strings from option parser into
        a dictionary with the keys all validated as a Population'''
        result = {}
        validNames = [pop.abbreviation for pop in self.populations]
        for entry in self.options.initial_migrations:
            toks = entry.split('_')
            for i in range(2):
                if toks[0] not in validNames:
                    raise ValueError('{} not a known population'
                                     .format(toks[0]))
            result['{}_{}'.format(toks[0], toks[1])] = float(toks[2])

        return result

    def get_later_migrations(self):
        '''convert list of strings from option parser into a
        dictionary.  Keys not validated as populations'''
        result = {}
        for entry in self.options.later_migrations:
            toks = entry.split('_')
            result['{}_{}'.format(toks[0], toks[1])] = float(toks[2])
        return result

    def get_migration_matrix(self):
        result = []
        migrations = self.get_initial_migrations()
        names = [pop.abbreviation for pop in self.populations]
        for r in names:
            row = []
            for c in names:
                row.append(migrations.get('{}_{}'.format(c, r), 0))
            result.append(row)
        return result

    def get_demographic_events(self):
        return sorted(self.events, key=lambda e: e.time)


class Tenn_demography(Base_demography):
    """Small change from base as base was made from Tenn"""

    def get_debug_configuration(self):
        return [pop.get_debug_configuration(includeRate=True)
                for pop in self.populations]


class Tenn_no_modern_migration(Tenn_demography):
    """Modified version of Tennison without modern migration"""
    def set_demographic_events(self):
        """set the list of demographic events in self.events"""
        ids = self.get_population_map()
        # NOTE: this is initial migrations to match legacy code
        # it is more logical to be later migrations.
        # would also want to change migration assignment after B
        migrations = self.get_initial_migrations()

        self.events = []

        self.events += [
            # turn off all migration from t=0 until T_ACL_GRW (5kya)
            msprime.MigrationRateChange(
                time=0,
                rate=0),

            # stop rapid population growth in AF
            msprime.PopulationParametersChange(
                time=self.T_ACL_GRW,
                initial_size=self.N_AF0,
                growth_rate=0,
                population_id=ids['AF']),

            # stop rapid population growth in EU
            msprime.PopulationParametersChange(
                time=self.T_ACL_GRW,
                growth_rate=self.r_EU1,
                population_id=ids['EU']),

            # stop rapid population growth in AS
            msprime.PopulationParametersChange(
                time=self.T_ACL_GRW,
                growth_rate=self.r_AS1,
                population_id=ids['AS']),


            # start modern migration parameters at 5kya
            msprime.MigrationRateChange(
                time=self.T_ACL_GRW,
                rate=migrations.get('AF_EU', 0),
                matrix_index=(ids['EU'], ids['AF'])),
            msprime.MigrationRateChange(
                time=self.T_ACL_GRW,
                rate=migrations.get('EU_AF', 0),
                matrix_index=(ids['AF'], ids['EU'])),
            msprime.MigrationRateChange(
                time=self.T_ACL_GRW,
                rate=migrations.get('AF_AS', 0),
                matrix_index=(ids['AS'], ids['AF'])),
            msprime.MigrationRateChange(
                time=self.T_ACL_GRW,
                rate=migrations.get('AS_AF', 0),
                matrix_index=(ids['AF'], ids['AS'])),
            msprime.MigrationRateChange(
                time=self.T_ACL_GRW,
                rate=migrations.get('EU_AS', 0),
                matrix_index=(ids['AS'], ids['EU'])),
            msprime.MigrationRateChange(
                time=self.T_ACL_GRW,
                rate=migrations.get('AS_EU', 0),
                matrix_index=(ids['EU'], ids['AS'])),

            # set AS popsize to AS0
            msprime.PopulationParametersChange(
                time=self.T_PULSE2-1,
                initial_size=self.N_AS0,
                growth_rate=0,
                population_id=ids['AS']),

            # Neand1 to EAS pulse of introgression
            msprime.MassMigration(
                time=self.T_PULSE2,
                source=ids['AS'],
                destination=ids['N1'],
                proportion=self.m_PULSE2),


            # set EU popsize at EU0
            msprime.PopulationParametersChange(
                time=self.T_EU_AS-1,
                initial_size=self.N_EU0,
                growth_rate=0,
                population_id=ids['EU']),

            # set AS popsize to AS0
            msprime.PopulationParametersChange(
                time=self.T_EU_AS-1,
                initial_size=self.N_AS0,
                growth_rate=0,
                population_id=ids['AS']),

            # AS merges into EU, now termed "B"
            msprime.MassMigration(
                time=self.T_EU_AS,
                source=ids['AS'],
                destination=ids['EU'],
                proportion=1.0)
        ]

        ids['B'] = ids['EU']
        migrations = self.get_later_migrations()

        self.events += [
            # set all migration rates to zero
            msprime.MigrationRateChange(
                time=self.T_EU_AS,
                rate=0),

            # Change rate m_AF_EU to m_AF_B m(source, dest) backward in time
            msprime.MigrationRateChange(
                time=self.T_EU_AS,
                rate=migrations.get('AF_B', 0),
                matrix_index=(ids['B'], ids['AF'])),
            msprime.MigrationRateChange(
                time=self.T_EU_AS,
                rate=migrations.get('B_AF', 0),
                matrix_index=(ids['AF'], ids['B'])),

            # set parameters of population "B"
            msprime.PopulationParametersChange(
                time=self.T_EU_AS,
                initial_size=self.N_B,
                growth_rate=0,
                population_id=ids['B']),

            # Neand1 to EUR_EAS pulse of introgression
            msprime.MassMigration(
                time=self.T_PULSE1,
                source=ids['B'],
                destination=0,
                proportion=self.m_PULSE1),

            # Population B merges into Africa at T_B
            msprime.MassMigration(
                time=self.T_B,
                source=ids['B'],
                destination=ids['AF'],
                proportion=1.0),

            # set all migration rates to zero
            msprime.MigrationRateChange(
                time=self.T_B,
                rate=0),

            # N_2 merges with N_1 at T_N1_N2
            msprime.MassMigration(
                time=self.T_N1_N2,
                source=ids['N2'],
                destination=ids['N1'],
                proportion=1.0),

            # set parameters of ancestral Neandertal population
            msprime.PopulationParametersChange(
                time=self.T_N1_N2,
                initial_size=self.N_N1,
                population_id=ids['N1']),

            # set parameters of ancestral modern human population
            msprime.PopulationParametersChange(
                time=self.T_AF,
                initial_size=self.N_A,
                population_id=ids['AF']),

            # DE merges with N1
            msprime.MassMigration(
                time=self.T_DE_N,
                source=ids['DE'],
                destination=ids['N1'],
                proportion=1.0),

            msprime.PopulationParametersChange(
                time=self.T_DE_N,
                initial_size=self.N_N1,
                population_id=ids['N1']),

            # Neandertals merge into modern human lineage at time T_MH_N
            msprime.MassMigration(
                time=self.T_MH_N,
                source=ids['N1'],
                destination=ids['AF'],
                proportion=1.0),

            # set parameters of ancetral hominin population
            msprime.PopulationParametersChange(
                time=self.T_MH_N,
                initial_size=self.N_A,
                population_id=ids['AF']),

            # Chimp lineage merges into ancestral hominin population
            # at time T_MH_CH
            msprime.MassMigration(
                time=self.T_MH_CH,
                source=ids['CH'],
                destination=ids['AF'],
                proportion=1.0),

            # set parameters of ancestral hominin population
            msprime.PopulationParametersChange(
                time=self.T_MH_CH,
                initial_size=self.N_A,
                population_id=ids['AF'])
        ]


class Tenn_pulsed_migration(Tenn_demography):
    """Modified version of Tennison with pulsed migration"""
    def set_demographic_events(self):
        """set the list of demographic events in self.events"""
        ids = self.get_population_map()

        # NOTE: this is initial migrations to match legacy code
        # it is more logical to be later migrations.
        # would also want to change migration assignment after B
        migrations = self.get_initial_migrations()

        self.events = []

        self.events += [
            # turn off all migration from t=0 until T_ACL_GRW (5kya)
            msprime.MigrationRateChange(
                time=0,
                rate=0),

            # stop rapid population growth in AF
            msprime.PopulationParametersChange(
                time=self.T_ACL_GRW,
                initial_size=self.N_AF0,
                growth_rate=0,
                population_id=ids['AF']),

            # stop rapid population growth in EU
            msprime.PopulationParametersChange(
                time=self.T_ACL_GRW,
                growth_rate=self.r_EU1,
                population_id=ids['EU']),

            # stop rapid population growth in AS
            msprime.PopulationParametersChange(
                time=self.T_ACL_GRW,
                growth_rate=self.r_AS1,
                population_id=ids['AS']),

            # set AS popsize to AS0
            msprime.PopulationParametersChange(
                time=self.T_PULSE2-1,
                initial_size=self.N_AS0,
                growth_rate=0,
                population_id=ids['AS']),

            # Neand1 to EAS pulse of introgression
            msprime.MassMigration(
                time=self.T_PULSE2,
                source=ids['AS'],
                destination=ids['N1'],
                proportion=self.m_PULSE2),

            # start modern migration parameters at 5kya
            msprime.MigrationRateChange(
                time=self.T_PULSE2,
                rate=(self.T_EU_AS - self.T_ACL_GRW) *
                migrations.get('AF_EU', 0),
                matrix_index=(ids['EU'], ids['AF'])),
            msprime.MigrationRateChange(
                time=self.T_PULSE2,
                rate=(self.T_EU_AS - self.T_ACL_GRW) *
                migrations.get('EU_AF', 0),
                matrix_index=(ids['AF'], ids['EU'])),
            msprime.MigrationRateChange(
                time=self.T_PULSE2,
                rate=(self.T_EU_AS - self.T_ACL_GRW) *
                migrations.get('AF_AS', 0),
                matrix_index=(ids['AS'], ids['AF'])),
            msprime.MigrationRateChange(
                time=self.T_PULSE2,
                rate=(self.T_EU_AS - self.T_ACL_GRW) *
                migrations.get('AS_AF', 0),
                matrix_index=(ids['AF'], ids['AS'])),
            msprime.MigrationRateChange(
                time=self.T_PULSE2,
                rate=(self.T_EU_AS - self.T_ACL_GRW) *
                migrations.get('EU_AS', 0),
                matrix_index=(ids['AS'], ids['EU'])),
            msprime.MigrationRateChange(
                time=self.T_PULSE2,
                rate=(self.T_EU_AS - self.T_ACL_GRW) *
                migrations.get('AS_EU', 0),
                matrix_index=(ids['EU'], ids['AS'])),

            # stop migration parameters
            msprime.MigrationRateChange(
                time=self.T_PULSE2+10,
                rate=0),

            # set EU popsize at EU0
            msprime.PopulationParametersChange(
                time=self.T_EU_AS-1,
                initial_size=self.N_EU0,
                growth_rate=0,
                population_id=ids['EU']),

            # set AS popsize to AS0
            msprime.PopulationParametersChange(
                time=self.T_EU_AS-1,
                initial_size=self.N_AS0,
                growth_rate=0,
                population_id=ids['AS']),

            # AS merges into EU, now termed "B"
            msprime.MassMigration(
                time=self.T_EU_AS,
                source=ids['AS'],
                destination=ids['EU'],
                proportion=1.0),
        ]

        ids['B'] = ids['EU']
        migrations = self.get_later_migrations()

        self.events += [
            # set all migration rates to zero
            msprime.MigrationRateChange(
                time=self.T_EU_AS,
                rate=0),

            # set parameters of population "B"
            msprime.PopulationParametersChange(
                time=self.T_EU_AS,
                initial_size=self.N_B,
                growth_rate=0,
                population_id=ids['B']),


            # start ancient migration parameters at 45kya
            msprime.MigrationRateChange(
                time=self.T_B_BN,
                rate=(self.T_B - self.T_EU_AS) *
                migrations.get('AF_B', 0),
                matrix_index=(ids['B'], ids['AF'])),
            msprime.MigrationRateChange(
                time=self.T_B_BN,
                rate=(self.T_B - self.T_EU_AS) *
                migrations.get('B_AF', 0),
                matrix_index=(ids['AF'], ids['B'])),
            # stop migration parameters
            msprime.MigrationRateChange(
                time=self.T_B_BN+10,
                rate=0),


            # Neand1 to EUR_EAS pulse of introgression
            msprime.MassMigration(
                time=self.T_PULSE1,
                source=ids['B'],
                destination=ids['N1'],
                proportion=self.m_PULSE1),

            # Population B merges into Africa at T_B
            msprime.MassMigration(
                time=self.T_B,
                source=ids['B'],
                destination=ids['AF'],
                proportion=1.0),

            # set all migration rates to zero
            msprime.MigrationRateChange(
                time=self.T_B,
                rate=0),

            # N_2 merges with N_1 at T_N1_N2
            msprime.MassMigration(
                time=self.T_N1_N2,
                source=ids['N2'],
                destination=ids['N1'],
                proportion=1.0),

            # set parameters of ancestral Neandertal population
            msprime.PopulationParametersChange(
                time=self.T_N1_N2,
                initial_size=self.N_N1,
                population_id=ids['N1']),

            # set parameters of ancestral modern human population
            msprime.PopulationParametersChange(
                time=self.T_AF,
                initial_size=self.N_A,
                population_id=ids['AF']),

            # DE merges with N1
            msprime.MassMigration(
                time=self.T_DE_N,
                source=ids['DE'],
                destination=ids['N1'],
                proportion=1.0),

            msprime.PopulationParametersChange(
                time=self.T_DE_N,
                initial_size=self.N_N1,
                population_id=ids['N1']),
            # Neandertals merge into modern human lineage at time T_MH_N
            msprime.MassMigration(
                time=self.T_MH_N,
                source=ids['N1'],
                destination=ids['AF'],
                proportion=1.0),

            # set parameters of ancetral hominin population
            msprime.PopulationParametersChange(
                time=self.T_MH_N,
                initial_size=self.N_A,
                population_id=ids['AF']),

            # Chimp lineage merges into ancestral hominin population
            # at time T_MH_CH
            msprime.MassMigration(
                time=self.T_MH_CH,
                source=ids['CH'],
                destination=ids['AF'],
                proportion=1.0),

            # set parameters of ancestral hominin population
            msprime.PopulationParametersChange(
                time=self.T_MH_CH,
                initial_size=self.N_A,
                population_id=ids['AF'])
        ]


class Sriram_demography(Base_demography):
    def set_constants(self):
        Base_demography.set_constants(self)

        # effective population size used for scaling
        self.N_A = 10000
        # Altai/Vindija lineage effective population size
        self.N_N1 = 10000
        # Eastern Neandertal effective population size
        self.N_N2 = 10000
        # Neandertal population bottleneck size
        self.N_N_BN = 100
        # Population size after EUR + ASN join
        self.N_B = 10000
        # Population size of EUR_ASN bottleneck
        # (following initial admixture event)
        self.N_B_BN = 100
        self.N_AF = 10000
        self.N_CH = 10000
        self.N_DE = 10000
        self.N_EU = 10000
        self.N_AS = 10000

        # Neandertal bottleneck            ### DIFFERENT FROM TENNESSEN
        self.T_N_BN = 150e3 / self.generation_time
        # out of Africa                    ### DIFFERENT FROM TENNESSEN
        self.T_B = 75e3 / self.generation_time
        # Neandertal introgression pulse 1 ### DIFFERENT FROM TENNESSEN
        self.T_PULSE1 = 50e3 / self.generation_time
        # A 3rd Bottleneck following shortly after initial admixture event
        # SET AT 1800 GEN
        self.T_B_BN = 45e3 / self.generation_time
        # European / East Asian split
        self.T_EU_AS = 23e3 / self.generation_time

    def set_populations(self):
        """set the self.populations variable
        a list of Populations for the simulation"""

        Base_demography.set_populations(self)

        # only change from tennison is all growth rates are 0
        for pop in self.populations:
            pop.rate = 0.0

    def set_demographic_events(self):
        ids = self.get_population_map()
        migrations = self.get_later_migrations()

        self.events = []

        self.events += [
            # Neand1 to EAS pulse of introgression
            msprime.MassMigration(
                time=self.T_PULSE2,
                source=ids['AS'],
                destination=ids['N1'],
                proportion=self.m_PULSE2),

            # AS merges into EU, now termed "B"
            msprime.MassMigration(
                time=self.T_EU_AS,
                source=ids['AS'],
                destination=ids['EU'],
                proportion=1.0)
        ]

        ids['B'] = ids['EU']

        self.events += [
            # set all migration rates to zero
            msprime.MigrationRateChange(
                time=self.T_EU_AS,
                rate=0),

            # migration between "B" and Africa begins
            msprime.MigrationRateChange(
                time=self.T_EU_AS,
                rate=migrations.get('AF_B', 0),
                matrix_index=(ids['B'], ids['AF'])),
            msprime.MigrationRateChange(
                time=self.T_EU_AS,
                rate=migrations.get('B_AF', 0),
                matrix_index=(ids['AF'], ids['B'])),

            # set parameters of population "B"
            msprime.PopulationParametersChange(
                time=self.T_EU_AS,
                initial_size=self.N_B,
                growth_rate=0,
                population_id=ids['B']),

            # population bottleneck of B begins
            msprime.PopulationParametersChange(
                time=self.T_B_BN,
                initial_size=self.N_B_BN,
                growth_rate=0,
                population_id=ids['B']),

            # population bottleneck of B ends shortly before initial admixture
            msprime.PopulationParametersChange(
                time=self.T_B_BN + 20,
                initial_size=self.N_B,
                growth_rate=0,
                population_id=ids['B']),

            # Neand1 to EUR_EAS pulse of introgression
            msprime.MassMigration(
                time=self.T_PULSE1,
                source=ids['B'],
                destination=ids['N1'],
                proportion=self.m_PULSE1),

            # Population B merges into Africa at T_B
            msprime.MassMigration(
                time=self.T_B,
                source=ids['B'],
                destination=ids['AF'],
                proportion=1.0),

            # set all migration rates to zero
            msprime.MigrationRateChange(
                time=self.T_B,
                rate=0),

            # N_2 merges with N_1 at T_N1_N2
            msprime.MassMigration(
                time=self.T_N1_N2,
                source=ids['N2'],
                destination=ids['N1'],
                proportion=1.0),

            # set parameters of ancestral Neandertal population
            msprime.PopulationParametersChange(
                time=self.T_N1_N2,
                initial_size=self.N_A,
                population_id=ids['N1']),

            # Neand1 has short bottleneck
            msprime.PopulationParametersChange(
                time=self.T_N_BN,
                initial_size=self.N_N_BN,
                growth_rate=0,
                population_id=ids['N1']),

            # Neand1 bottleneck ends
            msprime.PopulationParametersChange(
                time=self.T_N_BN + 20,
                initial_size=self.N_A,
                growth_rate=0,
                population_id=ids['N1']),

            # DE merges with N1
            msprime.MassMigration(
                time=self.T_DE_N,
                source=ids['DE'],
                destination=ids['N1'],
                proportion=1.0),
            msprime.PopulationParametersChange(
                time=self.T_DE_N,
                initial_size=self.N_A,
                population_id=ids['N1']),

            # Neandertals merge into modern human lineage at time T_MH_N
            msprime.MassMigration(
                time=self.T_MH_N,
                source=ids['N1'],
                destination=ids['AF'],
                proportion=1.0),

            # set parameters of ancetral hominin population
            msprime.PopulationParametersChange(
                time=self.T_MH_N,
                initial_size=self.N_A,
                population_id=ids['AF']),

            # Chimp lineage merges into ancestral hominin population
            # at time T_MH_CH
            msprime.MassMigration(
                time=self.T_MH_CH,
                source=ids['CH'],
                destination=ids['AF'],
                proportion=1.0),

            # set parameters of ancestral hominin population
            msprime.PopulationParametersChange(
                time=self.T_MH_CH,
                initial_size=self.N_A,
                population_id=ids['AF'])
        ]


class SplitPop_demography(Base_demography):
    def __init__(self, options):
        # add in new options for split population
        self.N_SP = options.split_population_size
        self.f_SP = options.split_population_proportion
        # use all normal options
        Base_demography.__init__(self, options)

    def set_constants(self):
        Base_demography.set_constants(self)

        # split population off of B
        self.T_SP = 54e3 / self.generation_time

        # rapid population growth at 5.1kya
        self.r_EU2 = 0.0195
        # rapid population growth at 5.1kya
        self.r_AS2 = 0.025
        # rapid popualtion growth at 5.1kya
        self.r_AF2 = 0.0166

        ### TESTING FST
        # Updae to Gravel 2011 Values
        self.r_EU1 = 0.0038
        # self.r_AS1 = 0.0031
        # self.r_EU1 = 0.0027
        # Gravel NIEHS vales
        # self.r_EU1 = 0.0040
        # self.r_AS1 = 0.0055

        # Update to Gravel 2011 Values ; Exons
        self.N_AF0 = 15388
        self.N_EU0 = 1620
        self.N_AS0 = 821
        self.N_B = 2758
        ## No effect
        # self.T_B = 60e3 / self.generation_time
        ## Keep These values
        self.T_AF = 316e3 / self.generation_time
        self.T_B = 98e3 / self.generation_time
        self.T_EU_AS = 28e3 / self.generation_time
        ## Increasing mutation rate dirves up number of variants
            ## but doesn't change Fst values much
        self.mutation_rate = 2.5e-8
        ## Decreasing Ne drize up the Fst value...
        # self.N_EU = 100e3
        # self.N_AS = 100e3
        # self.N_AF = 100e3
        ## Increase Ne drove down Fst values...pretty close to empirical
            ## AFR-EUR  AFR-EAS  EUR-EAS
            ## 0.157105 0.177657 0.0408926
        self.N_EU = 512e4
        self.N_AS = 640e4
        self.N_AF = 424e4

    def set_populations(self):
        self.populations = [
            Population(abbreviation='N1',
                       population_size=self.N_N1,
                       samples=self.S_N1,
                       generations=self.T_N1_SAMPLE,
                       long_name='Neand1'),

            Population(abbreviation='N2',
                       population_size=self.N_N2,
                       samples=self.S_N2,
                       generations=self.T_N2_SAMPLE,
                       long_name='Neand2'),

            Population(abbreviation='AF',
                       population_size=self.N_AF,
                       growth_rate=self.r_AF2,
                       samples=self.options.AF_sample_size,
                       long_name='AFR'),

            Population(abbreviation='EU',
                       population_size=self.N_EU,
                       growth_rate=self.r_EU2,
                       samples=self.options.EU_sample_size,
                       long_name='EUR'),

            Population(abbreviation='AS',
                       population_size=self.N_AS,
                       growth_rate=self.r_AS2,
                       samples=self.options.AS_sample_size,
                       long_name='ASN'),

            Population(abbreviation='CH',
                       population_size=self.N_CH,
                       long_name='Chimp'),

            Population(abbreviation='DE',
                       population_size=self.N_DE,
                       generations=50e3 / self.generation_time,
                       long_name='Deni'),

            Population(abbreviation='SP',
                       population_size=self.N_SP,
                       generations=50e3 / self.generation_time,
                       long_name='Split')

        ]

    def get_samples(self):
        """need to override default to remove last (SP) population"""
        result = []
        for i, p in enumerate(self.populations[:-1]):
            result += p.get_sample(i)
        return result

    def get_sample_indices(self):
        result = []
        for i, p in enumerate(self.populations[:-1]):
            result += p.get_indices(i)
        return result

    def set_demographic_events(self):
        ids = self.get_population_map()
        migrations = self.get_later_migrations()

        self.events = []

        self.events += [
            # stop rapid population growth in AF
            msprime.PopulationParametersChange(
                time=self.T_ACL_GRW,
                initial_size=self.N_AF0,
                growth_rate=0,
                population_id=ids['AF']),

            # stop rapid population growth in EU
            msprime.PopulationParametersChange(
                time=self.T_ACL_GRW,
                growth_rate=self.r_EU1,
                population_id=ids['EU']),

            # stop rapid population growth in AS
            msprime.PopulationParametersChange(
                time=self.T_ACL_GRW,
                growth_rate=self.r_AS1,
                population_id=ids['AS']),

            # # set AS popsize to AS0
            # msprime.PopulationParametersChange(
            #     time=self.T_PULSE2-1,
            #     initial_size=self.N_AS0,
            #     growth_rate=0,
            #     population_id=ids['AS']),

            # Neand1 to EAS pulse of introgression
            msprime.MassMigration(
                time=self.T_PULSE2,
                source=ids['AS'],
                destination=ids['N1'],
                proportion=self.m_PULSE2),

            # set EU popsize at EU0
            msprime.PopulationParametersChange(
                time=self.T_EU_AS-1,
                initial_size=self.N_EU0,
                growth_rate=0,
                population_id=ids['EU']),

            # set AS popsize to AS0
            msprime.PopulationParametersChange(
                time=self.T_EU_AS-1,
                initial_size=self.N_AS0,
                growth_rate=0,
                population_id=ids['AS']),

            # AS merges into EU, now termed "B"
            msprime.MassMigration(
                time=self.T_EU_AS,
                source=ids['AS'],
                destination=ids['EU'],
                proportion=1.0)
        ]

        ids['B'] = ids['EU']

        self.events += [
            # set all migration rates to zero
            msprime.MigrationRateChange(
                time=self.T_EU_AS,
                rate=0),

            # migration between "B" and Africa begins
            msprime.MigrationRateChange(
                time=self.T_EU_AS,
                rate=migrations.get('AF_B', 0.0),
                matrix_index=(ids['B'], ids['AF'])),
            msprime.MigrationRateChange(
                time=self.T_EU_AS,
                rate=migrations.get('B_AF', 0.0),
                matrix_index=(ids['AF'], ids['B'])),

            # set parameters of population "B"
            msprime.PopulationParametersChange(
                time=self.T_EU_AS,
                initial_size=self.N_B,
                growth_rate=0,
                population_id=ids['B']),

            # Population splits off of EUR_EAS, 10% migration rate to SplitPop
            msprime.MassMigration(
                time=self.T_SP,
                source=ids['B'],
                destination=ids['SP'],
                proportion=self.f_SP),

            # SplitPop is set at size Ne=100
            msprime.PopulationParametersChange(
                time=self.T_SP,
                initial_size=self.N_SP,
                growth_rate=0,
                population_id=ids['SP']),

            # Neand1 to SplitPop pulse of introgression
            msprime.MassMigration(
                time=self.T_PULSE1,
                source=ids['SP'],
                destination=ids['N1'],
                proportion=self.m_PULSE1),

            # SplitPop merges into EUR_EAS at 1kya before admixture
            msprime.MassMigration(
                time=self.T_PULSE1+40,
                source=ids['SP'],
                destination=ids['B'],
                proportion=1.0),

            # Population B merges into African at T_B
            msprime.MassMigration(
                time=self.T_B,
                source=ids['B'],
                destination=ids['AF'],
                proportion=1.0),

            # set all migration rates to zero
            msprime.MigrationRateChange(
                time=self.T_B,
                rate=0),

            # N_2 merges with N_1 at T_N1_N2
            msprime.MassMigration(
                time=self.T_N1_N2,
                source=ids['N2'],
                destination=ids['N1'],
                proportion=1.0),

            # set parameters of ancestral Neandertal population
            msprime.PopulationParametersChange(
                time=self.T_N1_N2,
                initial_size=self.N_N1,
                population_id=ids['N1']),

            # set parameters of ancestral modern human population
            msprime.PopulationParametersChange(
                time=self.T_AF,
                initial_size=self.N_A,
                population_id=ids['AF']),

            # DE merges with N1
            msprime.MassMigration(
                time=self.T_DE_N,
                source=ids['DE'],
                destination=ids['N1'],
                proportion=1.0),
            msprime.PopulationParametersChange(
                time=self.T_DE_N,
                initial_size=self.N_N1,
                population_id=ids['N1']),

            # Neandertals merge into modern human lineage at time T_MH_N
            msprime.MassMigration(
                time=self.T_MH_N,
                source=ids['N1'],
                destination=ids['AF'],
                proportion=1.0),

            # set parameters of ancetral hominin population
            msprime.PopulationParametersChange(
                time=self.T_MH_N,
                initial_size=self.N_A,
                population_id=ids['AF']),

            # Chimp lineage merges into ancestral hominin population
            # at time T_MH_CH
            msprime.MassMigration(
                time=self.T_MH_CH,
                source=ids['CH'],
                destination=ids['AF'],
                proportion=1.0),

            # set parameters of ancestral hominin population
            msprime.PopulationParametersChange(
                time=self.T_MH_CH,
                initial_size=self.N_A,
                population_id=ids['AF'])
        ]

    def get_debug_configuration(self):
        return [pop.get_debug_configuration(includeRate=True)
                for pop in self.populations]


class Out_of_africa_demography(Base_demography):
    def simulate(self, replicates):
        raise NotImplementedError

    def set_constants(self):
        Base_demography.set_constants(self)

        self.N_A = 7300
        self.N_B = 2100
        self.N_AF = 12300
        self.N_oAF = 1000
        self.N_EU0 = 1000
        self.N_AS0 = 510

        self.T_AF = 220e3 / self.generation_time
        # create a split between and African outgroup and main African group
        self.T_oAF_AF = 200e3 / self.generation_time
        self.T_B = 140e3 / self.generation_time
        self.T_EU_AS = 21.2e3 / self.generation_time
        # We need to work out the starting (diploid) population sizes based on
        # the growth rates provided for these two populations
        self.r_EU = 0.004
        self.r_AS = 0.0055
        self.N_EU = self.N_EU0 / math.exp(-self.r_EU * self.T_EU_AS)
        self.N_AS = self.N_AS0 / math.exp(-self.r_AS * self.T_EU_AS)

    def set_populations(self):
        self.populations = [
            Population(abbreviation='oAF',
                       population_size=self.N_oAF,
                       long_name='oAFR'),

            Population(abbreviation='AF',
                       population_size=self.N_AF,
                       samples=self.options.AF_sample_size,
                       long_name='AFR'),

            Population(abbreviation='EU',
                       population_size=self.N_EU,
                       growth_rate=self.r_EU,
                       samples=self.options.EU_sample_size,
                       long_name='EUR'),

            Population(abbreviation='AS',
                       population_size=self.N_AS,
                       growth_rate=self.r_AS,
                       samples=self.options.AS_sample_size,
                       long_name='ASN'),

            Population(abbreviation='CH',
                       population_size=self.N_CH,
                       long_name='Chimp'),
        ]

    def set_demographic_events(self):
        ids = self.get_population_map()
        migrations = self.get_later_migrations()

        self.events = []

        self.events += [

            # CEU and CHB merge into B with rate changes at T_EU_AS
            msprime.MassMigration(
                time=self.T_EU_AS,
                source=ids['AS'],
                destination=ids['EU'],
                proportion=1.0),

            msprime.MigrationRateChange(
                time=self.T_EU_AS,
                rate=0),

            msprime.MigrationRateChange(
                time=self.T_EU_AS,
                rate=migrations.get('AF_B', 0),
                matrix_index=(ids['EU'], ids['AF'])),

            msprime.MigrationRateChange(
                time=self.T_EU_AS,
                rate=migrations.get('B_AF', 0),
                matrix_index=(ids['AF'], ids['EU'])),

            msprime.PopulationParametersChange(
                time=self.T_EU_AS,
                initial_size=self.N_B,
                growth_rate=0,
                population_id=ids['EU'])
        ]

        ids['B'] = ids['EU']

        self.events += [
            # Population B merges into YRI at T_B
            msprime.MassMigration(
                time=self.T_B,
                source=ids['B'],
                destination=ids['AF'],
                proportion=1.0),

            msprime.MassMigration(
                time=self.T_oAF_AF,
                source=ids['oAF'],
                destination=ids['AF'],
                proportion=1.0),

            # Size changes to N_A at T_AF
            msprime.PopulationParametersChange(
                time=self.T_AF,
                initial_size=self.N_A,
                population_id=ids['AF'])
        ]

    def get_debug_configuration(self):
        return [pop.get_debug_configuration(includeRate=True)
                for pop in self.populations]


class Pre_out_of_africa_admixture(Sriram_demography):
    def set_constants(self):
        Sriram_demography.set_constants(self)

        self.T_PULSE_preOOA = 250e3 / self.generation_time


    def set_populations(self):
        """set the self.populations variable
        a list of Populations for the simulation"""

        Sriram_demography.set_populations(self)

        # only change from tennison is all growth rates are 0
        for pop in self.populations:
            pop.rate = 0.0

    def set_demographic_events(self):
        ids = self.get_population_map()
        migrations = self.get_later_migrations()

        self.events = []

        self.events += [
            # Neand1 to EAS pulse of introgression
            # msprime.MassMigration(
            #     time=self.T_PULSE2,
            #     source=ids['AS'],
            #     destination=ids['N1'],
            #     proportion=self.m_PULSE2),

            # AS merges into EU, now termed "B"
            msprime.MassMigration(
                time=self.T_EU_AS,
                source=ids['AS'],
                destination=ids['EU'],
                proportion=1.0)
        ]

        ids['B'] = ids['EU']

        self.events += [
            # set all migration rates to zero
            msprime.MigrationRateChange(
                time=self.T_EU_AS,
                rate=0),

            # migration between "B" and Africa begins
            msprime.MigrationRateChange(
                time=self.T_EU_AS,
                rate=migrations.get('AF_B', 0),
                matrix_index=(ids['B'], ids['AF'])),
            msprime.MigrationRateChange(
                time=self.T_EU_AS,
                rate=migrations.get('B_AF', 0),
                matrix_index=(ids['AF'], ids['B'])),

            # set parameters of population "B"
            msprime.PopulationParametersChange(
                time=self.T_EU_AS,
                initial_size=self.N_B,
                growth_rate=0,
                population_id=ids['B']),

            # population bottleneck of B begins
            msprime.PopulationParametersChange(
                time=self.T_B_BN,
                initial_size=self.N_B_BN,
                growth_rate=0,
                population_id=ids['B']),

            # population bottleneck of B ends shortly before initial admixture
            msprime.PopulationParametersChange(
                time=self.T_B_BN + 20,
                initial_size=self.N_B,
                growth_rate=0,
                population_id=ids['B']),

            # Neand1 to EUR_EAS pulse of introgression
            msprime.MassMigration(
                time=self.T_PULSE1,
                source=ids['B'],
                destination=ids['N1'],
                proportion=self.m_PULSE1),

            # Population B merges into Africa at T_B
            msprime.MassMigration(
                time=self.T_B,
                source=ids['B'],
                destination=ids['AF'],
                proportion=1.0),

            # set all migration rates to zero
            msprime.MigrationRateChange(
                time=self.T_B,
                rate=0),

            # N_2 merges with N_1 at T_N1_N2
            msprime.MassMigration(
                time=self.T_N1_N2,
                source=ids['N2'],
                destination=ids['N1'],
                proportion=1.0),

            # set parameters of ancestral Neandertal population
            msprime.PopulationParametersChange(
                time=self.T_N1_N2,
                initial_size=self.N_A,
                population_id=ids['N1']),

            # Neand1 has short bottleneck
            msprime.PopulationParametersChange(
                time=self.T_N_BN,
                initial_size=self.N_N_BN,
                growth_rate=0,
                population_id=ids['N1']),

            # Neand1 bottleneck ends
            msprime.PopulationParametersChange(
                time=self.T_N_BN + 20,
                initial_size=self.N_A,
                growth_rate=0,
                population_id=ids['N1']),

            # Human --> Neand admixture
            msprime.MassMigration(
                time=self.T_PULSE_preOOA,
                source=ids['N1'],
                destination=ids['AF'],
                proportion=self.m_PULSE2),

            # DE merges with N1
            msprime.MassMigration(
                time=self.T_DE_N,
                source=ids['DE'],
                destination=ids['N1'],
                proportion=1.0),
            msprime.PopulationParametersChange(
                time=self.T_DE_N,
                initial_size=self.N_A,
                population_id=ids['N1']),

            # Neandertals merge into modern human lineage at time T_MH_N
            msprime.MassMigration(
                time=self.T_MH_N,
                source=ids['N1'],
                destination=ids['AF'],
                proportion=1.0),

            # set parameters of ancetral hominin population
            msprime.PopulationParametersChange(
                time=self.T_MH_N,
                initial_size=self.N_A,
                population_id=ids['AF']),

            # Chimp lineage merges into ancestral hominin population
            # at time T_MH_CH
            msprime.MassMigration(
                time=self.T_MH_CH,
                source=ids['CH'],
                destination=ids['AF'],
                proportion=1.0),

            # set parameters of ancestral hominin population
            msprime.PopulationParametersChange(
                time=self.T_MH_CH,
                initial_size=self.N_A,
                population_id=ids['AF'])
        ]


class AFR_split_F4_demography(Base_demography):
    def set_constants(self):
        Sriram_demography.set_constants(self)

        self.T_PULSE_preOOA = 250e3 / self.generation_time

        # effective population size used for scaling
        self.N_A = 10000
        # Altai/Vindija lineage effective population size
        self.N_N1 = 10000
        # Eastern Neandertal effective population size
        self.N_N2 = 10000
        # Neandertal population bottleneck size
        self.N_N_BN = 100
        # Population size after EUR + ASN join
        self.N_B = 10000
        # Population size of EUR_ASN bottleneck
        # (following initial admixture event)
        self.N_B_BN = 100
        self.N_AF = 10000
        self.N_AF2 = 10000
        self.N_CH = 10000
        self.N_DE = 10000
        self.N_EU = 10000
        self.N_AS = 10000

        # Neandertal bottleneck            ### DIFFERENT FROM TENNESSEN
        self.T_N_BN = 150e3 / self.generation_time
        # out of Africa                    ### DIFFERENT FROM TENNESSEN
        self.T_B = 75e3 / self.generation_time
        # Neandertal introgression pulse 1 ### DIFFERENT FROM TENNESSEN
        self.T_PULSE1 = 50e3 / self.generation_time
        # A 3rd Bottleneck following shortly after initial admixture event
        # SET AT 1800 GEN
        self.T_B_BN = 45e3 / self.generation_time
        # European / East Asian split
        self.T_EU_AS = 23e3 / self.generation_time

    def set_populations(self):
        self.populations = [
            Population(abbreviation='N1',
                       population_size=self.N_N1,
                       samples=self.S_N1,
                       generations=self.T_N1_SAMPLE,
                       long_name='Neand1'),

            Population(abbreviation='N2',
                       population_size=self.N_N2,
                       samples=self.S_N2,
                       generations=self.T_N2_SAMPLE,
                       long_name='Neand2'),

            Population(abbreviation='AF',
                       population_size=self.N_AF,
                       growth_rate=0,
                       samples=self.options.AF_sample_size,
                       long_name='AFR1'),

            Population(abbreviation='AF2',
                       population_size=self.N_AF2,
                       growth_rate=0,
                       samples=self.options.AF_sample_size,
                       long_name='AFR2'),

            Population(abbreviation='EU',
                       population_size=self.N_EU,
                       growth_rate=0,
                       samples=self.options.EU_sample_size,
                       long_name='EUR'),

            Population(abbreviation='AS',
                       population_size=self.N_AS,
                       growth_rate=0,
                       samples=self.options.AS_sample_size,
                       long_name='ASN'),

            Population(abbreviation='CH',
                       population_size=self.N_CH,
                       long_name='Chimp'),

            Population(abbreviation='DE',
                       population_size=self.N_DE,
                       generations=50e3 / self.generation_time,
                       long_name='Deni')
        ]

    def set_demographic_events(self):
        ids = self.get_population_map()
        migrations = self.get_later_migrations()

        self.events = []

        self.events += [
            # Neand1 to EAS pulse of introgression
            # msprime.MassMigration(
            #     time=self.T_PULSE2,
            #     source=ids['AS'],
            #     destination=ids['N1'],
            #     proportion=self.m_PULSE2),

            # AS merges into EU, now termed "B"
            msprime.MassMigration(
                time=self.T_EU_AS,
                source=ids['AS'],
                destination=ids['EU'],
                proportion=1.0),

            # Merge Split AF2 into AF
            msprime.MassMigration(
                time=self.T_EU_AS,
                source=ids['AF2'],
                destination=ids['AF'],
                proportion=1.0)
        ]

        ids['B'] = ids['EU']

        self.events += [
            # set all migration rates to zero
            msprime.MigrationRateChange(
                time=self.T_EU_AS,
                rate=0),

            # migration between "B" and Africa begins
            msprime.MigrationRateChange(
                time=self.T_EU_AS,
                rate=migrations.get('AF_B', 0),
                matrix_index=(ids['B'], ids['AF'])),
            msprime.MigrationRateChange(
                time=self.T_EU_AS,
                rate=migrations.get('B_AF', 0),
                matrix_index=(ids['AF'], ids['B'])),

            # set parameters of population "B"
            msprime.PopulationParametersChange(
                time=self.T_EU_AS,
                initial_size=self.N_B,
                growth_rate=0,
                population_id=ids['B']),

            # population bottleneck of B begins
            msprime.PopulationParametersChange(
                time=self.T_B_BN,
                initial_size=self.N_B_BN,
                growth_rate=0,
                population_id=ids['B']),

            # population bottleneck of B ends shortly before initial admixture
            msprime.PopulationParametersChange(
                time=self.T_B_BN + 20,
                initial_size=self.N_B,
                growth_rate=0,
                population_id=ids['B']),

            # Neand1 to EUR_EAS pulse of introgression
            msprime.MassMigration(
                time=self.T_PULSE1,
                source=ids['B'],
                destination=ids['N1'],
                proportion=self.m_PULSE1),

            # Population B merges into Africa at T_B
            msprime.MassMigration(
                time=self.T_B,
                source=ids['B'],
                destination=ids['AF'],
                proportion=1.0),

            # set all migration rates to zero
            msprime.MigrationRateChange(
                time=self.T_B,
                rate=0),

            # N_2 merges with N_1 at T_N1_N2
            msprime.MassMigration(
                time=self.T_N1_N2,
                source=ids['N2'],
                destination=ids['N1'],
                proportion=1.0),

            # set parameters of ancestral Neandertal population
            msprime.PopulationParametersChange(
                time=self.T_N1_N2,
                initial_size=self.N_A,
                population_id=ids['N1']),

            # Neand1 has short bottleneck
            msprime.PopulationParametersChange(
                time=self.T_N_BN,
                initial_size=self.N_N_BN,
                growth_rate=0,
                population_id=ids['N1']),

            # Neand1 bottleneck ends
            msprime.PopulationParametersChange(
                time=self.T_N_BN + 20,
                initial_size=self.N_A,
                growth_rate=0,
                population_id=ids['N1']),

            # Human --> Neand admixture
            msprime.MassMigration(
                time=self.T_PULSE_preOOA,
                source=ids['N1'],
                destination=ids['AF'],
                proportion=self.m_PULSE2),

            # DE merges with N1
            msprime.MassMigration(
                time=self.T_DE_N,
                source=ids['DE'],
                destination=ids['N1'],
                proportion=1.0),
            msprime.PopulationParametersChange(
                time=self.T_DE_N,
                initial_size=self.N_A,
                population_id=ids['N1']),

            # Neandertals merge into modern human lineage at time T_MH_N
            msprime.MassMigration(
                time=self.T_MH_N,
                source=ids['N1'],
                destination=ids['AF'],
                proportion=1.0),

            # set parameters of ancetral hominin population
            msprime.PopulationParametersChange(
                time=self.T_MH_N,
                initial_size=self.N_A,
                population_id=ids['AF']),

            # Chimp lineage merges into ancestral hominin population
            # at time T_MH_CH
            msprime.MassMigration(
                time=self.T_MH_CH,
                source=ids['CH'],
                destination=ids['AF'],
                proportion=1.0),

            # set parameters of ancestral hominin population
            msprime.PopulationParametersChange(
                time=self.T_MH_CH,
                initial_size=self.N_A,
                population_id=ids['AF'])
        ]

    def get_debug_configuration(self):
        return [pop.get_debug_configuration(includeRate=True)
                for pop in self.populations]
