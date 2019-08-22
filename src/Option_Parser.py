import sys
from argparse import ArgumentParser


class admixture_option_parser(ArgumentParser):
    def __init__(self):
        description = "Simulate admixture"
        ArgumentParser.__init__(self, description=description)

        self.add_argument("-p", "--population",
                          dest="pop",
                          default="nonAfr",
                          choices=['EAS', 'EUR', 'nonAfr', 'AFR', 'modHum'],
                          help="call introg_haplotypes in AFR, EUR, EAS, "
                               "nonAfr, or all modern humans; "
                               "default=nonAfr")

        self.add_argument("-m", "--model",
                          dest="model",
                          default="Tenn",
                          choices=['Tenn', 'Sriram', 'SplitPop', 'OutOfAFR',
                                   'Tenn_nomod', 'Tenn_pulsed', 'PreOutOfAFR','AFRf4'],
                          help="specify which demographic model to use by "
                               "name; default=Tenn")

        self.add_argument("-s", "--seed",
                          type=int,
                          dest="seed",
                          default=1,
                          help="Set random seed for replicate chromosomes; "
                               "default=1")

        self.add_argument("-n", "--neand1_admixture_proportion",
                          type=float,
                          dest="n1_admix_prop",
                          default=0.02,
                          help="Set N1 admixture proportion; default=0.02")

        self.add_argument("-d", "--deni_admixture_propotion",
                          type=float,
                          dest="n2_admix_prop",
                          default=0.0,
                          help="Set N2 admixture proportion; default=0.0")

        self.add_argument("--Neand1_sample_size",
                          type=int,
                          dest="s_n1",
                          default=2,
                          help="Set Neand1 sample size(haplotypes; default=2")

        self.add_argument("--Neand2_sample_size",
                          type=int,
                          dest="s_n2",
                          default=2,
                          help="Set Neand2 sample size(haplotypes; default=2")

        self.add_argument("-t", "--time_N1_N2_split",
                          type=int,
                          dest="t_n1_n2",
                          default=145,
                          help="Set N1 N2 split time in kya; default=145")

        self.add_argument("--time_N1_sample",
                          type=float,
                          dest="t_n1_sample",
                          default=55,
                          help="set N1 (Vindija) sample time in kya; "
                               "default=55")

        self.add_argument("--time_N2_sample",
                          type=float,
                          dest="t_n2_sample",
                          default=125,
                          help="set N2 (Altai) sample time in kya; "
                               "default=125")

        self.add_argument("-l", "--length_chromosome",
                          type=float,
                          dest="length",
                          default=1e6,
                          help="Define length of simulated chromosome ; "
                               "default=1e6")

        self.add_argument("-e", "--european_sample_size",
                          type=int,
                          dest="EU_sample_size",
                          default=1006,
                          help="Set EU haploid sample size ; default=1006")

        self.add_argument("-a", "--asian_sample_size",
                          type=int,
                          dest="AS_sample_size",
                          default=1008,
                          help="Set AS haploid sample size ; default=1008")

        self.add_argument("-r", "--reference_sample_size",
                          type=int,
                          dest="AF_sample_size",
                          default=2,
                          help="Set AF haploid sample size ; default=2")

        self.add_argument("--split-population-size",
                          type=int,
                          dest="split_population_size",
                          default=100,
                          help="Set split population size in SplitPop model; "
                               "default=100")

        self.add_argument("--split-population-proportion",
                          type=float,
                          dest="split_population_proportion",
                          default=0.1,
                          help="Set split population proportion in SplitPop "
                               "model; default=0.1")

        self.add_argument("-g", "--initial-migration",
                          dest="initial_migrations",
                          action='append',
                          help="Add one or more migration values to initial "
                               "migration matrix, input as POP1_POP2_## to "
                               "set migration from POP1 to POP1 as ## "
                               "(e.g. AF_EU_15e-5 sets migration from "
                               "AF to EU as 15e-5) ;"
                               "default as symmetric migrations of "
                               "AF_EU_2.5e-5, AF_AS_0.78e-5, EU_AS_3.11e-5. "
                               "Setting a value overwrites only that default. "
                               "Only the last repeated population is retained")

        self.add_argument("-G", "--later-migration",
                          dest="later_migrations",
                          action='append',
                          help="Add one or more migration values to "
                               "migrations occuring during demographic "
                               "events. Same format as initial migrations ;"
                               "default symmetric migration of AF_B_15e-5. "
                               "Setting a value overwrites only that default.")

        # Note for the following 4 options:
        # If the flag is not set at all, variable will be None
        # If the flag is set without an argument,
        # variable will be '*' (non valid filename)
        # Otherwise will be the argument provided
        self.add_argument("--debug",
                          dest="debug_file",
                          nargs='?',
                          const='*',
                          default=None,
                          help="debug information filename.  If no argument "
                               "is given or no other filenames are provided, "
                               "will output to stdout.")

        self.add_argument("--haplo",
                          dest="haplo_file",
                          nargs='?',
                          const='*',
                          default=None,
                          help="haplotype call output filename.  If no "
                               "argument is given will output to stdout. If "
                               "not flagged, output will be suppressed.")

        self.add_argument("--ils",
                          dest="ils_file",
                          nargs='?',
                          const='*',
                          default=None,
                          help="ils call output filename.  If no "
                               "argument is given will output to stdout. If "
                               "not flagged, output will be suppressed.")

        self.add_argument("--options",
                          dest="option_file",
                          nargs='?',
                          const='*',
                          default=None,
                          help="Option information filename.  If no argument "
                               "is given will output to stdout. If not "
                               "flagged, output will be suppressed.")

        self.add_argument("--pi",
                          dest="pi_file",
                          nargs='?',
                          const='*',
                          default=None,
                          help="Pi statistic filename.  If no argument "
                               "is given will output to stdout. If not "
                               "flagged, output will be suppressed.")

        # if these are set a filename must be provided as
        # multiple outputs are generated
        self.add_argument("--vcf",
                          dest="vcf_file",
                          default=None,
                          help="vcf output base filename.  If argument is "
                               "BASE will generate BASE.vcf.gz")

        self.add_argument("--popfile",
                          dest="popfile_file",
                          default=None,
                          help="popfile output base filename.  If argument is "
                               "BASE will generate BASE.popfile")

        self.add_argument("--f4dstat",
                          dest="f4dstat_file",
                          default=None,
                          help="F4Dstat output base filename.  If argument "
                               "is BASE will generate eigenstratgeno.BASE.gz, "
                               "snp.BASE.gz, ind.BASE.gz, and "
                               "parfile.F4stat.BASE.gz")

        self.add_argument("--out-dir",
                          dest="out_dir",
                          default=None,
                          help="Base output directory.  When specified "
                               "other output options, all "
                               "output files are generated into the specified "
                               "directory using legacy formatting.")

    def parse_args(self, args=None, namespace=None):
        options = ArgumentParser.parse_args(self, args, namespace)

        default_initial_migrations = [
                'AF_EU_2.5e-5', 'EU_AF_2.5e-5',
                'AF_AS_0.78e-5', 'AS_AF_0.78e-5',
                'EU_AS_3.11e-5', 'AS_EU_3.11e-5']

        default_later_migrations = [
                'AF_B_15e-5', 'B_AF_15e-5']

        # override defaults, retain only last set values
        migration_dict = {}
        for m in default_initial_migrations:
            toks = m.split('_')
            migration_dict['_'.join(toks[0:2])] = toks[2]
        if options.initial_migrations is not None:
            for m in options.initial_migrations:
                toks = m.split('_')
                migration_dict['_'.join(toks[0:2])] = toks[2]

        options.initial_migrations = [
            k + '_' + v for k, v in migration_dict.items()]

        migration_dict = {}
        for m in default_later_migrations:
            toks = m.split('_')
            migration_dict['_'.join(toks[0:2])] = toks[2]
        if options.later_migrations is not None:
            for m in options.later_migrations:
                toks = m.split('_')
                migration_dict['_'.join(toks[0:2])] = toks[2]

        options.later_migrations = [
            k + '_' + v for k, v in migration_dict.items()]

        return options

    # taken from https://stackoverflow.com/questions/5943249 for testing
    def _get_action_from_name(self, name):
        """Given a name, get the Action instance registered with this parser.
        If only it were made available in the ArgumentError object. It is
        passed as it's first arg...
        """
        container = self._actions
        if name is None:
            return None
        for action in container:
            if '/'.join(action.option_strings) == name:
                return action
            elif action.metavar == name:
                return action
            elif action.dest == name:
                return action

    def error(self, message):
        exc = sys.exc_info()[1]
        if exc:
            exc.argument = self._get_action_from_name(exc.argument_name)
            raise exc
        ArgumentParser.error(self, message)
