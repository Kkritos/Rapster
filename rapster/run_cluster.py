'''
 Copyright (C) 2026  Konstantinos Kritos <kkritos1@jhu.edu>

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <https://www.gnu.org/licenses/>.

'''

import os
import sys

from . import __version__
from .constants import *
from .cluster_evolution import (
    initialize_cluster,
    restore_cluster,
    add_black_holes,
    _load_initial_bh_masses,
    form_neutron_stars,
    form_populations,
    form_pending_compact_objects,
    compute_cluster_properties,
    compute_timescales,
    form_binaries,
    evolve_interactions,
    compute_external_params,
    evolve_tdes,
    record_evolution,
    flush_records,
    compact_object_mass,
    update_cluster,
    print_status,
    write_output,
)
from .stellar_evolution import init_stellar_mass_sampler
from .compact_accretion import EOS_TABLES, accrete_gas
from .functions import set_escape_factor
from .ns_formation import DEFAULT_NS_FORMATION_TABLE
from .populations import read_populations_file, DEFAULT_N_MIX
from .restart import (RESTART_OPTIONS, checkpoint_path, save_checkpoint, load_checkpoint,
                      rows_before, write_restart_info)

class TeeStream:
    """Write to both a file and optionally to stdout.

    When print_to_screen is True, output goes to both the log file and the
    terminal. When False, output goes only to the log file.
    """

    def __init__(self, log_path, print_to_screen=True):
        self.log_file = open(log_path, 'w')
        self.stdout = sys.stdout
        self.print_to_screen = print_to_screen

    def write(self, text):
        self.log_file.write(text)
        self.log_file.flush()
        if self.print_to_screen:
            self.stdout.write(text)
            self.stdout.flush()

    def flush(self):
        self.log_file.flush()
        if self.print_to_screen:
            self.stdout.flush()

    def close(self):
        self.log_file.close()
        sys.stdout = self.stdout


def parse_args():
    """Parse command-line arguments and return a config dictionary.

    Returns:
        dict: Configuration dictionary with all simulation parameters including
            cluster properties (N, rh, Z, etc.), BH spin settings (s1g_max, SD),
            remnant mass prescription (RP), output file indicators, and I/O paths.
    """
    parser = argparse.ArgumentParser(description="Rapster input parameters")

    parser.add_argument('-N', '--number', type=float, metavar=' ', default=1000000, help='Initial number of stars')
    parser.add_argument('-r', '--half_mass_radius', type=float, metavar=' ', default=1.0, help='Initial half-mass radius [pc]')
    parser.add_argument('-mm', '--minimum_star_mass', type=float, metavar=' ', default=0.08, help='Smallest ZAMS mass [Msun]')
    parser.add_argument('-mM', '--maximum_star_mass', type=float, metavar=' ', default=150.0, help='Largest ZAMS mass [Msun]')
    parser.add_argument('-Z', '--metallicity', type=float, metavar=' ', default=0.001, help='Absolute metallicity')
    parser.add_argument('-z', '--cluster_formation_redshift', type=float, metavar=' ', default=3.0, help='Redshift of cluster formation')
    parser.add_argument('-n', '--central_stellar_density', type=float, metavar=' ', default=None, help='Initial central stellar number density [pc^-3]; default: Plummer value 0.53 N / r^3 (5.3e5 for N = 1e6, r = 1 pc)')
    parser.add_argument('-fb', '--binary_fraction', type=float, metavar=' ', default=0.1, help='Initial binary star fraction')
    parser.add_argument('-S', '--seed', type=int, metavar=' ', default=1234567890, help='Seed number')
    parser.add_argument('-dtm', '--minimum_time_step', type=float, metavar=' ', default=0.1, help='Minimum simulation time-step [Myr]')
    parser.add_argument('-dtM', '--maximum_time_step', type=float, metavar=' ', default=50.0, help='Maximum simulation time-step [Myr]')
    parser.add_argument('-tM', '--maximum_time', type=float, metavar=' ', default=14000.0, help='Maximum simulation time [Myr]')
    parser.add_argument('-fCO', '--max_compact_object_fraction', type=float, metavar=' ', default=fBH_max, help='The simulation stops when the mass in compact objects (single BHs and NSs and those in binaries, triples and BH-star pairs) exceeds this fraction of the stellar mass of the cluster')
    parser.add_argument('-wK', '--supernova_kick_parameter', type=float, metavar=' ', default=265.0, help='One-dimensional supernova kick parameter [km/s]')
    parser.add_argument('-K', '--natal_kick_prescription', type=int, metavar=' ', choices=[0, 1], default=1, help='Natal kick prescription (0 for fallback, 1 for momentum conservation kicks)')
    parser.add_argument('-R', '--galactocentric_radius', type=float, metavar=' ', default=8000.0, help='Initial galactocentric radius [pc]; 0 places the cluster at the galactic center (no galactic tides, no dynamical friction), e.g. for a nuclear star cluster')
    parser.add_argument('-GC', '--continue_at_galactic_center', type=int, metavar=' ', choices=[0, 1], default=0, help='When the cluster reaches the galactic center by dynamical friction: 0 to stop the simulation, 1 to continue it at the center without galactic tides')
    parser.add_argument('-vg', '--galactocentric_velocity', type=float, metavar=' ', default=220.0, help='Galactocentric circular velocity [km/s]')
    parser.add_argument('-s', '--spin_parameter', type=float, metavar=' ', default=0.0, help='Natal spin parameter of first generation (1g) BHs')
    parser.add_argument('-SD', '--spin_distribution', type=int, metavar=' ', choices=[0, 1, 2], default=0, help='Natal spin distribution model (0 for uniform, 1 for monochromatic, 2 for beta)')
    parser.add_argument('-P', '--print_information', type=int, metavar=' ', choices=[0, 1], default=1, help='Print runtime information (0 for no, 1 for yes)')
    parser.add_argument('-Mi', '--mergers_file_indicator', type=int, metavar=' ', choices=[0, 1], default=1, help='Export mergers file (0 for no, 1 for yes)')
    parser.add_argument('-MF', '--mergers_file_name', type=str, metavar=' ', default='mergers', help='Name of .txt output file with BBH merger source parameters')
    parser.add_argument('-Ei', '--evolution_file_indicator', type=int, metavar=' ', choices=[0, 1], default=1, help='Export evolution file (0 for no, 1 for yes)')
    parser.add_argument('-EF', '--evolution_file_name', type=str, metavar=' ', default='evolution', help='Name of .txt output file with time-dependent quantities')
    parser.add_argument('-Hi', '--hardening_file_indicator', type=int, metavar=' ', choices=[0, 1], default=1, help='Export hardening file (0 for no, 1 for yes)')
    parser.add_argument('-HF', '--hardening_file_name', type=str, metavar=' ', default='hardening', help='Name of .txt output file with BBH time evolution information')
    parser.add_argument('-BIi', '--blackholes_in_file_indicator', type=int, metavar=' ', choices=[0, 1], default=0, help='Use external BH file (0 for no, 1 for yes)')
    parser.add_argument('-BIF', '--blackholes_in_file_name', type=str, metavar=' ', default='input_BHs.npz', help='Name of .npz input file with initial BH masses')
    parser.add_argument('-BOi', '--blackholes_out_file_indicator', type=int, metavar=' ', choices=[0, 1], default=1, help='Export BH masses file (0 for no, 1 for yes)')
    parser.add_argument('-BOF', '--blackholes_out_file_name', type=str, metavar=' ', default='outputBHs', help='Name of .pkl file with the masses of all BHs in solar masses')
    parser.add_argument('-RP', '--remnant_mass_prescription', type=int, metavar=' ', choices=[0, 1, 2, 3], default=1, help='Remnant mass prescription (0 for SEVN delayed, 1 for Fryer+2012 delayed, 2 for SEVN rapid, 3 for Fryer+2012 rapid)')
    parser.add_argument('-NS', '--with_neutron_stars', type=int, metavar=' ', default=2, choices=[0, 1, 2, 3], help='include neutron stars (if =1 with monochromatic, =2 with bimodal or =3 with uniform mass distribution) else no (if =0)')
    parser.add_argument('-NSt', '--ns_formation_time_model', type=int, metavar=' ', default=0, choices=[0, 1, 2], help='NS formation-time model: 0 for all NSs at the BH formation time (default), 1 for times drawn from the COMPAS distribution p(t), 2 for times drawn from p(t | m_NS); the NS mass distribution is set by -NS')
    parser.add_argument('-WT', '--with_tdes', type=int, metavar=' ', choices=[0, 1], default=1, help='include tdes (if =1) else no (if =0): with 0, no TDEs of any kind (BH-star, BH-WD, BBH-star) and no white dwarfs')
    parser.add_argument('-Ti', '--tdes_file_indicator', type=int, metavar=' ', choices=[0, 1], default=1, help='Export tdes file (0 for no, 1 for yes)')
    parser.add_argument('-TF', '--tdes_file_name', type=str, metavar=' ', default='tdes', help='Name of .txt file containing tde parameters')
    parser.add_argument('-MBH', '--massive_black_hole_mass', type=float, metavar=' ', default=0, help='mass of the seed massive BH (if >0)')
    parser.add_argument('-sBH', '--massive_black_hole_spin', type=float, metavar=' ', default=0, help='spin of the seed massive BH (from 0 to 1)')
    parser.add_argument('-RF', '--results_folder_name', type=str, metavar=' ', default='Results', help='Name of the folder where output files will be exported')
    parser.add_argument('-BMD', '--bh_mass_distribution', type=int, metavar=' ', choices=[0, 1, 2, 3], default=0, help='Initial BH mass distribution (0 for Kroupa+collapse, 1 for uniform, 2 for Salpeter power law, 3 for log-uniform)')
    parser.add_argument('-mBH1gMin', '--min_1g_bh_mass', type=float, metavar=' ', default=3.0, help='Minimum 1g BH mass for uniform/Salpeter sampling [Msun]')
    parser.add_argument('-mBH1gMax', '--max_1g_bh_mass', type=float, metavar=' ', default=60.0, help='Maximum 1g BH mass for uniform/Salpeter sampling [Msun]')
    parser.add_argument('-RMP', '--random_mass_pairing_2body_3body', type=int, metavar=' ', choices=[0, 1], default=0, help='Use uniform random pairing for 3bb and 2-body capture instead of mass-weighted (0 for no, 1 for yes)')
    parser.add_argument('-plot', '--generate_plots', type=int, metavar=' ', choices=[0, 1], default=0, help='Generate diagnostic plots after simulation (0 for no, 1 for yes)')
    parser.add_argument('-analyze', '--analyze_results', type=int, metavar=' ', choices=[0, 1], default=0, help='Print analysis summary after simulation (0 for no, 1 for yes)')
    parser.add_argument('-fA', '--accreted_fraction', type=float, metavar=' ', default=0.5, help='Fraction of a disrupted star accreted by the compact object')
    parser.add_argument('-mb', '--mass_bias_power', type=float, metavar=' ', default=0.0, help='Mass bias power index p for drawing stars from IMF*m^p for TDEs')
    parser.add_argument('-EoS', '--equation_of_state', type=str, metavar=' ', default='APR', choices=sorted(EOS_TABLES), help='Neutron Star equation of state; one of ' + ', '.join(sorted(EOS_TABLES)))
    parser.add_argument('-RK', '--recoil_kick_model', type=int, metavar=' ', choices=[0, 1], default=0, help='GW recoil kick model (0 for Gerosa & Kesden 2016, 1 for gwModel_kick_prec_flow from Islam & Wadekar 2025)')
    parser.add_argument('-sfe', '--star_formation_efficiency', type=float, metavar=' ', default=1.0, help='Star formation efficiency epsilon in (0,1]; sets initial residual gas mass. Default 1.0 = no gas (recovers gas-free Rapster)')
    parser.add_argument('-tge', '--gas_expulsion_time', type=float, metavar=' ', default=10.0, help='Gas expulsion timescale [Myr]')
    parser.add_argument('-fEdd', '--eddington_ratio_cap', type=float, metavar=' ', default=1.0, help='Eddington ratio ceiling for gas accretion onto compact objects (1.0 = hard Eddington cap)')
    parser.add_argument('-cs', '--gas_sound_speed', type=float, metavar=' ', default=10.0, help='Gas sound speed [km/s]')
    parser.add_argument('-fesc', '--escape_velocity_factor', type=float, metavar=' ', default=2.0, help='Escape-velocity prefactor f in v_esc = f*sqrt(<v_star^2>+<v_BH^2>); 2 = rms escape speed (default)')
    parser.add_argument('-cxi', '--bh_temperature_prefactor', type=float, metavar=' ', default=1.0, help='Prefactor c_xi of the BH-to-star temperature ratio xi = c_xi q^(3/5) Q^(2/5) (lnL_BH/lnL)^(-2/5) (Breen & Heggie 2013); r_h,BH scales as 1/c_xi')
    parser.add_argument('-POP', '--populations_file', type=str, metavar=' ', default=None, help='Text file of stellar populations born after t = 0, one per line: t_birth [Myr] N [Z] [SFE] [x = r_h,k/r_h] [n_mix] (README, "Multiple stellar populations")')
    parser.add_argument('-POPi', '--populations_file_indicator', type=int, metavar=' ', choices=[0, 1], default=1, help='Export populations file when there are several populations (0 for no, 1 for yes)')
    parser.add_argument('-POPF', '--populations_file_name', type=str, metavar=' ', default='populations', help='Name of .txt output file with the properties of each stellar population')
    parser.add_argument('-CKi', '--checkpoint_indicator', type=int, metavar=' ', choices=[0, 1], default=1, help='Write checkpoint.pkl, to restart the simulation from, every -CKdt Myr and at the end (0 for no, 1 for yes)')
    parser.add_argument('-CKdt', '--checkpoint_interval', type=float, metavar=' ', default=500.0, help='Simulated time between checkpoints [Myr]')
    parser.add_argument('-RS', '--restart_from', type=str, metavar=' ', default=None, help='Restart from the checkpoint in this results folder (or this checkpoint file); the restarted run writes to a new folder (-RF). README, "Checkpoints and restarts"')
    parser.add_argument('-addBH', '--add_black_holes_file', type=str, metavar=' ', default=None, help='Restart only: .npz file of BHs to add (array mBH [Msun], optional array sBH of spins)')
    parser.add_argument('-addN', '--add_number', type=float, metavar=' ', default=None, help='Restart only: number of stars of a new stellar population born at the restart')
    parser.add_argument('-addZ', '--add_metallicity', type=float, metavar=' ', default=None, help='Restart only: metallicity of the new population (default: that of the cluster)')
    parser.add_argument('-addSFE', '--add_star_formation_efficiency', type=float, metavar=' ', default=1.0, help='Restart only: star formation efficiency of the new population')
    parser.add_argument('-addx', '--add_radius_ratio', type=float, metavar=' ', default=1.0, help='Restart only: half-mass radius of the new population relative to that of the cluster (1: well mixed)')
    parser.add_argument('-addnmix', '--add_mixing_time', type=float, metavar=' ', default=DEFAULT_N_MIX, help='Restart only: mixing time of the new population if segregated [half-mass relaxation times]')

    args = parser.parse_args()

    # options given on the command line (the others take the defaults; for a restart,
    # those of the run restarted from):
    for action in parser._actions:
        action.default = argparse.SUPPRESS
    given = set(vars(parser.parse_args()))

    # Initial central density: unless given with -n, the Plummer-profile value for the
    # chosen N and r_h, n = 3 (1.305)^3 N / (4 pi r_h^3) ~= 0.53 N / r_h^3, written as
    # 5.3e5 (N/1e6) (r_h/pc)^-3 so that the default N = 1e6, r_h = 1 pc gives exactly
    # 5.3e5 pc^-3 as before. (Previously the default stayed 5.3e5 whatever -N and -r.)
    if args.central_stellar_density is None:
        if args.half_mass_radius > 0 and np.isfinite(args.half_mass_radius) and np.isfinite(args.number):
            args.central_stellar_density = 5.3e5 * (args.number / 1e6) / args.half_mass_radius**3
        else:
            args.central_stellar_density = 5.3e5  # invalid -N/-r are reported below

    if args.blackholes_in_file_indicator and not os.path.isfile(args.blackholes_in_file_name):
        parser.error(f'initial black-hole input file not found: {args.blackholes_in_file_name}')

    float_parameters = (
        'number', 'half_mass_radius', 'minimum_star_mass', 'maximum_star_mass',
        'metallicity', 'cluster_formation_redshift', 'central_stellar_density',
        'binary_fraction', 'minimum_time_step', 'maximum_time_step', 'maximum_time',
        'supernova_kick_parameter', 'galactocentric_radius', 'galactocentric_velocity',
        'spin_parameter', 'massive_black_hole_mass', 'massive_black_hole_spin',
        'min_1g_bh_mass', 'max_1g_bh_mass', 'accreted_fraction', 'mass_bias_power',
        'star_formation_efficiency', 'gas_expulsion_time', 'eddington_ratio_cap',
        'gas_sound_speed', 'escape_velocity_factor', 'bh_temperature_prefactor',
        'max_compact_object_fraction',
    )
    for name in float_parameters:
        if not np.isfinite(getattr(args, name)):
            parser.error(f'{name} must be finite')

    positive_parameters = (
        'number', 'half_mass_radius', 'minimum_star_mass', 'central_stellar_density',
        'minimum_time_step', 'maximum_time_step', 'maximum_time',
        'galactocentric_velocity', 'gas_expulsion_time',
        'gas_sound_speed', 'max_compact_object_fraction',
    )
    for name in positive_parameters:
        if getattr(args, name) <= 0:
            parser.error(f'{name} must be positive')
    if args.galactocentric_radius < 0:
        parser.error('galactocentric_radius must be non-negative (0: at the galactic center)')

    if args.maximum_star_mass <= args.minimum_star_mass:
        parser.error('maximum_star_mass must exceed minimum_star_mass')
    if args.maximum_time_step < args.minimum_time_step:
        parser.error('maximum_time_step must be greater than or equal to minimum_time_step')
    if args.cluster_formation_redshift < 0 or args.metallicity < 0:
        parser.error('cluster_formation_redshift and metallicity must be non-negative')
    if not 0 <= args.binary_fraction <= 1:
        parser.error('binary_fraction must be in [0, 1]')
    if args.supernova_kick_parameter < 0 or args.eddington_ratio_cap < 0:
        parser.error('supernova_kick_parameter and eddington_ratio_cap must be non-negative')
    if not 0 <= args.accreted_fraction <= 1:
        parser.error('accreted_fraction must be in [0, 1]')
    if not 0 <= args.spin_parameter <= 1 or not 0 <= args.massive_black_hole_spin <= 1:
        parser.error('spin parameters must be in [0, 1]')
    if args.massive_black_hole_mass < 0:
        parser.error('massive_black_hole_mass must be non-negative')
    if not 0 <= args.seed <= 2**32 - 1:
        parser.error('seed must be in [0, 2**32 - 1]')

    if args.bh_mass_distribution != 0 and (
        args.min_1g_bh_mass <= 0 or args.max_1g_bh_mass <= args.min_1g_bh_mass
    ):
        parser.error('1g black-hole mass limits must be positive and increasing')

    if not (0 < args.star_formation_efficiency <= 1):
        parser.error('star_formation_efficiency (-sfe) must be in (0, 1]')

    if args.ns_formation_time_model > 0:
        if args.with_neutron_stars == 0:
            parser.error('ns_formation_time_model (-NSt) > 0 requires neutron stars (-NS 1, 2 or 3)')
        if not os.path.isfile(DEFAULT_NS_FORMATION_TABLE):
            parser.error(f'NS formation-time table not found: {DEFAULT_NS_FORMATION_TABLE}')

    if args.escape_velocity_factor <= 0:
        parser.error('escape_velocity_factor (-fesc) must be positive')

    if args.bh_temperature_prefactor <= 0:
        parser.error('bh_temperature_prefactor (-cxi) must be positive')

    if args.checkpoint_interval <= 0 or not np.isfinite(args.checkpoint_interval):
        parser.error('checkpoint_interval (-CKdt) must be positive')

    add_options = {'add_black_holes_file', 'add_number', 'add_metallicity', 'add_star_formation_efficiency',
                   'add_radius_ratio', 'add_mixing_time'}
    if args.restart_from is None:
        if given & add_options:
            parser.error('-addBH, -addN, -addZ, -addSFE, -addx and -addnmix are options of a restart (-RS); '
                         'use -POP for populations born later in a new run, and -BIi/-BIF for initial BHs')
    else:
        if not os.path.isfile(checkpoint_path(args.restart_from)):
            parser.error(f'checkpoint not found: {checkpoint_path(args.restart_from)}')
        option_names = {a.dest: a.option_strings[0] for a in parser._actions if a.option_strings}
        not_allowed = sorted(option_names[d] for d in given - RESTART_OPTIONS)
        if not_allowed:
            parser.error('a restart continues with the parameters of the run restarted from; it cannot set '
                         + ', '.join(not_allowed) + ' (README, "Checkpoints and restarts")')
        if os.path.realpath(os.path.join(os.getcwd(), args.results_folder_name)) == \
                os.path.realpath(os.path.dirname(os.path.abspath(checkpoint_path(args.restart_from)))):
            parser.error('the restarted run must write to a new results folder (-RF)')
        if (given & add_options) - {'add_black_holes_file', 'add_number'} and 'add_number' not in given:
            parser.error('-addZ, -addSFE, -addx and -addnmix describe a population added with -addN')
        if args.add_number is not None and not (np.isfinite(args.add_number) and args.add_number >= 1):
            parser.error('add_number (-addN) must be at least 1')
        if args.add_metallicity is not None and not (np.isfinite(args.add_metallicity) and args.add_metallicity >= 0):
            parser.error('add_metallicity (-addZ) must be non-negative')
        if not 0 < args.add_star_formation_efficiency <= 1:
            parser.error('add_star_formation_efficiency (-addSFE) must be in (0, 1]')
        if not (np.isfinite(args.add_radius_ratio) and args.add_radius_ratio > 0):
            parser.error('add_radius_ratio (-addx) must be positive')
        if not (np.isfinite(args.add_mixing_time) and args.add_mixing_time > 0):
            parser.error('add_mixing_time (-addnmix) must be positive')
        if args.add_black_holes_file is not None and not os.path.isfile(args.add_black_holes_file):
            parser.error(f'BH file not found: {args.add_black_holes_file}')

    populations = []
    if args.populations_file is not None:
        if not os.path.isfile(args.populations_file):
            parser.error(f'populations file not found: {args.populations_file}')
        try:
            populations = read_populations_file(args.populations_file, args.metallicity)
        except ValueError as err:
            parser.error(str(err))

    config = {
        'N': args.number,
        'rh': args.half_mass_radius,
        'm_min': args.minimum_star_mass,
        'm_max': args.maximum_star_mass,
        'Z': args.metallicity,
        'seed': args.seed,
        'dt_min': args.minimum_time_step,
        'dt_max': args.maximum_time_step,
        't_max': args.maximum_time,
        'wSN_kick': args.supernova_kick_parameter,
        'R_gal': args.galactocentric_radius,
        'v_gal': args.galactocentric_velocity,
        'zCl_form': args.cluster_formation_redshift,
        'n_star': args.central_stellar_density,
        'NKP': args.natal_kick_prescription,
        'mergers_file': args.mergers_file_name,
        'evolution_file': args.evolution_file_name,
        'hardening_file': args.hardening_file_name,
        'fb': args.binary_fraction,
        'print_info': args.print_information,
        'Bi': args.blackholes_in_file_indicator,
        'input_BH_file': args.blackholes_in_file_name,
        's1g_max': args.spin_parameter,
        'SD': args.spin_distribution,
        'Mi': args.mergers_file_indicator,
        'Ei': args.evolution_file_indicator,
        'Hi': args.hardening_file_indicator,
        'BOi': args.blackholes_out_file_indicator,
        'BOF': args.blackholes_out_file_name,
        'RP': args.remnant_mass_prescription,
        'with_NSs': args.with_neutron_stars,
        'NSt': args.ns_formation_time_model,
        'with_tdes': args.with_tdes,
        'Ti': args.tdes_file_indicator,
        'tdes_file': args.tdes_file_name,
        'M_BH0': args.massive_black_hole_mass,
        's_BH0': args.massive_black_hole_spin,
        'results_folder_name': args.results_folder_name,
        'BMD': args.bh_mass_distribution,
        'min_1g_bh_mass': args.min_1g_bh_mass,
        'max_1g_bh_mass': args.max_1g_bh_mass,
        'random_pairing': bool(args.random_mass_pairing_2body_3body),
        'generate_plots': bool(args.generate_plots),
        'analyze_results': bool(args.analyze_results),
        'f_accreted': args.accreted_fraction,
        'mass_bias_power': args.mass_bias_power,
        'EoS': args.equation_of_state,
        'recoil_kick_model': args.recoil_kick_model,
        'SFE': args.star_formation_efficiency,
        't_ge': args.gas_expulsion_time,
        'fCO_max': args.max_compact_object_fraction,
        'continue_at_GC': args.continue_at_galactic_center,
        'f_Edd': args.eddington_ratio_cap,
        'c_s': args.gas_sound_speed,
        'f_esc': args.escape_velocity_factor,
        'c_xi': args.bh_temperature_prefactor,
        'populations': populations,
        'populations_path': args.populations_file,
        'POPi': args.populations_file_indicator,
        'populations_file': args.populations_file_name,
        'CKi': args.checkpoint_indicator,
        'CKdt': args.checkpoint_interval,
        'restart_from': args.restart_from,
        'add_BH_file': args.add_black_holes_file,
        'add_population': (None if args.add_number is None else dict(
            N=args.add_number, Z=args.add_metallicity, SFE=args.add_star_formation_efficiency,
            x=args.add_radius_ratio, n_mix=args.add_mixing_time)),
        'given': sorted(given),
    }

    return config


# settings of a run (as opposed to parameters of the cluster) that a restart takes from its
# own command line, given or not:
_RUN_SETTINGS = ('results_folder_name', 'print_info', 'generate_plots', 'analyze_results',
                 'CKi', 'CKdt', 'restart_from', 'add_BH_file', 'add_population', 'given')

# parameters that a restart takes from its command line only when given (option -> key):
_RESTART_OVERRIDES = {'maximum_time': 't_max', 'seed': 'seed',
                      'max_compact_object_fraction': 'fCO_max',
                      'continue_at_galactic_center': 'continue_at_GC'}


def restart_config(saved, current):
    """Configuration of a restarted run: that of the run restarted from (saved), with the
    run settings and the parameters given on the restart's command line (current)."""
    config = dict(saved)
    for key in _RUN_SETTINGS:
        config[key] = current[key]
    for option, key in _RESTART_OVERRIDES.items():
        if option in current['given']:
            config[key] = current[key]
    config['rh_restart'] = current['rh'] if 'half_mass_radius' in current['given'] else None
    config['populations_path'] = current['populations_path']
    return config


def _next_checkpoint_time(t, interval):
    """First multiple of the checkpoint interval after time t."""
    return (np.floor(t / interval) + 1) * interval


def _start_restart(checkpoint, config):
    """State of a restarted run, with the changes asked for on its command line."""
    t = float(checkpoint['state']['t'])
    new_pops = []
    if config['populations_path'] is not None:
        try:
            new_pops = read_populations_file(config['populations_path'], checkpoint['config']['Z'])
        except ValueError as err:
            raise SystemExit(f'rapster: error: {err}')
        early = [p['t_birth'] for p in new_pops if not p['t_birth'] > t]
        if early:
            raise SystemExit('rapster: error: populations of a restart (-POP) must be born after the restart '
                             '(t = {:.6g} Myr); got t_birth = {}'.format(t, ', '.join('{:.6g}'.format(x) for x in early)))
    add = config['add_population']
    if add is not None:
        new_pops.insert(0, dict(t_birth=t, N=add['N'], Z=checkpoint['config']['Z'] if add['Z'] is None else add['Z'],
                                SFE=add['SFE'], x=add['x'], n_mix=add['n_mix']))

    state = restore_cluster(checkpoint, config, with_populations=bool(new_pops))
    if checkpoint['records'].get('populations') is None and state.get('populations') is not None:
        # populations.txt starts at this restart: give it the rows of the evolution rows of
        # this run from before the restart (one population then, the whole cluster), as an
        # uninterrupted run with -POP would
        evolution, record = state['evolution'], state['populations']
        rows = np.frombuffer(evolution.data, dtype=np.float64).reshape(-1, evolution.ncols)
        first = 1 if evolution.removed == 0 else 0          # (the placeholder row)
        for i in range(first, rows.shape[0]):
            for row in state['pops'].rows(rows[i, 1], state['m_avg0'], rows[i, 5]):
                record.append(row)
            state['population_marks'].append((evolution.removed + i + 1, len(record)))
        del state['population_marks'][:-4]
    changes = {}

    # a new seed of the random numbers (otherwise the generator continues where it was):
    if 'seed' in config['given']:
        np.random.seed(config['seed'])
        state['seed'] = config['seed']
        changes['seed'] = config['seed']

    # a new half-mass radius (the central density follows from it, as during the run):
    if config['rh_restart'] is not None:
        changes['r_h'] = [float(state['rh']), float(config['rh_restart'])]
        state['rh'] = config['rh_restart']

    # BHs added:
    if config['add_BH_file'] is not None:
        archive = np.load(config['add_BH_file'], allow_pickle=False)
        try:
            key = 'mBH' if 'mBH' in archive.files else 'mBH_ini'
            if key not in archive.files:
                raise SystemExit('rapster: error: the BH file of -addBH must contain an array mBH of masses [Msun]')
            mBH = np.asarray(archive[key], dtype=float).ravel()
            sBH = np.asarray(archive['sBH'], dtype=float).ravel() if 'sBH' in archive.files else None
        finally:
            archive.close()
        if mBH.size == 0 or not np.all(np.isfinite(mBH)) or np.any(mBH <= 0):
            raise SystemExit('rapster: error: BH masses of -addBH must be finite and positive')
        if sBH is not None and (sBH.shape != mBH.shape or not np.all((sBH >= 0) & (sBH <= 1))):
            raise SystemExit('rapster: error: BH spins of -addBH must be in [0, 1], one per mass')
        add_black_holes(state, config, mBH, sBH)
        changes['BHs_added'] = int(mBH.size)

    # stellar populations to be born (the one of -addN at the restart time):
    if new_pops:
        pops = state['pops']
        pops.scheduled = sorted(pops.scheduled + new_pops, key=lambda p: p['t_birth'])
        config['populations'] = list(config.get('populations', ())) + new_pops
        changes['populations_added'] = new_pops

    return state, changes


def main():
    """Main entry point for running a Rapster cluster simulation."""

    config = parse_args()

    checkpoint = None
    if config['restart_from'] is not None:
        checkpoint = load_checkpoint(config['restart_from'])
        config = restart_config(checkpoint['config'], config)
        if not checkpoint['state']['t'] < config['t_max']:
            raise SystemExit('rapster: error: the checkpoint is at t = {:.6g} Myr, not before the maximum time '
                             't_max = {:.6g} Myr: give a later -tM'.format(checkpoint['state']['t'], config['t_max']))

    # set up output directory and logging:
    results_dir = os.path.join(os.getcwd(), config['results_folder_name'])
    os.makedirs(results_dir, exist_ok=True)
    log_path = os.path.join(results_dir, 'log.txt')
    tee = TeeStream(log_path, print_to_screen=bool(config['print_info']))
    sys.stdout = tee

    # start initialization clock:
    initialization_time_initial = time.time()
    print('RAPSTER VERSION:', __version__)
    init_stellar_mass_sampler(config['mass_bias_power'])
    set_escape_factor(config['f_esc'])
    if checkpoint is None:
        print('INITIALIZING...')
        print('Initial central stellar density: {:.6g} pc^-3'.format(config['n_star']))
        state = initialize_cluster(config)
    else:
        source = os.path.abspath(checkpoint_path(config['restart_from']))
        print('RESTARTING FROM {} (RAPSTER VERSION {}) AT t = {:.6g} Myr'.format(
            source, checkpoint['version'], checkpoint['state']['t']))
        state, changes = _start_restart(checkpoint, config)
        for key, value in changes.items():
            print('RESTART CHANGE:', key, value)
        restart_info = {
            'rapster_version': __version__, 'restarted_from': source,
            'checkpoint_version': checkpoint['version'], 't_restart': checkpoint['state']['t'],
            'rows_before': rows_before(checkpoint['records']),
            'snapshots_before': checkpoint['n_snapshots'],
            'file_names': {'mergers': config['mergers_file'], 'evolution': config['evolution_file'],
                           'hardening': config['hardening_file'], 'tdes': config['tdes_file'],
                           'populations': config.get('populations_file', 'populations'),
                           'outputBHs': config['BOF']},
            'changes': changes,
        }
        write_restart_info(results_dir, restart_info)

    print('END OF INITIALIZATION. RUNTIME:', "{:.3g}".format(np.abs(time.time() - initialization_time_initial)), 's')
    print('\n')

    print('SIMULATING...')
    print('\n')

    # start simulation clock:
    simulation_time_initial = time.time()

    # Simulation:
    # (R_gal = 0: at the galactic center, see -R and -GC; the simulation also stops when the
    #  compact objects exceed a fraction fCO_max of the stellar mass, see -fCO)
    # checkpoints (restart.py), every CKdt Myr of simulated time:
    next_checkpoint = _next_checkpoint_time(state['t'], config['CKdt'])

    while state['t']<config['t_max'] and state['R_gal']>=0 and state['Mcl']>0 and compact_object_mass(state)<config['fCO_max']*state['Mcl']:

        # checkpoint at the start of the step:
        if config['CKi'] == 1 and state['t'] >= next_checkpoint:
            save_checkpoint(state, config, __version__)
            next_checkpoint = _next_checkpoint_time(state['t'], config['CKdt'])

        # start local clock:
        local_time_initial = time.time()

        # stellar populations born now (-POP), and BHs and NSs of later populations that
        # form now:
        form_populations(state, config)
        form_pending_compact_objects(state, config)

        # form NSs whose formation time has been reached (NSt > 0):
        form_neutron_stars(state, config)

        # compute cluster and BH subsystem properties:
        keep_going = compute_cluster_properties(state, config)
        if not keep_going:
            break

        # compute dynamical timescales and adaptive timestep:
        compute_timescales(state, config)

        # binary BH formation (3bb, 2-body capture, exchanges):
        form_binaries(state, config)

        # BBH evolution, pair-pair interactions, triple evolution:
        evolve_interactions(state, config)

        # external parameters (Jacobi radius, escape rate, dynamical friction):
        compute_external_params(state, config)

        # tidal disruption events:
        evolve_tdes(state, config)

        # append evolution record:
        record_evolution(state)

        # write finished rows of the evolution, hardening and TDE records to their files:
        flush_records(state)

        # gas accretion onto compact objects (drains the gas reservoir):
        accrete_gas(state, config)

        # cluster structural evolution (mass loss, expansion, time update):
        keep_going = update_cluster(state, config)
        if not keep_going:
            break

        # print runtime information:
        print_status(state, config, local_time_initial, simulation_time_initial)

    else:
        # the loop condition ended the run (the other stops print their reason where
        # they occur): say which condition
        if not state['t'] < config['t_max']:
            print('MAXIMUM SIMULATION TIME REACHED (t_max = {:.6g} Myr)'.format(config['t_max']))
        elif not compact_object_mass(state) < config['fCO_max']*state['Mcl']:
            print('COMPACT OBJECTS REACHED A FRACTION fCO_max = {:.3g} OF THE CLUSTER MASS'.format(config['fCO_max']))
        elif not state['R_gal'] >= 0:
            print('CLUSTER REACHED GALAXY CENTER')
        elif not state['Mcl'] > 0:
            print('CLUSTER DISSOLVED')
        # checkpoint at the end, to continue the run later (e.g. to a larger t_max):
        if config['CKi'] == 1:
            save_checkpoint(state, config, __version__)

    write_output(state, config)

    if checkpoint is not None:
        restart_info['merger_suffix'] = state['merger_suffix']
        restart_info['t_end'] = state['t']
        write_restart_info(results_dir, restart_info)

    if config['NSt']>0:
        print('NEUTRON STARS: formed', state['N_NS_formed'], '| retained', state['N_NS_retained'],
              '| still to form at end of run', state['tNS_pending'].size + int(state['nCO_pending'].sum()))
    if state['pops'].size > 1 or config.get('populations'):
        print('STELLAR POPULATIONS: born', state['pops'].size, '| BHs and NSs of later populations formed',
              state['N_pop_CO_formed'], '| retained', state['N_pop_CO_retained'],
              '| still to form at end of run', state['tCO_pending'].size)

    print('END OF SIMULATION. RUNTIME:', "{:.3g}".format(np.abs(time.time() - simulation_time_initial)), 's')
    print('\n')

    # print analysis summary if requested (imported here: only needed on request):
    if config['analyze_results']:
        from .analyze_cluster import analyze_cluster
        analyze_cluster(results_dir)

    # generate diagnostic plots if requested (matplotlib is imported only then):
    if config['generate_plots']:
        from .plot_cluster import generate_all_plots
        generate_all_plots(results_dir)

    # close log file and restore stdout:
    tee.close()


if __name__ == "__main__":
    main()

# End of file.
