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
    form_neutron_stars,
    compute_cluster_properties,
    compute_timescales,
    form_binaries,
    evolve_interactions,
    compute_external_params,
    evolve_tdes,
    record_evolution,
    update_cluster,
    print_status,
    write_output,
)
from .stellar_evolution import init_stellar_mass_sampler
from .compact_accretion import EOS_TABLES, accrete_gas
from .functions import set_escape_factor
from .ns_formation import DEFAULT_NS_FORMATION_TABLE

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
    parser.add_argument('-wK', '--supernova_kick_parameter', type=float, metavar=' ', default=265.0, help='One-dimensional supernova kick parameter [km/s]')
    parser.add_argument('-K', '--natal_kick_prescription', type=int, metavar=' ', choices=[0, 1], default=1, help='Natal kick prescription (0 for fallback, 1 for momentum conservation kicks)')
    parser.add_argument('-R', '--galactocentric_radius', type=float, metavar=' ', default=8000.0, help='Initial galactocentric radius [pc]')
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
    parser.add_argument('-WT', '--with_tdes', type=int, metavar=' ', choices=[0, 1], default=1, help='include tdes (if =1) else no (if =0)')
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
    parser.add_argument('-fge', '--gas_expulsion_tcross', type=float, metavar=' ', default=5.0, help='Gas expulsion timescale in units of the initial crossing time')
    parser.add_argument('-fEdd', '--eddington_ratio_cap', type=float, metavar=' ', default=1.0, help='Eddington ratio ceiling for gas accretion onto compact objects (1.0 = hard Eddington cap)')
    parser.add_argument('-cs', '--gas_sound_speed', type=float, metavar=' ', default=10.0, help='Gas sound speed [km/s]')
    parser.add_argument('-fesc', '--escape_velocity_factor', type=float, metavar=' ', default=2.0, help='Escape-velocity prefactor f in v_esc = f*sqrt(<v_star^2>+<v_BH^2>); 2 = rms escape speed (default)')
    parser.add_argument('-cxi', '--bh_temperature_prefactor', type=float, metavar=' ', default=1.0, help='Prefactor c_xi of the BH-to-star temperature ratio xi = c_xi q^(3/5) Q^(2/5) (lnL_BH/lnL)^(-2/5) (Breen & Heggie 2013); r_h,BH scales as 1/c_xi')

    args = parser.parse_args()

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
        'star_formation_efficiency', 'gas_expulsion_tcross', 'eddington_ratio_cap',
        'gas_sound_speed', 'escape_velocity_factor', 'bh_temperature_prefactor',
    )
    for name in float_parameters:
        if not np.isfinite(getattr(args, name)):
            parser.error(f'{name} must be finite')

    positive_parameters = (
        'number', 'half_mass_radius', 'minimum_star_mass', 'central_stellar_density',
        'minimum_time_step', 'maximum_time_step', 'maximum_time',
        'galactocentric_radius', 'galactocentric_velocity', 'gas_expulsion_tcross',
        'gas_sound_speed',
    )
    for name in positive_parameters:
        if getattr(args, name) <= 0:
            parser.error(f'{name} must be positive')

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
        'f_ge': args.gas_expulsion_tcross,
        'f_Edd': args.eddington_ratio_cap,
        'c_s': args.gas_sound_speed,
        'f_esc': args.escape_velocity_factor,
        'c_xi': args.bh_temperature_prefactor,
    }

    return config


def main():
    """Main entry point for running a Rapster cluster simulation."""

    config = parse_args()

    # set up output directory and logging:
    results_dir = os.path.join(os.getcwd(), config['results_folder_name'])
    os.makedirs(results_dir, exist_ok=True)
    log_path = os.path.join(results_dir, 'log.txt')
    tee = TeeStream(log_path, print_to_screen=bool(config['print_info']))
    sys.stdout = tee

    # start initialization clock:
    initialization_time_initial = time.time()
    print('RAPSTER VERSION:', __version__)
    print('INITIALIZING...')
    print('Initial central stellar density: {:.6g} pc^-3'.format(config['n_star']))

    init_stellar_mass_sampler(config['mass_bias_power'])
    set_escape_factor(config['f_esc'])
    state = initialize_cluster(config)

    print('END OF INITIALIZATION. RUNTIME:', "{:.3g}".format(np.abs(time.time() - initialization_time_initial)), 's')
    print('\n')

    print('SIMULATING...')
    print('\n')

    # start simulation clock:
    simulation_time_initial = time.time()

    # Simulation:
    while state['t']<config['t_max'] and state['R_gal']>0 and state['Mcl']>0 and state['BH'].m.sum()<fBH_max*state['Mcl']:

        # start local clock:
        local_time_initial = time.time()

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

        # gas accretion onto compact objects (drains the gas reservoir):
        accrete_gas(state, config)

        # cluster structural evolution (mass loss, expansion, time update):
        keep_going = update_cluster(state, config)
        if not keep_going:
            break

        # print runtime information:
        print_status(state, config, local_time_initial, simulation_time_initial)

    write_output(state, config)

    if config['NSt']>0:
        print('NEUTRON STARS: formed', state['N_NS_formed'], '| retained', state['N_NS_retained'],
              '| still to form at end of run', state['tNS_pending'].size)

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
