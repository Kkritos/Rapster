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

from .constants import *
from .cluster_evolution import (
    initialize_cluster,
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
    parser.add_argument('-n', '--central_stellar_density', type=float, metavar=' ', default=5.3e5, help='Initial central stellar number density [pc^-3]')
    parser.add_argument('-fb', '--binary_fraction', type=float, metavar=' ', default=0.1, help='Initial binary star fraction')
    parser.add_argument('-S', '--seed', type=int, metavar=' ', default=1234567890, help='Seed number')
    parser.add_argument('-dtm', '--minimum_time_step', type=float, metavar=' ', default=0.1, help='Minimum simulation time-step [Myr]')
    parser.add_argument('-dtM', '--maximum_time_step', type=float, metavar=' ', default=50.0, help='Maximum simulation time-step [Myr]')
    parser.add_argument('-tM', '--maximum_time', type=float, metavar=' ', default=14000.0, help='Maximum simulation time [Myr]')
    parser.add_argument('-wK', '--supernova_kick_parameter', type=float, metavar=' ', default=265.0, help='One-dimensional supernova kick parameter [km/s]')
    parser.add_argument('-K', '--natal_kick_prescription', type=int, metavar=' ', default=1, help='Natal kick prescription (0 for fallback, 1 for momentum conservation kicks)')
    parser.add_argument('-R', '--galactocentric_radius', type=float, metavar=' ', default=8000.0, help='Initial galactocentric radius [pc]')
    parser.add_argument('-vg', '--galactocentric_velocity', type=float, metavar=' ', default=220.0, help='Galactocentric circular velocity [km/s]')
    parser.add_argument('-s', '--spin_parameter', type=float, metavar=' ', default=0.0, help='Natal spin parameter of first generation (1g) BHs')
    parser.add_argument('-SD', '--spin_distribution', type=int, metavar=' ', default=0, help='Natal spin distribution model (0 for uniform, 1 for monochromatic, 2 for beta)')
    parser.add_argument('-P', '--print_information', type=int, metavar=' ', default=1, help='Print runtime information (0 for no, 1 for yes)')
    parser.add_argument('-Mi', '--mergers_file_indicator', type=int, metavar=' ', default=1, help='Export mergers file (0 for no, 1 for yes)')
    parser.add_argument('-MF', '--mergers_file_name', type=str, metavar=' ', default='mergers', help='Name of .txt output file with BBH merger source parameters')
    parser.add_argument('-Ei', '--evolution_file_indicator', type=int, metavar=' ', default=1, help='Export evolution file (0 for no, 1 for yes)')
    parser.add_argument('-EF', '--evolution_file_name', type=str, metavar=' ', default='evolution', help='Name of .txt output file with time-dependent quantities')
    parser.add_argument('-Hi', '--hardening_file_indicator', type=int, metavar=' ', default=1, help='Export hardening file (0 for no, 1 for yes)')
    parser.add_argument('-HF', '--hardening_file_name', type=str, metavar=' ', default='hardening', help='Name of .txt output file with BBH time evolution information')
    parser.add_argument('-BIi', '--blackholes_in_file_indicator', type=int, metavar=' ', default=0, help='Use external BH file (0 for no, 1 for yes)')
    parser.add_argument('-BIF', '--blackholes_in_file_name', type=str, metavar=' ', default='inputBHs.npz', help='Name of .npz input file with initial BH masses')
    parser.add_argument('-BOi', '--blackholes_out_file_indicator', type=int, metavar=' ', default=1, help='Export BH masses file (0 for no, 1 for yes)')
    parser.add_argument('-BOF', '--blackholes_out_file_name', type=str, metavar=' ', default='outputBHs', help='Name of .pkl file with the masses of all BHs in solar masses')
    parser.add_argument('-RP', '--remnant_mass_prescription', type=int, metavar=' ', default=1, help='Remnant mass prescription (0 for SEVN delayed, 1 for Fryer+2012 delayed, 2 for SEVN rapid, 3 for Fryer+2012 rapid)')
    parser.add_argument('-NS', '--with_neutron_stars', type=int, metavar=' ', default=1, help='include neutron stars (if =1) else no (if =0)')
    parser.add_argument('-WT', '--with_tdes', type=int, metavar=' ', default=1, help='include tdes (if =1) else no (if =0)')
    parser.add_argument('-Ti', '--tdes_file_indicator', type=int, metavar=' ', default=1, help='Export tdes file (0 for no, 1 for yes)')
    parser.add_argument('-TF', '--tdes_file_name', type=str, metavar=' ', default='tdes', help='Name of .txt file containing tde parameters')
    parser.add_argument('-MBH', '--massive_black_hole_mass', type=float, metavar=' ', default=0, help='mass of the seed massive BH (if >0)')
    parser.add_argument('-sBH', '--massive_black_hole_spin', type=float, metavar=' ', default=0, help='spin of the seed massive BH (from 0 to 1)')
    parser.add_argument('-RF', '--results_folder_name', type=str, metavar=' ', default='Results', help='Name of the folder where output files will be exported')
    parser.add_argument('-RMP', '--random_mass_pairing_2body_3body', type=int, metavar=' ', default=0, help='Use uniform random pairing for 3bb and 2-body capture instead of mass-weighted (0 for no, 1 for yes)')

    args = parser.parse_args()

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
        'with_tdes': args.with_tdes,
        'Ti': args.tdes_file_indicator,
        'tdes_file': args.tdes_file_name,
        'M_BH0': args.massive_black_hole_mass,
        's_BH0': args.massive_black_hole_spin,
        'results_folder_name': args.results_folder_name,
        'random_pairing': bool(args.random_mass_pairing_2body_3body),
    }

    return config


def main():
    """Main entry point for running a Rapster cluster simulation."""

    config = parse_args()

    # start initialization clock:
    initialization_time_initial = time.time()
    print('INITIALIZING...')

    state = initialize_cluster(config)

    print('END OF INITIALIZATION. RUNTIME:', "{:.3g}".format(np.abs(time.time() - initialization_time_initial)), 's')
    print('\n')

    print('SIMULATING...')
    print('\n')

    # start simulation clock:
    simulation_time_initial = time.time()

    # Simulation:
    while state['t']<config['t_max'] and state['R_gal']>0 and state['Mcl']>0 and state['mBH'].sum()<fBH_max*state['Mcl']:

        # start local clock:
        local_time_initial = time.time()

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

        # cluster structural evolution (mass loss, expansion, time update):
        keep_going = update_cluster(state, config)
        if not keep_going:
            break

        # print runtime information:
        print_status(state, config, local_time_initial, simulation_time_initial)

    write_output(state, config)

    print('END OF SIMULATION. RUNTIME:', "{:.3g}".format(np.abs(time.time() - simulation_time_initial)), 's')
    print('\n')


if __name__ == "__main__":
    main()

# End of file.
