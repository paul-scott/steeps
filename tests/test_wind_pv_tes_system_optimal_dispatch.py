#! /bin/env python

from __future__ import division
import unittest
import argparse  # Add argparse for command line argument parsing

from solartherm import simulation
import DyMat

import os
import numpy as np
import time
from math import pi, ceil
import matplotlib.pyplot as plt

class TestScheduler(unittest.TestCase):
    def setUp(self):
        fn = f'{st_folder}/examples/WindPVSimpleSystemOptimalDispatch.mo'
        sim = simulation.Simulator(fn)
        sim.compile_model()
        sim.compile_sim(args=['-s'])
        tinit = time.time()
        sim.simulate(start=0, stop='30d', step='1h', solver='dassl', nls='homotopy', tolerance='1e-06', args=['-noEventEmit'])
        print('Simulation time: %4.2f s' % (time.time() - tinit))

    def test_sched(self):
        df = DyMat.DyMatFile('WindPVSimpleSystemOptimalDispatch_res.mat')
        Q_flow_dis = df.data('Q_flow_dis') / 1.e6
        Q_flow_chg = df.data('renewable_input.electricity') / 1.e6
        optimalDispatch = df.data('optimalDispatch')
        P_elec_net = df.data('renewable_input.P_elec_net') / 1.e6
        blk_state = df.data('blk_state')
        t_startup_next_fun = df.data('t_startup_next') / 3600.
        runtime = df.data('runtime')
        pre_blk_state = df.data('pre_blk_state')
        E_max = df.data('E_max')[0]
        E = df.data('E') / E_max * 100
        times = df.abscissa('Q_flow_dis')[0] / 3600. / 24.
        t_storage = df.data('t_storage')[0]
        renewable_multiple = df.data('renewable_multiple')[0]
        heater_multiple = df.data('heater_multiple')[0]
        pv_fraction = df.data('pv_fraction')[0]

        fig, axes = plt.subplots(3, 1, figsize=(30 / 2.54, 22.5 / 2.54), sharex=False)
        line1, = axes[0].plot(times, Q_flow_dis, ls='--', marker='', color='tab:red', markevery=5, label='Process Heat Input [MWt]', zorder=2.5)
        line2, = axes[0].plot(times, P_elec_net, ls='-', marker='', color='tab:blue', label='Renewable Input [MWe]', zorder=2, alpha=0.5)
        line3, = axes[0].plot(times, optimalDispatch, ls='-', marker='', color='tab:orange', label='Optimal Dispatch [MWt]', zorder=2.1)
        axes[0].set_xlabel('time (d)')
        axes[0].set_ylabel('Power [MW]')

        max_Q_flow_chg = np.amax(Q_flow_chg)
        upper_limit = ceil(max_Q_flow_chg / 10) * 10 + 10
        axes[0].set_ylim([0, upper_limit])

        axe2 = axes[0].twinx()
        line4, = axe2.plot(times, E, ls='--', color='tab:green', label='Storage Level [%]', alpha=0.7)
        axe2.set_ylabel('Storage Level [%]')
        axe2.set_ylim([0, 100])

        # Combine legends from both axes
        lines = [line1, line2, line3, line4]
        labels = [line.get_label() for line in lines]

        # Place legend outside the plot area
        plt.legend(lines, labels, loc='lower left', bbox_to_anchor=(0, 1.02, 1, 0.2), ncol=4)

        axes[1].plot(times, blk_state, ls='-', marker='', color='tab:red', markevery=2, label='Controller state', zorder=2.5)
        axes[1].plot(times, pre_blk_state, ls='-', marker='*', color='tab:blue', markevery=2, label='Pre Controller state', zorder=2.3)
        axes[1].set_xlabel('time (d)')
        axes[1].set_yticks([1, 2, 3, 4])
        axes[1].set_yticklabels(['Off', 'Ramp-up', 'On', 'Ramp-down'])
        axes[1].legend(loc='upper left')

        if args.plot_runtime:  # Use the parsed argument
            bar_width = (times[1] - times[0]) * 0.5  # Adjust the width to prevent overlap
            axes[2].bar(times, runtime, width=bar_width, edgecolor='k', linewidth=0.5)
            axes[2].set_xlabel('time (d)')
            axes[2].set_ylabel('Runtime (s)')
        else:
            axes[2].plot(times, t_startup_next_fun, ls='-', marker='', color='tab:orange', markevery=2, label='t_startup_next_fun', zorder=2.0)
            axes[2].set_xlabel('time (d)')
            axes[2].set_ylabel('Next startup time (h)')

        days = int(times[-1])
        plt.tight_layout()
        filename = f'windpv_annual_TES{t_storage:.0f}_RM{renewable_multiple:.1f}_HM{heater_multiple:.1f}_PV{pv_fraction:.1f}_{days}d_res'
        plt.savefig(f'{st_folder}/examples/{filename}.png', dpi=300)
        plt.show()

        # Deleting simulation files
        os.system(f'mv WindPVSimpleSystemOptimalDispatch_res.mat {st_folder}/examples/{filename}.mat')
        os.system('rm WindPVSimpleSystemOptimalDispatch*')

        csv = np.c_[times * 24. * 3600., Q_flow_dis, Q_flow_chg, optimalDispatch, P_elec_net, blk_state, E]
        np.savetxt("WindPVSimpleSystemOptimalDispatch_res.csv", csv, delimiter=",", header="times,Q_flow_dis,Q_flow_chg,optimalDispatch,P_elec_net,blk_state,E")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run the TestScheduler and plot results.')
    parser.add_argument('--plot_runtime', action='store_true', help='If set, plots runtime instead of next startup time.')
    args = parser.parse_args()

    script_dir = os.path.dirname(__file__)  # Get the directory where the script is located
    st_folder = os.path.abspath(os.path.join(script_dir, os.pardir))  # Set st_folder as the parent directory
    commands = f'cd {st_folder} && scons && scons install && cd tests'
    os.system(commands)
    unittest.main(argv=['first-arg-is-ignored'])  # unittest will parse its own arguments
