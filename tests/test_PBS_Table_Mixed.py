#! /bin/env python

from __future__ import division
import unittest

import cleantest
from solartherm import simulation
from solartherm import postproc

from math import pi
import os

class TestReference(unittest.TestCase):
	def setUp(self):
		fn = '../examples/PBS_Surround_SCO2NREL_Table_Mixed.mo'
		sim = simulation.Simulator(fn)
		sim.compile_model()
		sim.compile_sim(args=['-s'])
		sim.simulate(start=0, stop='1y', step='300s', tolerance='1e-03', maxStep='60s', initStep='60s', solver='dassl', nls='mixed', args=['-noEventEmit'])
		self.res = postproc.SimResultElec(sim.res_fn)
		self.perf = self.res.calc_perf()

	def test_sched(self):
		# Note these are set to the values for what is thought to be a working
		# version.  They are not validated against anything or independently
		# calculated.
		self.assertTrue(abs(self.perf[0]- 746586.0)/746586.0 <0.01) # epy
		#self.assertTrue(abs(self.perf[1]-  158.11)/ 158.11<0.01) # LCOE
		#self.assertTrue(abs(self.perf[2]- 45.63)/45.63<0.01) # Capacity factor
		cleantest.clean('PBS_Surround_SCO2NREL_Table_Mixed')




if __name__ == '__main__':
	unittest.main()
