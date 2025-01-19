#! /bin/env python

from __future__ import division
import unittest

import cleantest
from solartherm import simulation
import DyMat

from math import pi
import os

class TestReference(unittest.TestCase):
	def setUp(self):
		fn = '../examples/WindPVAnnularTESsystem.mo'
		sim = simulation.Simulator(fn)
		sim.compile_model()
		sim.compile_sim(args=['-s'])
		sim.simulate(start=0, stop='1y', step='300s', tolerance='1e-06', maxStep='60s', initStep='60s', solver='dassl', nls='mixed', args=['-noEventEmit'])
		self.mat = DyMat.DyMatFile(sim.res_fn)

	def test_sched(self):
		# Note these are set to the values for what is thought to be a working
		# version.  They are not validated against anything or independently
		# calculated.
		ref_annual_output = 770491.67 # Reference annnual thermal output
		annual_output=self.mat.data('E_supplied')[-1]
		annual_output = annual_output/(1e6*3600) # Convert from J/year to MWh/year
		self.assertTrue(abs(annual_output- ref_annual_output)/ref_annual_output<0.01)
		cleantest.clean('WindPVAnnularTESsystem')




if __name__ == '__main__':
	unittest.main()
