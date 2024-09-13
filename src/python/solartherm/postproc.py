from __future__ import division, print_function,unicode_literals

import solartherm.finances as fin
import DyMat
import xml.etree.ElementTree as ET
import re
import numpy as np

class SimResult(object):
	def __init__(self, fn, init_fn=None):
		"""
		fn: results file (*_res*.mat)
		init_fn: input paramet file (*_init*.xml)
		"""
		self.fn = fn
		self.init_fn = init_fn
		self.mat = None
		self.units = None
		self.load_res()
	
	def load_res(self):
		self.mat = DyMat.DyMatFile(self.fn)
		self.load_units()

	def load_units(self):
		self.units = None
		init_fn = self.init_fn

		if init_fn is None:
			# Find the .xml file with the name matching the .mat file.
			res_re = re.compile('(\S+)_res_?(\S+)?\.mat')
			ps = res_re.match(self.fn)
			if ps is not None:
				if ps.groups()[1] is not None:
					init_fn = ps.groups()[0] + '_init_' + ps.groups()[1] + '.xml'
				else:
					init_fn = ps.groups()[0] + '_init.xml'

		# Can fail
		et = None
		if init_fn is not None:
			try:
				et = ET.parse(init_fn)
			except IOError:
				pass

		if et is not None:
			root = et.getroot()

			self.units = {}
			for n in self.mat.names():
				ns = str(n)
				node = root.find('*ScalarVariable[@name=\''+ns+'\']/*[@unit]')
				self.units[ns] = '' if node is None else node.attrib['unit']
	
	def get_names(self):
		return self.mat.names()

	def get_time(self, name):
		return self.mat.abscissa(name, valuesOnly=True)

	def get_values(self, name):
		return self.mat.data(name)

	def get_unit(self, name):
		return self.units[name]

	def lower_ind(self, ab, t):
		"""Get index for point just below or equal to the requested time.
		"""
		# Note that values are provided before and after events for quantities
		assert ab[0] <= t <= ab[-1], "Time outside range"

		il = int(len(ab)*t/(ab[-1] - ab[0]))
		il = min(len(ab)-2, max(0, il))

		# FIXME make use of 'import bisect' here to increase speed
		while t < ab[il]:
			il -= 1
		while t > ab[il+1]:
			il += 1

		return il
	
	def closest(self, names, t):
		"""Closest point in interval for a list of variables placed in 'names' list as strings.
		"""
		if isinstance(names,str):
			names = [names]

		ab = self.mat.abscissa(names[0], valuesOnly=True)

		il = self.lower_ind(ab, t)
		iu = il + 1

		if t <= ((ab[iu] + ab[il])/2):
			return np.array([self.mat.data(name)[il] for name in names])
		else:
			return np.array([self.mat.data(name)[iu] for name in names])

	def interpolate(self, names, t):
		"""Linear interpolation of point for a list of variables placed in 'names' list as strings.
		"""
		if isinstance(names,str):
			names = [names]

		ab = self.mat.abscissa(names[0], valuesOnly=True)

		il = self.lower_ind(ab, t)
		iu = il + 1	

		vl = np.array([self.mat.data(name)[il] for name in names])
		vu = np.array([self.mat.data(name)[iu] for name in names])

		if ab[il] == ab[iu]:
			return vl
		else:
			return (vu - vl)*(t - ab[il])/(ab[iu] - ab[il]) + vl

	def integrate(self, names, t0, t1):
		"""Integration of linear interpolation (trapezoidal rule) over interval for a list of variables placed in 'names' list as strings.
		"""
		if isinstance(names,str):
			names = [names]

		ab = self.mat.abscissa(names[0], valuesOnly=True)

		val = np.empty((ab.shape[0],len(names)))
		for i,name in enumerate(names):
			val[:,i] = self.mat.data(name)

		il = self.lower_ind(ab, t0)
		iu = self.lower_ind(ab, t1) + 1

		vsum = np.zeros(len(names))
		for i in range(il, iu):
			vl = val[i,:]
			vu = val[i+1,:]

			tl = ab[i]
			tu = ab[i+1]

			if tl == tu:
				continue

			if tl < t0:
				vl = (vu - vl)*(t0 - tl)/(tu - tl) + vl
				tl = t0
			if tu > t1:
				vu = (vu - vl)*(t1 - tl)/(tu - tl) + vl
				tu = t1

			vsum += 0.5*(vu + vl)*(tu - tl)

		return vsum

	def sample(self, names, step):
		"""Sampling of a list of variables placed in 'names' list as strings,
		with 'step' timestep instead of the timestep in the original _res.mat file.
		"""
		if isinstance(names,str):
			names = [names]

		ab = self.mat.abscissa(names[0], valuesOnly=True)
		n = int((ab[-1] - ab[0])/step)

		t = np.empty(n)
		v = np.empty((n,len(names)))
		for i in range(n):
			t0 = step*i + ab[0]
			t1 = step*(i + 1) + ab[0]
			t[i] = (t0 + t1)/2
			v[i,:] = self.integrate(names, t0, t1)/step

		return t, v

	def constrained_optimisation(self):
		"""
		Showing if the optimisation process is constrained or not.
		"""
		try:
			constrained = self.mat.data('constrained') # trigger for constrained optimisation.
			distance_v = self.mat.data('distance') # distance added to a constant offset as penalty when a constraint is not respected.
			distance = distance_v[-1]

				# Note that 'constrained' must be defined in the Modelica model if there is a constraint, e.g. if T>T_0 then constrained=true; false otherwise.
			if constrained[-1] == 1.0:
				constr = True # i.e. the optimisation is constrained.
			else:
				constr = False # i.e. the optimisation is NOT constrained.
		except KeyError:
			constr = False
			distance = 0.

		return constr, distance

class SimResultElec(SimResult):
	def calc_perf(self, peaker=False):
		"""Calculate the solar power plant performance.
		Some of the metrics will be returned as none if simulation runtime is
		not a multiple of a year.

		peaker: bool, True: to calculate performance of a peaker plant
		"""
		var_names = self.get_names()
		assert('E_elec' in var_names), "For a levelised cost of electricity calculation, It is expected to see E_elec variable in the results file!"

		eng_t = self.mat.abscissa('E_elec', valuesOnly=True) # Time [s]
		eng_v = self.mat.data('E_elec') # Cumulative electricity generated [J]
		cap_v = self.mat.data('C_cap') # Capital costs [$]
		om_y_v = self.mat.data('C_year') # O&M costs per year [$/year]
		om_p_v = self.mat.data('C_prod') # O&M costs per production per year [$/J/year]
		disc_v = self.mat.data('r_disc') # Discount rate [-]
		life_v = self.mat.data('t_life') # Plant lifetime [year]
		cons_v = self.mat.data('t_cons') # Construction time [year]
		name_v = self.mat.data('P_name') # Generator nameplate [W]
		rev_v = self.mat.data('R_spot') # Cumulative revenue [$]

		dur = eng_t[-1] - eng_t[0] # Time duration [s]
		years = dur/31536000 # number of years of simulation [year]
		# Only provide certain metrics if runtime is a multiple of a year
		close_to_year = years > 0.5 and abs(years - round(years)) <= 0.01

		epy = fin.energy_per_year(dur, eng_v[-1]) # Energy expected in a year [J]
		srev = rev_v[-1] # spot market revenue [$]
		lcoe = None # Levelised cost of electricity
		capf = None # Capacity factor
		if close_to_year: 
			if peaker:
				tod_v=self.mat.data('TOD_W')
				tod_factor=tod_v[-1]/eng_v[-1]
				lcoe = fin.lcoe_p(cap_v[0], om_y_v[0] + om_p_v[0]*epy, disc_v[0],
						int(life_v[0]), int(cons_v[0]), epy, tod_factor)
				capf = fin.capacity_factor(name_v[0], tod_v[-1])

			else:
				lcoe = fin.lcoe_r(cap_v[0], om_y_v[0] + om_p_v[0]*epy, disc_v[0],
						int(life_v[0]), int(cons_v[0]), epy)
				capf = fin.capacity_factor(name_v[0], epy)

		# Convert to useful units
		epy = epy/(1e6*3600) # Convert from J/year to MWh/year
		if close_to_year: 
			lcoe = lcoe*1e6*3600 # Convert from $/J to $/MWh
			capf = 100*capf

		return [epy, lcoe, capf, srev,]


	def cost_breakdown(self):
		"""Calculate costs breakdown for the solar power plant"""
		eng_t = self.mat.abscissa('E_elec', valuesOnly=True) # Time [s]
		eng_v = self.mat.data('E_elec') # Cumulative electricity generated [J]
		disc_v = self.mat.data('r_disc') # Discount rate [-]
		life_v = self.mat.data('t_life') # Plant lifetime [year]

		C_field_v = self.mat.data('C_field') # Solar field capital cost [$]
		C_tower_v = self.mat.data('C_tower') # Tower capital cost [$]
		C_receiver_v = self.mat.data('C_receiver') # Receiver capital cost [$]
		C_storage_v = self.mat.data('C_storage') # Storage tanks capital cost [$]
		C_block_v = self.mat.data('C_block') # Power block capital cost [$]
		C_bop_v = self.mat.data('C_bop') # Balance of plant cost [$]
		C_land_v = self.mat.data('C_land') # Land cos [$]

		C_cap_v = self.mat.data('C_cap') # Capital costs [$]

		om_y_v = self.mat.data('C_year') # Fixed O&M costs per year [$/year]
		om_p_v = self.mat.data('C_prod') # Variable O&M costs per production per year [$/J/year]

		dur = eng_t[-1] - eng_t[0] # Time duration [s]
		epy = fin.energy_per_year(dur, eng_v[-1]) # Energy expected in a year [J]

		C_cap = C_cap_v[0] * 1e-3 # Total capital investment (TCI) [k$]
		C_cap_ann = fin.annualised_capital_cost(C_cap, disc_v[0], int(life_v[0])) # Annualised TCI [k$/year]
		C_year = (om_y_v[0] + om_p_v[0]*epy) * 1e-3 # Total operational costs [k$/year]

		C_cap_bd_n = ['Solar field', 'Tower', 'Receiver', 'Storage', 'PB', 'BOP', 'Land'] # Capital cost components name
		C_cap_bd_u = 'k$' # Capital cost components unit
		C_cap_bd_v = [C_field_v[0]*1e-3, C_tower_v[0]*1e-3, C_receiver_v[0]*1e-3, C_storage_v[0]*1e-3, C_block_v[0]*1e-3, C_bop_v[0]*1e-3, C_land_v[0]*1e-3] # Capital cost breakdown [k$]

		C_op_bd_n = ['Fixed O&M', 'variable O&M'] # Operational cost components name
		C_op_bd_u = 'k$/year' # Operational cost components unit
		C_op_bd_v = [om_y_v[0]*1e-3, om_p_v[0]*epy*1e-3] # Operational cost breakdown [k$/year]

		C_ann_bd_n = ['Total capital investment', 'Operational'] # Annualised cost breakdown names
		C_ann_bd_u = 'k$/year' # Annualised cost breakdown unit
		C_ann_bd_v = [C_cap_ann, C_year] # Annualised cost breakdown [k$/year]

		return C_cap_bd_n, C_cap_bd_u, C_cap_bd_v, C_op_bd_n, C_op_bd_u, C_op_bd_v, C_ann_bd_n, C_ann_bd_u, C_ann_bd_v

	perf_n = ['epy', 'lcoe', 'capf', 'srev']
	perf_u = ['MWh/year', '$/MWh', '%', '$']


class SimResultFuel(SimResult):
	def calc_perf(self):
		"""Calculate solar fuels plant performance.
		Some of the metrics will be returned as none if simulation runtime is
		not a multiple of a year.
		"""
		var_names = self.get_names()

		assert('V_fuel' in var_names), "For a fuel_calc calculation, it is expected to see V_fuel variable in the results file!"

		V_fuel_t = self.mat.abscissa('V_fuel', valuesOnly=True) # Time [s]
		V_fuel_v = self.mat.data('V_fuel') # Cumulative produced fuel [m3]
		C_cap_v = self.mat.data('C_cap') # Capital costs [$]
		C_labor_v = self.mat.data('C_labor') # Labor cost [$/year]
		C_catalyst_v = self.mat.data('C_catalyst') # Catalysts cost [$/year]
		C_om_v = self.mat.data('C_om') # Maintenance cost [$/year]
		C_water_v = self.mat.data('C_water') # Cost of water [$/year]
		C_algae_v = self.mat.data('C_algae') # Cost of algae [$/year]
		C_H2_v = self.mat.data('C_H2')  # Cost of hydrogen [$/year]
		C_CO2_v = self.mat.data('C_CO2')  # Cost of CO2 emissions [$/year]
		C_O2_v = self.mat.data('C_O2')  # Cost of Oxygen to sell [$/year]
		C_elec_v = self.mat.data('C_elec') # Cost of electricity consumption [$/year]
		disc_v = self.mat.data('r_disc') # Discount rate [-]
		infl_v = self.mat.data('r_i') # Inflation rate [-]
		life_v = self.mat.data('t_life') # Plant lifetime [year]
		cons_v = self.mat.data('t_cons') # Construction time [year]
		name_v = self.mat.data('v_flow_fuel_des') # Nominal fuel volumetric flow rate [m3/s]
		rev_v = self.mat.data('R_spot') # cumulative revenue

		C_op_v = C_water_v + C_algae_v + C_H2_v + C_CO2_v - C_O2_v + C_elec_v # Operating costs [$/year]
		C_year = C_labor_v[0] + C_catalyst_v[0] + C_om_v[0] + C_op_v[-1] # Total operational costs [$/year]

		dur = V_fuel_t[-1] - V_fuel_t[0] # Time duration [s]
		years = dur/31536000 # number of years of simulation [year]
		# Only provide certain metrics if runtime is a multiple of a year
		close_to_year = years > 0.5 and abs(years - round(years)) <= 0.01

		fpy = fin.fuel_per_year(dur, V_fuel_v[-1]) # Fuel produced in a year [m3/year]
		srev = rev_v[-1] # Spot market revenue [$/year]
		lcof = None # Levelised cost of fuel
		capf = None # Capacity factor
		if close_to_year: 
			lcof = fin.lcof_r(C_cap_v[0], C_year, disc_v[0], int(life_v[0]), int(cons_v[0]), fpy)
			capf = fin.capacity_factor_f(name_v[0], fpy)

		# Convert to useful units
		fpy = fpy*1e3 # convert from m3/year to L/year
		if close_to_year: 
			lcof = lcof/1e3 # convert from $/m3 to $/L
			capf = 100*capf

		return [fpy, lcof, capf, srev,]

	def cost_breakdown(self):
		"""Calculate costs breakdown for the solar fuels plant"""
		disc_v = self.mat.data('r_disc') # Discount rate [-]
		life_v = self.mat.data('t_life') # Plant lifetime [year]

		C_field_v = self.mat.data('C_field') # Solar field capital cost [$]
		C_tower_v = self.mat.data('C_tower') # Tower capital cost [$]
		C_rx_v = self.mat.data('C_rx') # Receiver capital cost [$]
		C_st_v = self.mat.data('C_st') # Storage tanks capital cost [$]
		C_ft_v = self.mat.data('C_ds_rx') # Downstream (e.g. Fischer tropsch or methanol synthesis) reactor cost [$]
		C_land_v = self.mat.data('C_land') # Land cos [$]

		C_cap_v = self.mat.data('C_cap') # Capital costs [$]

		C_labor_v = self.mat.data('C_labor') # Labor cost [$/year]
		C_catalyst_v = self.mat.data('C_catalyst') # Catalysts cost [$/year]
		C_om_v = self.mat.data('C_om') # Maintenance cost [$/year]

		C_water_v = self.mat.data('C_water') # Cost of water [$/year]
		C_algae_v = self.mat.data('C_algae') # Cost of algae [$/year]
		C_H2_v = self.mat.data('C_H2')  # Cost of hydrogen [$/year]
		C_CO2_v = self.mat.data('C_CO2')  # Cost of CO2 emissions [$/year]
		C_O2_v = self.mat.data('C_O2')  # Cost of Oxygen to sell [$/year]
		C_elec_v = self.mat.data('C_elec') # Cost of electricity consumption [$/year]

		C_op_v = C_water_v + C_algae_v + C_H2_v + C_CO2_v - C_O2_v + C_elec_v # Operating costs [$/year]

		C_cap = C_cap_v[0] * 1e-6 # Total capital investment (TCI) [M$]
		C_cap_ann = fin.annualised_capital_cost(C_cap, disc_v[0], int(life_v[0])) # Annualised TCI [M$/year]
		C_year = (C_labor_v[0] + C_catalyst_v[0] + C_om_v[0] + C_op_v[-1]) * 1e-6 # Total operational costs [M$/year]

		C_cap_bd_n = ['Solar field', 'Tower', 'Reactors', 'Storage', 'FT', 'Land'] # Capital cost components name
		C_cap_bd_u = 'M$' # Capital cost components unit
		C_cap_bd_v = [C_field_v[0]*1e-6, C_tower_v[0]*1e-6, C_rx_v[0]*1e-6, C_st_v[0]*1e-6, C_ft_v[0]*1e-6, C_land_v[0]*1e-6] # Capital cost breakdown [M$]

		C_op_bd_n = ['Labor', 'Catalysts', 'O&M', 'Water', 'Algae', 'H2', 'CO2', 'Electricity'] # Operational cost components name
		C_op_bd_u = 'M$/year' # Operational cost components unit
		C_op_bd_v = [C_labor_v[0]*1e-6, C_catalyst_v[0]*1e-6, C_om_v[0]*1e-6, C_water_v[-1]*1e-6,
			C_algae_v[-1]*1e-6, C_H2_v[-1]*1e-6, C_CO2_v[-1]*1e-6, C_elec_v[-1]*1e-6] # Operational cost breakdown [M$/year]

		C_ann_bd_n = ['Total capital investment', 'Operational'] # Annualised cost breakdown names
		C_ann_bd_u = 'M$/year' # Annualised cost breakdown unit
		C_ann_bd_v = [C_cap_ann, C_year] # Annualised cost breakdown [M$/year]

		return C_cap_bd_n, C_cap_bd_u, C_cap_bd_v, C_op_bd_n, C_op_bd_u, C_op_bd_v, C_ann_bd_n, C_ann_bd_u, C_ann_bd_v

	perf_n = ['fpy', 'lcof', 'capf', 'srev']
	perf_u = ['L/year', '$/L', '%', '$']


class SimResultHeat(SimResult):
	def calc_perf(self):
		"""Calculate LCOH of renewable systems
		Some of the metrics will be returned as none if simulation runtime is
		not a multiple of a year.
		"""
		var_names = self.get_names()
		assert('E_supplied' in var_names), "For a levelised cost of heat calculation, It is expected to see E_elec variable in the results file!"

		eng_t = self.mat.abscissa('E_supplied', valuesOnly=True) # Time [s]
		eng_v = self.mat.data('E_supplied') # Cumulative electricity generated [J]
		cap_v = self.mat.data('FCI_total') # Capital costs [$]
		om_y_v = [0.] #self.mat.data('C_year') # O&M costs per year [$/year]
		om_p_v =[0.]# self.mat.data('C_prod') # O&M costs per production per year [$/J/year]
		disc_v = [0.064] #self.mat.data('r_disc') # Discount rate [-]
		life_v = [25] #self.mat.data('t_life') # Plant lifetime [year]
		cons_v = [0] #self.mat.data('t_cons') # Construction time [year]
		#name_v = self.mat.data('P_name') # Generator nameplate [W]
		#rev_v = self.mat.data('R_spot') # Cumulative revenue [$]

		dur = eng_t[-1] - eng_t[0] # Time duration [s]
		years = dur/31536000 # number of years of simulation [year]
		# Only provide certain metrics if runtime is a multiple of a year
		close_to_year = years > 0.5 and abs(years - round(years)) <= 0.01

		epy = fin.energy_per_year(dur, eng_v[-1]) # Energy expected in a year [J]
		#srev = rev_v[-1] # spot market revenue [$]
		lcoh = None # Levelised cost of heat
		capf = self.mat.data('Capacity_Factor')[-1] # Capacity factor

		lcoh = fin.lcoe_r(cap_v[0], om_y_v[0] + om_p_v[0]*epy, disc_v[0],int(life_v[0]), int(cons_v[0]), epy)
		#capf = fin.capacity_factor(name_v[0], epy)

		# Convert to useful units
		epy = epy/(1e6*3600) # Convert from J/year to MWh/year
		if close_to_year: 
			lcoh = lcoh*1e6*3600 # Convert from $/J to $/MWh
			capf = 100*capf

		return [epy, lcoh, capf]


	def cost_breakdown(self):
		"""Calculate costs breakdown for the renewable system"""
		eng_t = self.mat.abscissa('E_supplied', valuesOnly=True) # Time [s]
		eng_v = self.mat.data('E_supplied') # Cumulative electricity generated [J]
		disc_v = [0.064] #self.mat.data('r_disc') # Discount rate [-]
		life_v = [25] #self.mat.data('t_life') # Plant lifetime [year]

		C_filler_v = self.mat.data('FCI_filler') # filler capital cost [$]
		C_insulation_v = self.mat.data('FCI_insulation') # tank insulation capital cost [$]
		C_tank_v = self.mat.data('FCI_tank') # tank capital cost [$]
		C_heater_v = self.mat.data('FCI_heater') # heater capital cost [$]
		C_HX_v = self.mat.data('FCI_HX') # HX cost [$]
		C_blower1_v = self.mat.data('FCI_blower1') # Blower 1 cost [$]
		C_blower2_v = self.mat.data('FCI_blower2') # Blower 2 cost [$]
	
		C_pipe1_v = self.mat.data('FCI_piping1') # piping cost 1 [$]
		C_pipe2_v = self.mat.data('FCI_piping2') # piping cost 2 [$]
		C_pipe3_v = self.mat.data('FCI_piping3') # piping cost 3 [$]
		C_pipe4_v = self.mat.data('FCI_piping4') # piping cost 4 [$]
		C_pipe5_v = self.mat.data('FCI_piping5') # piping cost 5 [$]
		C_PV_v = self.mat.data('FCI_PV') # PV cost [$]
		C_wind_v = self.mat.data('FCI_wind') # wind cost [$]

		C_cap_v = self.mat.data('FCI_total') # Capital costs [$]

		C_HBS=C_filler_v[0]+C_insulation_v[0]+C_tank_v[0]
		C_blowers=C_blower1_v[0]+C_blower2_v[0]
		C_pip=C_pipe1_v[0]+C_pipe2_v[0]+C_pipe3_v[0]+C_pipe4_v[0]+C_pipe5_v[0]


		om_y_v = [0]# self.mat.data('C_year') # Fixed O&M costs per year [$/year]
		om_p_v = [0]# self.mat.data('C_prod') # Variable O&M costs per production per year [$/J/year]

		dur = eng_t[-1] - eng_t[0] # Time duration [s]
		epy = fin.energy_per_year(dur, eng_v[-1]) # Energy expected in a year [J]

		C_cap = C_cap_v[0] * 1e-3 # Total capital investment (TCI) [k$]
		C_cap_ann = fin.annualised_capital_cost(C_cap, disc_v[0], int(life_v[0])) # Annualised TCI [k$/year]
		C_year = (om_y_v[0] + om_p_v[0]*epy) * 1e-3 # Total operational costs [k$/year]

		C_cap_bd_n = ['HBS', 'Heater', 'HX', 'Blower1', 'Blower2', 'Piping', 'PV', 'Wind'] # Capital cost components name
		C_cap_bd_u = 'k$' # Capital cost components unit
		C_cap_bd_v = [C_HBS*1e-3, C_heater_v[0]*1e-3, C_HX_v[0]*1e-3, C_blower1_v[0]*1e-3, C_blower2_v[0]*1e-3, C_pip*1e-3, C_PV_v[0]*1e-3, C_wind_v[0]*1e-3] # Capital cost breakdown [k$]

		C_op_bd_n = ['Fixed O&M', 'variable O&M'] # Operational cost components name
		C_op_bd_u = 'k$/year' # Operational cost components unit
		C_op_bd_v = [om_y_v[0]*1e-3, om_p_v[0]*epy*1e-3] # Operational cost breakdown [k$/year]

		C_ann_bd_n = ['Total capital investment', 'Operational'] # Annualised cost breakdown names
		C_ann_bd_u = 'k$/year' # Annualised cost breakdown unit
		C_ann_bd_v = [C_cap_ann, C_year] # Annualised cost breakdown [k$/year]

		return C_cap_bd_n, C_cap_bd_u, C_cap_bd_v, C_op_bd_n, C_op_bd_u, C_op_bd_v, C_ann_bd_n, C_ann_bd_u, C_ann_bd_v

	perf_n = ['epy', 'lcoe', 'capf']
	perf_u = ['MWh/year', '$/MWh', '%']

class DecisionMaker(object):
	"""
	This class is used to select the final optimum point amongs a list of optimal solutions in a multi-objective optimisation problem.
	"""
	def __init__(self, cand, fitness):
		self.cand = cand
		self.fitness = fitness

	def linmap(self):
		""" LINMAP decision-making method for multi-objective optimisation. This method is used to select a final optimal point amongs a list of optimal solutions.
		It is assumed that obj1(e.g. lcoe) was minimised and obj2(e.g. epy) was maximised.
		cand is the list of optimal individuals (design variables) and fitness is the list of optimal fitness (objective functions), where each objective set (i.e. obj1 and obj2) are nested in the list as a tuple
		"""
		obj1 = np.array([vv[0] for vv in self.fitness]) # Objective 1 that has been minimised
		obj2 = np.array([vv[1] for vv in self.fitness]) # Objective 2 that has been maximised

		obj1_n = obj1 / np.sqrt(np.sum(obj1**2)) # Normalised version of objective 1
		obj2_n = obj2 / np.sqrt(np.sum(obj2**2)) # Normalised version of objective 2

		obj1_n_min = obj1_n.min() # The minimum value of normalised objective 1 (i.e the x coordinate for the ideal point)
		i1, = np.where(obj1_n==obj1_n_min) # x coordinate index of the ideal point

		obj2_n_max = obj2_n.max() # The maximum value of normalised objective 2 (i.e the y coordinate for the ideal point)
		i2, = np.where(obj2_n==obj2_n_max) # y coordinate index of the ideal point

		d_plus = np.sqrt((obj1_n - obj1_n[i1[0]])**2 + (obj2_n - obj2_n[i2[0]])**2) # Distance of each optimal solution from the ideal point
		d_plus_min = d_plus.min() # Minimum distance form the ideal point
		ind_opt, = np.where(d_plus==d_plus_min) # Index of the point that has the minimum distance from the ideal point

		best_ind = self.cand[ind_opt[0]] # Final optimal solution
		best_fitness = list(self.fitness[ind_opt[0]]) # Final optimal objective

		return best_ind, best_fitness

	def topsis(self):
		""" TOPSIS decision-making method for multi-objective optimisation. This method is used to select a final optimal point amongs a list of optimal solutions.
		It is assumed that obj1(e.g. lcoe) was minimised and obj2(e.g. epy) was maximised.
		cand is the list of optimal individuals (design variables) and fitness is the list of optimal fitness (objective functions), where each objective set (i.e. obj1 and obj2) are nested in the list as a tuple
		"""
		obj1 = np.array([vv[0] for vv in self.fitness]) # Objective 1 that has been minimised
		obj2 = np.array([vv[1] for vv in self.fitness])  # Objective 2 that has been maximised

		obj1_n = obj1 / np.sqrt(np.sum(obj1**2))  # Normalised version of objective 1
		obj2_n = obj2 / np.sqrt(np.sum(obj2**2))  # Normalised version of objective 12 

		obj1_n_min = obj1_n.min() # The minimum value of normalised objective 1 (i.e the x coordinate for the ideal point)
		i1, = np.where(obj1_n==obj1_n_min) # x coordinate index of the ideal point

		obj2_n_max = obj2_n.max() # The maximum value of normalised objective 2 (i.e the y coordinate for the ideal point)
		i2, = np.where(obj2_n==obj2_n_max) # y coordinate index of the ideal point

		obj1_n_max = obj1_n.max() # The maximum value of normalised objective 1 (i.e the x coordinate for the mon-ideal point)
		ni1, = np.where(obj1_n==obj1_n_max) # x coordinate index of the non-ideal point

		obj2_n_min = obj2_n.min() # The minimum value of normalised objective 2 (i.e the y coordinate for the non-ideal point)
		ni2, = np.where(obj2_n==obj2_n_min) # y coordinate index of the non-ideal point

		d_plus = np.sqrt((obj1_n - obj1_n[i1[0]])**2 + (obj2_n - obj2_n[i2[0]])**2) # Distance of each optimal solution from the ideal point
		d_minus = np.sqrt((obj1_n - obj1_n[ni1[0]])**2 + (obj2_n - obj2_n[ni2[0]])**2) # Distance of each optimal solution from the non-ideal point

		cl = d_minus/(d_plus+d_minus)
		cl_max = cl.max()

		ind_opt, = np.where(cl==cl_max) # Index of the point that has the minimum distance from the ideal point and maximum distance from the non-ideal point
		best_ind = self.cand[ind_opt[0]] # Final optimal solution
		best_fitness = list(self.fitness[ind_opt[0]]) # Final optimal objective

		return best_ind, best_fitness

class CSVResult(object):
	"""Results from a CSV file.
	Assumes first row is a header.
	Assumes first column is time.
	"""
	def __init__(self, fn, delim=','):
		self.fn = fn
		self.delim = delim
		self.names = []
		self.data = {}
		self.units = {}
		self.load_res()
	
	def load_res(self):
		f = open(self.fn)
		header = f.readline().split(self.delim)
		for label in header:
			l = label.split('(')
			name = l[0].strip()
			self.names.append(name)
			if len(l) > 1:
				self.units[name] = l[1].strip().strip(')')
			else:
				self.units[name] = ''
			self.data[name] = []

		for line in f.readlines():
			l = line.split(self.delim)
			for i, val in enumerate(l):
				self.data[self.names[i]].append(float(val))

	def get_names(self):
		return self.names

	def get_time(self, name):
		return self.data[self.names[0]] # assume first

	def get_values(self, name):
		return self.data[name]

	def get_unit(self, name):
		return self.units[name]
