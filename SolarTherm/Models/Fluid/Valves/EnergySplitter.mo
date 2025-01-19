within SolarTherm.Models.Fluid.Valves;
model EnergySplitter "Splits input between two output ports"
	import SI = Modelica.SIunits;
	input Real frac(min=0, max=1) "Fraction to output port 1";
	SolarTherm.Models.Fluid.Interfaces.EnergyPort p_i;
	SolarTherm.Models.Fluid.Interfaces.EnergyPort p_o1;
	SolarTherm.Models.Fluid.Interfaces.EnergyPort p_o2;
equation
	p_o1.P = -frac*p_i.P;
	p_o2.P = -(1 - frac)*p_i.P;
end EnergySplitter;