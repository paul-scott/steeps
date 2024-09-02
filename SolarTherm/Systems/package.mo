within SolarTherm;
package Systems

record Parameters
	extends Modelica.Icons.Record;
	parameter Real m_DRI = 1;
end Parameters;

model HydrogenStorage
	extends SolarTherm.Icons.StorageModel;
	parameter SI.Efficiency eff_in = 0.99 "Efficiency of storing hydrogen";
	parameter SI.Efficiency eff_out = 0.99 "Efficiency of discharging hydrogen";

	Modelica.Blocks.Interfaces.RealInput E_in "H2 energy input" 
		annotation(
		Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), 
		iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

	Modelica.Blocks.Interfaces.RealOutput E_out "H2 energy output" 
		annotation(
		Placement(visible = true, transformation(origin = {100, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), 
		iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

    SI.Energy E "H2 Stored Chemical Energy";

equation

	der(E) = E_in * eff_in - E_out / eff_out;

end HydrogenStorage;
model Electrolyser
	extends SolarTherm.Icons.ResistiveHeater;
	Modelica.Blocks.Interfaces.RealInput P_elec_in "Electrical power input" 
		annotation(
		Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), 
		iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

	Modelica.Blocks.Interfaces.RealOutput P_thermal_out "Electrical power input" 
		annotation(
		Placement(visible = true, transformation(origin = {100, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), 
		iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

end Electrolyser;
model PumpHydroStorage
	import SI = Modelica.SIunits;
	parameter SI.Energy E_max = 1000 "Maximum energy";
	parameter SI.Energy E_0 = 0 "Starting energy";
	parameter SI.Efficiency eta_phes "Design point storage efficiency";

	Modelica.Blocks.Interfaces.RealInput P_elec_in "Thermal power input" 
		annotation(
		Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), 
		iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

	Modelica.Blocks.Interfaces.RealOutput P_elec_out "Thermal power output" 
		annotation(
		Placement(visible = true, transformation(origin = {80, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), 
		iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

	Modelica.Blocks.Interfaces.RealOutput Level annotation(
		Placement(visible = true, transformation(origin = {80, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 0), 
		iconTransformation(origin = {100, -60}, extent = {{-10,-10},{10, 10}}, rotation = 0)));

	SI.Energy E(start=E_0, fixed=true, min=0, max=E_max) "Energy in tank";
	SI.HeatFlowRate P_elec_loss "Heat loss";
initial equation
	Level = E_0/E_max*100;
equation
	der(E) = P_elec_in - P_elec_out - P_elec_loss;
	P_elec_loss = P_elec_in*eta_phes;
	Level = E/E_max*100;

annotation(
        Icon(
            coordinateSystem(
                extent= {{-100,-100},{100,100}},
                preserveAspectRatio= false
            ),
            graphics= {
                Rectangle(
                    extent= {{652.18,-492.16},{675.31,-501.96}},
                    fillColor= {128,128,128},
                    fillPattern= FillPattern.Solid,
                    lineThickness= 0.25,
                    origin= {-650.36,523.24}
                ),
                Polygon(
                    fillColor= {128,179,255},
                    fillPattern= FillPattern.Solid,
                    lineThickness= 0.25,
                    origin= {-650.36,523.24},
                    points= {{570.85, -466.79}, {618.41, -466.76}, {615.31, -478.33}, {600.98, -478.33}, {600.59, -515.48}, {626.66, -546.72}, {633.97, -546.72}, {643.72, -536.74}, {647.53, -540.58}, {659.75, -540.58}, {663.31, -536.80}, {673.02, -546.51}, {707.18, -546.51}, {707.18, -538.40}, {741.22, -538.40}, {741.22, -568.89}, {719.34, -568.89}, {706.62, -557.95}, {672.46, -557.95}, {667.55, -567.33}, {640.09, -567.55}, {635.63, -556.39}, {621.79, -556.39}, {585.62, -520.22}, {585.62, -478.33}, {573.94, -478.33}}
                ),
                Polygon(
                    fillColor= {222,170,135},
                    fillPattern= FillPattern.Solid,
                    lineThickness= 0.25,
                    origin= {-650.36,523.24},
                    points= {{559.50, -455.26}, {567.76, -455.26}, {573.94, -478.33}, {585.62, -478.33}, {585.62, -520.22}, {621.79, -556.39}, {635.63, -556.39}, {640.09, -567.55}, {667.55, -567.33}, {672.46, -557.95}, {706.62, -557.95}, {719.34, -568.89}, {741.22, -568.89}, {741.22, -591.22}, {559.50, -591.22}}
                ),
                Polygon(
                    fillColor= {222,170,135},
                    fillPattern= FillPattern.Solid,
                    lineThickness= 0.25,
                    origin= {-650.36,523.24},
                    points= {{615.31, -478.33}, {621.51, -455.18}, {634.51, -455.18}, {643.91, -460.78}, {647.35, -487.13}, {659.85, -502.08}, {675.43, -502.08}, {680.35, -507.11}, {686.89, -524.62}, {707.18, -530.28}, {707.18, -546.51}, {673.02, -546.51}, {663.31, -536.80}, {643.78, -536.80}, {633.97, -546.72}, {626.66, -546.72}, {600.59, -518.96}, {600.98, -478.33}}
                ),
                Polygon(
                    fillColor= {204,204,204},
                    fillPattern= FillPattern.Solid,
                    lineThickness= 0.25,
                    origin= {-650.36,523.24},
                    points= {{671.34, -523.32}, {671.34, -544.83}, {663.31, -536.80}, {663.31, -523.38}}
                ),
                Polygon(
                    fillColor= {204,204,204},
                    fillPattern= FillPattern.Solid,
                    lineThickness= 0.25,
                    origin= {-650.36,523.24},
                    points= {{643.72, -536.74}, {647.53, -540.58}, {659.75, -540.58}, {663.31, -536.80}, {663.31, -523.38}, {643.72, -523.32}}
                ),
                Polygon(
                    fillColor= {153,153,153},
                    fillPattern= FillPattern.Solid,
                    lineThickness= 0.25,
                    origin= {-650.36,523.24},
                    points= {{645.13, -547.97}, {645.13, -557.44}, {662.16, -557.44}, {662.16, -547.97}}
                ),
                Line(
                    origin= {-650.36,523.24},
                    points= {{653.64, -547.97}, {653.64, -540.58}},
                    thickness= 0.25
                ),
                Polygon(
                    fillColor= {204,204,204},
                    fillPattern= FillPattern.Solid,
                    lineThickness= 0.25,
                    origin= {-650.36,523.24},
                    points= {{635.69, -523.26}, {635.69, -544.76}, {643.72, -536.74}, {643.72, -523.32}}
                ),
                Polygon(
                    fillColor= {153,153,153},
                    fillPattern= FillPattern.Solid,
                    lineThickness= 0.25,
                    origin= {-650.36,523.24},
                    points= {{639.70, -523.29}, {643.72, -509.43}, {663.31, -509.43}, {667.33, -523.35}}
                ),
                Ellipse(
                    extent= {{566.26,-495.52},{577.15,-500.38}},
                    fillColor= {160,90,44},
                    fillPattern= FillPattern.Solid,
                    origin= {-650.36,523.24},
                    pattern= LinePattern.None
                ),
                Ellipse(
                    extent= {{577.97,-529.25},{588.86,-534.12}},
                    fillColor= {160,90,44},
                    fillPattern= FillPattern.Solid,
                    origin= {-650.36,523.24},
                    pattern= LinePattern.None
                ),
                Ellipse(
                    extent= {{562.66,-517.19},{573.56,-522.06}},
                    fillColor= {160,90,44},
                    fillPattern= FillPattern.Solid,
                    origin= {-650.36,523.24},
                    pattern= LinePattern.None
                ),
                Ellipse(
                    extent= {{568.92,-559.62},{579.82,-564.49}},
                    fillColor= {160,90,44},
                    fillPattern= FillPattern.Solid,
                    origin= {-650.36,523.24},
                    pattern= LinePattern.None
                ),
                Ellipse(
                    extent= {{566.84,-578.40},{577.73,-583.27}},
                    fillColor= {160,90,44},
                    fillPattern= FillPattern.Solid,
                    origin= {-650.36,523.24},
                    pattern= LinePattern.None
                ),
                Ellipse(
                    extent= {{593.73,-554.52},{604.63,-559.39}},
                    fillColor= {160,90,44},
                    fillPattern= FillPattern.Solid,
                    origin= {-650.36,523.24},
                    pattern= LinePattern.None
                ),
                Ellipse(
                    extent= {{598.83,-577.70},{609.73,-582.57}},
                    fillColor= {160,90,44},
                    fillPattern= FillPattern.Solid,
                    origin= {-650.36,523.24},
                    pattern= LinePattern.None
                ),
                Ellipse(
                    extent= {{625.49,-567.97},{636.39,-572.84}},
                    fillColor= {160,90,44},
                    fillPattern= FillPattern.Solid,
                    origin= {-650.36,523.24},
                    pattern= LinePattern.None
                ),
                Ellipse(
                    extent= {{649.38,-580.26},{660.27,-585.12}},
                    fillColor= {160,90,44},
                    fillPattern= FillPattern.Solid,
                    origin= {-650.36,523.24},
                    pattern= LinePattern.None
                ),
                Ellipse(
                    extent= {{686.47,-569.82},{697.37,-574.69}},
                    fillColor= {160,90,44},
                    fillPattern= FillPattern.Solid,
                    origin= {-650.36,523.24},
                    pattern= LinePattern.None
                ),
                Ellipse(
                    extent= {{672.79,-581.65},{683.69,-586.52}},
                    fillColor= {160,90,44},
                    fillPattern= FillPattern.Solid,
                    origin= {-650.36,523.24},
                    pattern= LinePattern.None
                ),
                Ellipse(
                    extent= {{604.86,-498.65},{615.76,-503.51}},
                    fillColor= {160,90,44},
                    fillPattern= FillPattern.Solid,
                    origin= {-650.36,523.24},
                    pattern= LinePattern.None
                ),
                Ellipse(
                    extent= {{620.86,-516.27},{631.75,-521.13}},
                    fillColor= {160,90,44},
                    fillPattern= FillPattern.Solid,
                    origin= {-650.36,523.24},
                    pattern= LinePattern.None
                ),
                Ellipse(
                    extent= {{621.55,-492.15},{632.45,-497.02}},
                    fillColor= {160,90,44},
                    fillPattern= FillPattern.Solid,
                    origin= {-650.36,523.24},
                    pattern= LinePattern.None
                ),
                Ellipse(
                    extent= {{629.90,-469.20},{640.80,-474.07}},
                    fillColor= {160,90,44},
                    fillPattern= FillPattern.Solid,
                    origin= {-650.36,523.24},
                    pattern= LinePattern.None
                ),
                Ellipse(
                    extent= {{640.56,-499.57},{651.46,-504.44}},
                    fillColor= {160,90,44},
                    fillPattern= FillPattern.Solid,
                    origin= {-650.36,523.24},
                    pattern= LinePattern.None
                ),
                Ellipse(
                    extent= {{684.15,-536.90},{695.05,-541.77}},
                    fillColor= {160,90,44},
                    fillPattern= FillPattern.Solid,
                    origin= {-650.36,523.24},
                    pattern= LinePattern.None
                ),
                Line(
                    origin= {-650.36,523.24},
                    points= {{653.51, -509.43}, {659.75, -501.96}},
                    thickness= 0.25
                ),
                Line(
                    origin= {-650.36,523.24},
                    points= {{662.49, -492.16}, {662.49, -487.46}},
                    thickness= 0.25
                ),
                Line(
                    origin= {-650.36,523.24},
                    points= {{668.06, -492.16}, {668.06, -487.46}},
                    thickness= 0.25
                ),
                Text(
		origin = {0, 0}, 
		lineColor = {0, 0, 255}, 
		extent = {{-100, -110}, {100, -70}}, 
		textString = "%name")
            }
        )
    );
end PumpHydroStorage;
model GridInputSplit
	extends SolarTherm.Icons.GridInput;
	parameter String pv_file "File with the reference PV farm output";
	parameter String wind_file "File with the reference Wind farm output";
	parameter Modelica.SIunits.Power P_elec_min = 25e6 "Minimum power input";
	parameter Modelica.SIunits.Power P_elec_max = 100e6 "Maximum power input";
	parameter Modelica.SIunits.Efficiency pv_fraction = 0.5 "PV fraction of renewable input at design";
	parameter Modelica.SIunits.Power P_elec_pv_ref_size = 50e6 "PV farm reference size";
	parameter Modelica.SIunits.Power P_elec_wind_ref_size = 50e6 "Wind farm reference size";
	parameter Real renewable_multiple = 2 "Oversizing factor with respect to the heat input";
	final parameter Modelica.SIunits.Power P_elec_max_pv = pv_fraction*renewable_multiple*P_elec_max "Maximum PV capacity";
	final parameter Modelica.SIunits.Power P_elec_max_wind = (1-pv_fraction)*renewable_multiple*P_elec_max "Maximum Wind farm capacity";

	Modelica.Blocks.Sources.CombiTimeTable P_elec_ref_pv(
		fileName = pv_file, 
		tableName = "Power", 
		tableOnFile = true,
		smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative);
	Modelica.Blocks.Sources.CombiTimeTable P_elec_rec_wind(
		fileName = wind_file, 
		tableName = "Power", 
		tableOnFile = true,
		smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative);

	Modelica.Blocks.Sources.RealExpression P_elec_net_switch1(y = P_elec_net) annotation(
		Placement(visible = true, transformation(origin = {-50, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

	Modelica.Blocks.Sources.RealExpression P_elec_curtail_switch1(y = min(P_schedule,P_elec_net)) annotation(
		Placement(visible = true, transformation(origin = {-50, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

	Modelica.Blocks.Logical.Switch switch1 annotation(
		Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));

	Modelica.Blocks.Interfaces.BooleanInput curtail annotation(
		Placement(visible = true, 
		transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), 
		iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

	Modelica.Blocks.Interfaces.RealInput P_schedule annotation(
		Placement(visible = true, 
		transformation(origin = {-100, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0), 
		iconTransformation(origin = {-100, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

	Modelica.Blocks.Interfaces.RealOutput P_elec_out annotation(
		Placement(visible = true, 
		transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), 
		iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

	Modelica.Blocks.Interfaces.BooleanOutput on_renewable annotation (
		Placement(transformation(origin={0,-114}, extent={{20,-20},{-20,20}},rotation=90),
		iconTransformation(origin = {0,-114}, extent = {{10,-10},{-10,10}}, rotation = 90)));

	Modelica.SIunits.Power P_elec_net;
	Modelica.SIunits.Power P_elec_raw;
	Modelica.SIunits.Power P_elec_raw_pv;
	Modelica.SIunits.Power P_elec_raw_wind;

initial equation
	on_renewable = P_elec_raw > P_elec_min;

equation
	P_elec_raw_pv = P_elec_max_pv / P_elec_pv_ref_size * P_elec_ref_pv.y[1];
	P_elec_raw_wind = P_elec_max_wind / P_elec_wind_ref_size * P_elec_rec_wind.y[1];
	
	P_elec_raw = min(P_elec_raw_pv + P_elec_raw_wind, P_elec_max);
	
	on_renewable = P_elec_raw > P_elec_min;
	
	P_elec_net = if on_renewable then P_elec_raw else 0.0;

	connect(P_elec_curtail_switch1.y, switch1.u1) annotation(
		Line(points = {{-39, -30}, {-25.5, -30}, {-25.5, -8}, {-12, -8}}, color = {0, 0, 127}));
	connect(curtail, switch1.u2) annotation(
		Line(points = {{-100, 0}, {-12, 0}}, color = {255, 0, 255}));
	connect(P_elec_net_switch1.y, switch1.u3) annotation(
		Line(points = {{-38, 30}, {-26, 30}, {-26, 8}, {-12, 8}}, color = {0, 0, 127}));
	connect(switch1.y, P_elec_out) annotation(
		Line(points = {{12, 0}, {100, 0}}, color = {0, 0, 127}));

annotation(
	Documentation(info="<html>
	<p>
	<b>GridInput</b> models the input of electrical power from renewable sources, such as photovoltaic (PV) and wind farms, into an electrical grid. The model calculates the net electrical power input based on reference power data, curtailment, and user-defined parameters.
	</p>
	<p>
	The <b>GridInput</b> model has the following connectors and parameters:
	</p>
	<ul>
	<li> Parameters:
	<ul>
	<li> <b>pv_file</b>: File with the reference PV farm output.</li>
	<li> <b>wind_file</b>: File with the reference Wind farm output.</li>
	<li> <b>P_elec_min</b>: Minimum power input, in Watts. Default: 25e6 W.</li>
	<li> <b>P_elec_max</b>: Maximum power input, in Watts. Default: 100e6 W.</li>
	<li> <b>pv_fraction</b>: PV fraction of renewable input at design. Default: 0.5.</li>
	<li> <b>P_elec_pv_ref_size</b>: PV farm reference size, in Watts. Default: 50e6 W.</li>
	<li> <b>P_elec_wind_ref_size</b>: Wind farm reference size, in Watts. Default: 50e6 W.</li>
	</ul>
	</li>
	<li> Connectors:
	<ul>
	<li> <b>P_elec_net</b>: Net electrical power output to the grid.</li>
	<li> <b>curtail</b>: Input for curtailment control.</li>
	<li> <b>P_schedule</b>: Input for controlling power scheduling.</li>
	<li> <b>P_elec_out1</b>: Output for the electrical power input to the grid.</li>
	<li> <b>on_renewable</b>: Output indicating whether the power input is from renewable sources.</li>
	</ul>
	</li>
	</ul>
	</html>", revisions="<html>
	<ul>
	<li><i>September 2023</i> by <a href=\"mailto:armando.fontalvo@anu.edu.au\">Armando Fontalvo</a>:<br>
	Created documentation for GridInput.</li>
	</ul>
	</html>"));

end GridInputSplit;
model Beneficiation
Modelica.Blocks.Interfaces.RealInput E_in "Energy input" 
	annotation(
	Placement(visible = true, transformation(origin = {-100, 50}, extent = {{-20, -20}, {20, 20}}, rotation = 0), 
	iconTransformation(origin = {-100, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

Modelica.Blocks.Interfaces.RealInput E_schedule annotation (
	Placement(visible = true,transformation(origin={-50,108},extent={{-16,-16},{16,16}},rotation=-90), 
	iconTransformation(origin={-101,1},extent={{-11,-11},{11,11}},rotation=0)));

Modelica.Blocks.Interfaces.RealOutput y annotation(
      Placement(visible = true, transformation(origin = {100, 1.9984e-15}, extent = {{-20, -20}, {20, 20}}, rotation = 0), 
      iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

equation

E_in + E_schedule = 0;

annotation(
      Icon(graphics = {Rectangle(origin = {0, -20}, lineColor = {252, 175, 62}, fillColor = {252, 233, 79}, fillPattern = FillPattern.VerticalCylinder, extent = {{-60, 60}, {60, -60}}), Rectangle(origin = {0, 45}, fillColor = {211, 215, 207}, fillPattern = FillPattern.Solid, extent = {{-62, 5}, {62, -5}}), Rectangle(origin = {47, -88}, lineColor = {252, 175, 62}, fillColor = {252, 175, 62}, fillPattern = FillPattern.Solid, extent = {{-17, 8}, {17, -8}}), Rectangle(origin = {-27, 65}, fillColor = {211, 215, 207}, fillPattern = FillPattern.Solid, extent = {{-23, 15}, {23, -15}}), Rectangle(origin = {-25, 85}, lineColor = {252, 175, 62}, fillColor = {252, 233, 79}, fillPattern = FillPattern.VerticalCylinder, extent = {{-36, 5}, {36, -5}}), Rectangle(origin = {-55, 70}, fillColor = {211, 215, 207}, fillPattern = FillPattern.Solid, extent = {{-5, 10}, {5, -10}}), Rectangle(origin = {0, 65}, fillColor = {211, 215, 207}, fillPattern = FillPattern.Solid, extent = {{-4, 15}, {4, -15}}), Rectangle(origin = {0, 25}, lineColor = {193, 125, 17}, fillColor = {252, 233, 79}, fillPattern = FillPattern.VerticalCylinder, extent = {{-4, 15}, {4, -15}}), Rectangle(origin = {0, 9}, fillPattern = FillPattern.Solid, extent = {{-8, 1}, {8, -1}}), Rectangle(origin = {0, -20}, extent = {{-60, 60}, {60, -60}}), Rectangle(origin = {47, -88}, fillColor = {252, 175, 62}, extent = {{-17, 8}, {17, -8}}), Rectangle(origin = {-25, 85}, fillColor = {252, 233, 79}, extent = {{-36, 5}, {36, -5}}), Rectangle(origin = {65, -88}, fillPattern = FillPattern.Solid, extent = {{-1, 10}, {1, -10}}), Rectangle(origin = {-62, -73}, fillColor = {252, 175, 62}, fillPattern = FillPattern.Solid, extent = {{-2, 7}, {2, -7}}),
      Text(origin = {-10, 14}, lineColor = {0, 0, 255}, extent = {{-169, -114}, {171, -154}}, textString = "%name")}));

end Beneficiation;

model FBReactor
Modelica.Blocks.Interfaces.RealInput u1	annotation(
	Placement(visible = true, transformation(origin = {-100, 50}, extent = {{-20, -20}, {20, 20}}, rotation = 0), 
	iconTransformation(origin = {-100, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

Modelica.Blocks.Interfaces.RealInput u2 annotation (
	Placement(visible = true,transformation(origin={-100,-8.88178e-16},extent={{-20, -20},{20, 20}},rotation=0), 
	iconTransformation(origin={-99,1},extent={{-11,-11},{11,11}},rotation=0)));

Modelica.Blocks.Interfaces.RealInput u3	annotation(
	Placement(visible = true, transformation(origin = {-100, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 0), 
	iconTransformation(origin = {-100, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

Modelica.Blocks.Interfaces.RealOutput y annotation(
      Placement(visible = true, transformation(origin = {100, 1.9984e-15}, extent = {{-20, -20}, {20, 20}}, rotation = 0), 
      iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

equation

u1 + u2 = 0;

annotation(
      Icon(graphics = {
      Ellipse(origin = {0, 79}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-40, -21}, {40, 21}}), 
      Rectangle(lineColor = {85, 87, 83}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-40, 80}, {40, -80}}), 
      Polygon(origin = {0, -90}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, points = {{-44, 10}, {44, 10}, {40, -10}, {-40, -10}, {-44, 10}}), Rectangle(origin = {0, -58}, fillPattern = FillPattern.Solid, extent = {{-40, 2}, {40, -2}}), Rectangle(origin = {0, -18}, fillPattern = FillPattern.Solid, extent = {{-40, 2}, {40, -2}}), Rectangle(origin = {0, 22}, fillPattern = FillPattern.Solid, extent = {{-40, 2}, {40, -2}}), Rectangle(origin = {0, 32}, lineColor = {204, 0, 0}, fillColor = {245, 121, 0}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-40, 8}, {40, -8}}), Rectangle(origin = {0, -8}, lineColor = {204, 0, 0}, fillColor = {245, 121, 0}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-40, 8}, {40, -8}}), Rectangle(origin = {0, -48}, lineColor = {204, 0, 0}, fillColor = {245, 121, 0}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-40, 8}, {40, -8}}),
      Text(origin = {-10, 14}, lineColor = {0, 0, 255}, 	extent = {{-149, -114}, {151, -154}}, textString = "%name")}));

end FBReactor;

model HydrogenBurner
Modelica.Blocks.Interfaces.RealInput E_in "Energy input" 
	annotation(
	Placement(visible = true, transformation(origin = {-100, 50}, extent = {{-20, -20}, {20, 20}}, rotation = 0), 
	iconTransformation(origin = {-98, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

Modelica.Blocks.Interfaces.RealOutput y annotation(
      Placement(visible = true, transformation(origin = {100, 1.9984e-15}, extent = {{-20, -20}, {20, 20}}, rotation = 0), 
      iconTransformation(origin = {100, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

equation

E_in + E_schedule = 0;

annotation(
      Icon(graphics = {Rectangle(origin = {-60, -10}, fillColor = {239, 41, 41}, fillPattern = FillPattern.Solid, extent = {{-40, 50}, {40, -50}}), Rectangle(origin = {28, 20}, fillColor = {85, 87, 83}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-40, 20}, {40, -20}}), Rectangle(origin = {-16, 20}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-4, 24}, {4, -24}}), Rectangle(origin = {72, 20}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, extent = {{-4, 24}, {4, -24}}), Rectangle(origin = {88, 20}, lineColor = {211, 215, 207}, fillColor = {85, 87, 83}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-12, 20}, {12, -20}}), Rectangle(origin = {88, 20}, fillColor = {85, 87, 83}, extent = {{-12, 20}, {12, -20}}),
      Text(origin = {-10, 14}, lineColor = {0, 0, 255}, 	extent = {{-149, -114}, {151, -154}}, textString = "%name")
      }));

end HydrogenBurner;

model SmelterBOF
Modelica.Blocks.Interfaces.RealInput u1
	annotation(
	Placement(visible = true, transformation(origin = {-100, 50}, extent = {{-20, -20}, {20, 20}}, rotation = 0), 
	iconTransformation(origin = {-94, 80}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));

Modelica.Blocks.Interfaces.RealInput u2
	annotation(
	Placement(visible = true, transformation(origin = {-100, 20}, extent = {{-20, -20}, {20, 20}}, rotation = 0), 
	iconTransformation(origin = {-94, 20}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));

Modelica.Blocks.Interfaces.RealOutput y annotation(
      Placement(visible = true, transformation(origin = {100, 1.9984e-15}, extent = {{-20, -20}, {20, 20}}, rotation = 0), 
      iconTransformation(origin = {97, -61}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));

equation

y = u1;

annotation(
      Icon(graphics = {Rectangle(origin = {-38, 52}, fillColor = {186, 189, 182}, fillPattern = FillPattern.Solid, extent = {{-50, 40}, {50, -40}}), 
      Rectangle(origin = {-38, 38}, lineColor = {245, 121, 0}, fillColor = {204, 0, 0}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-50, 26}, {50, -26}}), 
      Rectangle(origin = {-82, 72}, fillColor = {136, 138, 133}, fillPattern = FillPattern.VerticalCylinder, extent = {{-2, 28}, {2, -28}}), 
      Rectangle(origin = {-70, 72}, fillColor = {136, 138, 133}, fillPattern = FillPattern.VerticalCylinder, extent = {{-2, 28}, {2, -28}}), 
      Rectangle(origin = {-56, 72}, fillColor = {136, 138, 133}, fillPattern = FillPattern.VerticalCylinder, extent = {{-2, 28}, {2, -28}}), 
      Rectangle(origin = {-44, 72}, fillColor = {136, 138, 133}, fillPattern = FillPattern.VerticalCylinder, extent = {{-2, 28}, {2, -28}}), 
      Rectangle(origin = {-30, 72}, fillColor = {136, 138, 133}, fillPattern = FillPattern.VerticalCylinder, extent = {{-2, 28}, {2, -28}}), 
      Rectangle(origin = {-18, 72}, fillColor = {136, 138, 133}, fillPattern = FillPattern.VerticalCylinder, extent = {{-2, 28}, {2, -28}}), 
      Rectangle(origin = {-6, 72}, fillColor = {136, 138, 133}, fillPattern = FillPattern.VerticalCylinder, extent = {{-2, 28}, {2, -28}}), 
      Rectangle(origin = {6, 72}, fillColor = {136, 138, 133}, fillPattern = FillPattern.VerticalCylinder, extent = {{-2, 28}, {2, -28}}), 
      Rectangle(origin = {-38, 38}, fillColor = {204, 0, 0}, extent = {{-50, 26}, {50, -26}}), 
      Polygon(origin = {50, -60}, fillColor = {136, 138, 133}, fillPattern = FillPattern.Solid, points = {{-30, 40}, {10, 40}, {30, 20}, {30, -34}, {10, -40}, {-30, -40}, {-50, -32}, {-50, 20}, {-30, 40}}), 
      Polygon(origin = {50, -60}, fillColor = {186, 189, 182}, fillPattern = FillPattern.Solid, points = {{-30, 36}, {10, 36}, {26, 20}, {26, -30}, {10, -36}, {-30, -36}, {-46, -28}, {-46, 20}, {-30, 36}}), 
      Polygon(origin = {50, -60}, lineColor = {204, 0, 0}, fillColor = {245, 121, 0}, fillPattern = FillPattern.HorizontalCylinder, points = {{-28, -6}, {8, -6}, {26, -6}, {26, -30}, {10, -36}, {-30, -36}, {-46, -28}, {-46, -6}, {-28, -6}}), 
      Rectangle(origin = {40, -34}, fillColor = {186, 189, 182}, fillPattern = FillPattern.VerticalCylinder, extent = {{-12, 24}, {12, -24}}), 
      Line(origin = {20, -15}, points = {{-8, 65}, {40, 65}, {40, 21}, {-40, 21}, {-40, -65}, {-22, -65}, {-20, -65}}, color = {32, 74, 135}, pattern = LinePattern.DashDotDot),
      Text(origin = {68, 143},lineColor = {0, 0, 255}, extent = {{-50, -61}, {51, -82}}, textString = "Smelter"),
      Text(origin = {82, 60},lineColor = {0, 0, 255}, extent = {{-40, -61}, {41, -82}}, textString = "BOF")
      }, 
      coordinateSystem(extent = {{-100, -100}, {100, 100}})),
      Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}})));

end SmelterBOF;

model GenericControl
	extends SolarTherm.Icons.Control;

	// Parameters
	parameter Real level_on = 10 "Fraction of Emax to start discharge";
	parameter Real level_off = 5 "Fraction of Emax to stop discharge";
	parameter Real level_curtailment_on = 99 "Fraction of Emax to start curtailment";
	parameter Real level_curtailment_off = 96 "Fraction of Emax to stop curtailment";

	parameter SI.Time t_process_ramp_up = 15*60 "Delay until power block starts";
	parameter SI.Time t_process_ramp_dw = 10*60 "Delay until power block shuts off";

	parameter Integer ramp_order = 1 "ramping filter order";

	parameter SI.HeatFlowRate Q_schedule_des = 5e7;

	// Connectors
	Modelica.Blocks.Interfaces.RealInput Q_schedule "Scheduled discharge heat flow rate" 
		annotation(
		Placement(visible = true, transformation(origin = {-108, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 0), 
		iconTransformation(origin = {-108, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
	
	Modelica.Blocks.Interfaces.RealInput storage_level "Instantaneous storage level" 
		annotation(
		Placement(visible = true, transformation(origin = {-108, 60}, extent = {{-20, -20}, {20, 20}}, rotation = 0), 
		iconTransformation(origin = {-108, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

	Modelica.Blocks.Interfaces.RealOutput Q_discharging "Heat flow rate discharging the storage" 
		annotation(
		Placement(visible = true, transformation(origin = {108, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), 
		iconTransformation(origin = {108, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
	
	Modelica.Blocks.Interfaces.BooleanOutput curtailment "Curtailment signal to power input" 
		annotation(Placement(transformation(origin = {0, -107}, extent = {{10, -10}, {-10, 10}}, rotation = 90)));

	SolarTherm.Utilities.Transition.Ramp ramp_up_proc(ramp_order=ramp_order, t_dur= t_process_ramp_up, up=true);
	SolarTherm.Utilities.Transition.Ramp ramp_down_proc(ramp_order=ramp_order, t_dur= t_process_ramp_dw, up=false);

	SI.Time  t_blk_w_now "Time of power block current warm-up event";
	SI.Time  t_blk_w_next "Time of power block next warm-up event";
	SI.Time  t_blk_c_now "Time of power block current cool-down event";
	SI.Time  t_blk_c_next "Time of power block next cool-down event";

	Integer process_state(min=1, max=4) "Power block state";
	Real fr_ramp_blk (min=0.0, max=1.0) "ramping transition rate for the power block";

algorithm
	when process_state == 2 and Q_schedule <= 0 then
		process_state := 1; // turn off (or stop ramping) due to no demand
	elsewhen process_state == 2 and storage_level <= level_off then
		process_state := 1; // turn off (or stop ramping) due to empty tank
	elsewhen process_state == 3 and Q_schedule <= 0 and t_process_ramp_dw > 0 then
		process_state := 4; // ramp down due to no demand
	elsewhen process_state == 3 and Q_schedule <= 0 and t_process_ramp_dw <= 0 then
		process_state := 1; // turn off (no ramp-down) due to no demand
	elsewhen process_state == 3 and storage_level <= level_off and t_process_ramp_dw > 0 then
		process_state := 4; // ramp down due to empty tank
	elsewhen process_state == 3 and storage_level <= level_off and t_process_ramp_dw <= 0 then
		process_state := 1; // turn off (no ramp down) due to empty tank
	elsewhen process_state == 2 and time >= t_blk_w_next then
		process_state := 3; // operational, ramp-up completed
	elsewhen process_state == 1 and Q_schedule > 0 and storage_level >= level_on  and t_process_ramp_up > 0 then
		process_state := 2; // ramp up, demand and tank has capacity
	elsewhen process_state == 1 and Q_schedule > 0 and storage_level >= level_on  and t_process_ramp_up <= 0 then
		process_state := 3; // operational (no ramp-up)
	elsewhen process_state == 4 and time >= t_blk_c_next then
		process_state := 1; // turn off after the ramp-down is complete
	end when;

	when process_state == 2 then
		t_blk_w_now := time;
		t_blk_w_next := time + t_process_ramp_up;
	end when;

	when process_state == 4 then
		t_blk_c_now := time;
		t_blk_c_next := time + t_process_ramp_dw;
	end when;

	when storage_level > level_curtailment_on then
		curtailment := true;
	elsewhen storage_level < level_curtailment_off then
		curtailment := false;
	end when;

	if process_state == 2 then
		fr_ramp_blk := if ramp_order == 0 then 0.0 else abs(ramp_up_proc.y);
	elseif process_state == 4 then
		fr_ramp_blk := if ramp_order == 0 then 0.0 else abs(ramp_down_proc.y);
	else
		fr_ramp_blk := 0;
	end if;

initial equation
	process_state = 3;
	t_blk_w_now = 0;
	t_blk_w_next = 0;
	t_blk_c_now = 0;
	t_blk_c_next = 0;

	if storage_level > level_curtailment_off then
		curtailment = true;
	elseif storage_level < level_curtailment_off then
		curtailment = false;
	else
		curtailment = true;
	end if;

equation

	ramp_up_proc.x = t_blk_w_now;
	ramp_down_proc.x = t_blk_c_now;

	if process_state <=1 then
		Q_discharging = 0;
	elseif process_state == 2 then
		Q_discharging = if ramp_order == 0 then Q_schedule else fr_ramp_blk * Q_schedule;
	elseif process_state == 4 then
		Q_discharging = fr_ramp_blk * Q_schedule;
	else
		Q_discharging = Q_schedule;
	end if;

annotation(
    Icon(graphics = {Text(origin = {-10, 254}, lineColor = {0, 0, 255}, extent = {{-149, -114}, {151, -154}}, textString = "%name")}));
end GenericControl;
model EnergyMerger "Merges input between two output ports"
	import SI = Modelica.SIunits;
	input Real frac(min=0, max=1) "Fraction to output port 1";
	SolarTherm.Models.Fluid.Interfaces.EnergyPort p_i;
	SolarTherm.Models.Fluid.Interfaces.EnergyPort p_o1;
	SolarTherm.Models.Fluid.Interfaces.EnergyPort p_o2;
 Modelica.Blocks.Interfaces.RealInput u1 annotation(
      Placement(visible = true, transformation(origin = {-100, 20}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-16, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
 Modelica.Blocks.Interfaces.RealInput u2 annotation(
      Placement(visible = true, transformation(origin = {-100, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-16, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
 Modelica.Blocks.Interfaces.RealOutput y annotation(
      Placement(visible = true, transformation(origin = {100, 1.9984e-15}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {16, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

equation
	p_o1.P = -frac*p_i.P;
	p_o2.P = -(1 - frac)*p_i.P;
annotation(
	Icon(graphics = {Rectangle(fillColor = {150, 150, 150}, fillPattern = FillPattern.Solid, extent = {{-5, 100}, {5, -100}})}, coordinateSystem(extent = {{-40, -100}, {40, 100}})),
      Diagram(coordinateSystem(extent = {{-40, -100}, {40, 100}})));

end EnergyMerger;
model EnergySplitter "Merges input between two output ports"
	import SI = Modelica.SIunits;
	input Real frac(min=0, max=1) "Fraction to output port 1";
	SolarTherm.Models.Fluid.Interfaces.EnergyPort p_i;
	SolarTherm.Models.Fluid.Interfaces.EnergyPort p_o1;
	SolarTherm.Models.Fluid.Interfaces.EnergyPort p_o2;
 Modelica.Blocks.Interfaces.RealInput u1 annotation(
      Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-16, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
 Modelica.Blocks.Interfaces.RealOutput y1 annotation(
      Placement(visible = true, transformation(origin = {100, 20}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {16, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
 Modelica.Blocks.Interfaces.RealOutput y2 annotation(
      Placement(visible = true, transformation(origin = {100, -20}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {16, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

equation
	p_o1.P = -frac*p_i.P;
	p_o2.P = -(1 - frac)*p_i.P;
annotation(
	Icon(graphics = {Rectangle(fillColor = {150, 150, 150}, fillPattern = FillPattern.Solid, extent = {{-5, 100}, {5, -100}})}, coordinateSystem(extent = {{-40, -100}, {40, 100}})),
      Diagram(coordinateSystem(extent = {{-40, -100}, {40, 100}})));

end EnergySplitter;

end Systems;