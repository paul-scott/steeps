within SolarTherm.Media.SolidParticles.FeO_utilities;
function k_T "Thermal conductivity of Fe3O4 as a function of temperature"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Temperature T "Temperature (K)";
	output Modelica.SIunits.ThermalConductivity k "Thermal conductivity (W/mK)";
protected
    constant Modelica.SIunits.Temperature T_data[122] = {298.15, 300.00, 310.00, 320.00, 330.00, 340.00, 350.00, 360.00, 370.00, 380.00, 390.00, 400.00, 410.00, 420.00, 430.00, 440.00, 450.00, 460.00, 470.00, 480.00, 490.00, 500.00, 510.00, 520.00, 530.00, 540.00, 550.00, 560.00, 570.00, 580.00, 590.00, 600.00, 610.00, 620.00, 630.00, 640.00, 650.00, 660.00, 670.00, 680.00, 690.00, 700.00, 710.00, 720.00, 730.00, 740.00, 750.00, 760.00, 770.00, 780.00, 790.00, 800.00, 810.00, 820.00, 830.00, 840.00, 850.00, 860.00, 870.00, 880.00, 890.00, 900.00, 910.00, 920.00, 930.00, 940.00, 950.00, 960.00, 970.00, 980.00, 990.00, 1000.00, 1010.00, 1020.00, 1030.00, 1040.00, 1050.00, 1060.00, 1070.00, 1080.00, 1090.00, 1100.00, 1110.00, 1120.00, 1130.00, 1140.00, 1150.00, 1160.00, 1170.00, 1180.00, 1190.00, 1200.00, 1210.00, 1220.00, 1230.00, 1240.00, 1250.00, 1260.00, 1270.00, 1280.00, 1290.00, 1300.00, 1310.00, 1320.00, 1330.00, 1340.00, 1350.00, 1360.00, 1370.00, 1380.00, 1390.00, 1400.00, 1410.00, 1420.00, 1430.00, 1440.00, 1450.00, 1460.00, 1470.00, 1480.00, 1490.00, 1500.00};
    constant Modelica.SIunits.ThermalConductivity k_data[122] = {5.282, 5.268, 5.194, 5.123, 5.056, 4.991, 4.929, 4.869, 4.811, 4.755, 4.702, 4.650, 4.600, 4.552, 4.505, 4.460, 4.416, 4.373, 4.332, 4.292, 4.253, 4.215, 4.178, 4.142, 4.107, 4.073, 4.040, 4.008, 3.976, 3.946, 3.916, 3.886, 3.858, 3.830, 3.802, 3.775, 3.749, 3.724, 3.699, 3.674, 3.650, 3.626, 3.603, 3.580, 3.558, 3.536, 3.515, 3.494, 3.473, 3.453, 3.433, 3.414, 3.394, 3.376, 3.357, 3.339, 3.340, 3.348, 3.357, 3.366, 3.375, 3.384, 3.393, 3.403, 3.412, 3.421, 3.430, 3.440, 3.449, 3.459, 3.468, 3.478, 3.487, 3.497, 3.507, 3.516, 3.526, 3.536, 3.546, 3.556, 3.566, 3.576, 3.586, 3.597, 3.607, 3.617, 3.628, 3.638, 3.649, 3.659, 3.670, 3.680, 3.691, 3.702, 3.713, 3.724, 3.735, 3.746, 3.757, 3.768, 3.780, 3.791, 3.802, 3.814, 3.826, 3.837, 3.849, 3.861, 3.873, 3.884, 3.896, 3.909, 3.921, 3.933, 3.945, 3.958, 3.970, 3.983, 3.995, 4.008, 4.021, 4.033};
algorithm
	k := Utilities.Interpolation.Interpolate1D(T_data,k_data,T);
end k_T;