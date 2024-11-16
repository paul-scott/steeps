within SolarTherm.Models.Chemistry.Property_Tables;

record Fe2O3_H2O
  constant Modelica.SIunits.Temperature T_table[122] = {298.15, 300.00, 310.00, 320.00, 330.00, 340.00, 350.00, 360.00, 370.00, 380.00, 390.00, 400.00, 410.00, 420.00, 430.00, 440.00, 450.00, 460.00, 470.00, 480.00, 490.00, 500.00, 510.00, 520.00, 530.00, 540.00, 550.00, 560.00, 570.00, 580.00, 590.00, 600.00, 610.00, 620.00, 630.00, 640.00, 650.00, 660.00, 670.00, 680.00, 690.00, 700.00, 710.00, 720.00, 730.00, 740.00, 750.00, 760.00, 770.00, 780.00, 790.00, 800.00, 810.00, 820.00, 830.00, 840.00, 850.00, 860.00, 870.00, 880.00, 890.00, 900.00, 910.00, 920.00, 930.00, 940.00, 950.00, 960.00, 970.00, 980.00, 990.00, 1000.00, 1010.00, 1020.00, 1030.00, 1040.00, 1050.00, 1060.00, 1070.00, 1080.00, 1090.00, 1100.00, 1110.00, 1120.00, 1130.00, 1140.00, 1150.00, 1160.00, 1170.00, 1180.00, 1190.00, 1200.00, 1210.00, 1220.00, 1230.00, 1240.00, 1250.00, 1260.00, 1270.00, 1280.00, 1290.00, 1300.00, 1310.00, 1320.00, 1330.00, 1340.00, 1350.00, 1360.00, 1370.00, 1380.00, 1390.00, 1400.00, 1410.00, 1420.00, 1430.00, 1440.00, 1450.00, 1460.00, 1470.00, 1480.00, 1490.00, 1500.00} "Absolute temperature (K)";
  constant Modelica.SIunits.SpecificEnthalpy h_table[122] = {0.00, 1833.00, 11864.00, 22104.00, 32551.00, 43207.00, 54070.00, 65142.00, 76422.00, 87909.00, 99605.00, 111509.00, 123517.00, 135525.00, 147532.00, 159540.00, 171548.00, 183556.00, 195564.00, 207572.00, 219580.00, 231588.00, 243595.00, 255603.00, 267611.00, 279619.00, 291627.00, 303635.00, 315643.00, 327651.00, 339658.00, 351666.00, 363674.00, 375682.00, 387690.00, 399698.00, 411706.00, 423713.00, 435721.00, 447729.00, 459737.00, 471745.00, 483753.00, 495761.00, 507769.00, 519776.00, 531784.00, 543792.00, 555800.00, 567808.00, 579816.00, 591824.00, 603831.00, 615839.00, 627847.00, 639855.00, 651863.00, 663871.00, 675879.00, 687887.00, 699894.00, 711902.00, 723910.00, 735918.00, 747926.00, 759934.00, 771942.00, 783950.00, 795957.00, 807965.00, 819973.00, 831981.00, 843989.00, 855997.00, 868005.00, 880012.00, 892020.00, 904028.00, 916036.00, 928044.00, 940052.00, 952060.00, 964068.00, 976075.00, 988083.00, 1000091.00, 1012099.00, 1024107.00, 1036115.00, 1048123.00, 1060130.00, 1072138.00, 1084146.00, 1096154.00, 1108162.00, 1120170.00, 1132178.00, 1144186.00, 1156193.00, 1168201.00, 1180209.00, 1192217.00, 1204225.00, 1216233.00, 1228241.00, 1240249.00, 1252256.00, 1264264.00, 1276272.00, 1288280.00, 1300288.00, 1312296.00, 1324304.00, 1336311.00, 1348319.00, 1360327.00, 1372335.00, 1384343.00, 1396351.00, 1408359.00, 1420367.00, 1432374.00} "Specific enthalpy (J/kg)";
  constant Modelica.SIunits.SpecificHeatCapacityAtConstantPressure cp_table[122] = {988.88, 992.73, 1013.54, 1034.34, 1055.15, 1075.95, 1096.76, 1117.57, 1138.37, 1159.18, 1179.98, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79, 1200.79} "Specific heat capacity at constant pressure (J/kgK)";
  constant Modelica.SIunits.SpecificEntropy s_table[122] = {668.68, 674.80, 707.70, 740.20, 772.35, 804.16, 835.65, 866.84, 897.74, 928.37, 958.75, 988.89, 1018.54, 1047.48, 1075.73, 1103.34, 1130.32, 1156.71, 1182.54, 1207.82, 1232.58, 1256.84, 1280.62, 1303.93, 1326.81, 1349.25, 1371.28, 1392.92, 1414.17, 1435.06, 1455.59, 1475.77, 1495.62, 1515.14, 1534.35, 1553.26, 1571.88, 1590.21, 1608.27, 1626.06, 1643.59, 1660.87, 1677.90, 1694.70, 1711.26, 1727.60, 1743.71, 1759.62, 1775.32, 1790.81, 1806.11, 1821.21, 1836.13, 1850.86, 1865.42, 1879.80, 1894.01, 1908.05, 1921.94, 1935.66, 1949.23, 1962.64, 1975.91, 1989.04, 2002.02, 2014.86, 2027.57, 2040.14, 2052.58, 2064.90, 2077.09, 2089.16, 2101.11, 2112.94, 2124.65, 2136.26, 2147.75, 2159.13, 2170.40, 2181.57, 2192.64, 2203.61, 2214.47, 2225.24, 2235.92, 2246.50, 2256.98, 2267.38, 2277.69, 2287.91, 2298.04, 2308.09, 2318.05, 2327.94, 2337.74, 2347.46, 2357.11, 2366.68, 2376.17, 2385.59, 2394.93, 2404.20, 2413.40, 2422.54, 2431.60, 2440.59, 2449.52, 2458.38, 2467.18, 2475.91, 2484.58, 2493.19, 2501.74, 2510.22, 2518.65, 2527.02, 2535.33, 2543.58, 2551.78, 2559.92, 2568.00, 2576.04} "Specific entropy (J/kgK)";
  constant Modelica.SIunits.Density rho_table[122] = fill(4248.54,122) "Density (kg/m3)";
  constant Modelica.SIunits.ThermalConductivity k_table[122] = {13.13, 13.04, 12.61, 12.20, 11.81, 11.45, 11.11, 10.79, 10.49, 10.20, 9.93, 9.68, 9.43, 9.20, 8.98, 8.77, 8.57, 8.38, 8.20, 8.02, 7.85, 7.69, 7.54, 7.39, 7.25, 7.11, 6.98, 6.85, 6.73, 6.61, 6.49, 6.38, 6.28, 6.17, 6.07, 5.98, 5.88, 5.79, 5.70, 5.62, 5.54, 5.45, 5.38, 5.30, 5.23, 5.15, 5.08, 5.02, 4.95, 4.89, 4.82, 4.76, 4.70, 4.64, 4.59, 4.53, 4.48, 4.42, 4.37, 4.32, 4.27, 4.23, 4.18, 4.13, 4.09, 4.04, 4.00, 3.99, 3.98, 3.97, 3.96, 3.95, 3.94, 3.93, 3.92, 3.91, 3.90, 3.89, 3.89, 3.88, 3.87, 3.86, 3.85, 3.84, 3.83, 3.82, 3.81, 3.80, 3.79, 3.79, 3.78, 3.77, 3.76, 3.75, 3.74, 3.73, 3.72, 3.72, 3.71, 3.70, 3.69, 3.68, 3.67, 3.67, 3.66, 3.65, 3.64, 3.63, 3.63, 3.62, 3.61, 3.60, 3.59, 3.59, 3.58, 3.57, 3.56, 3.55, 3.55, 3.54, 3.53, 3.52} "Thermal conductivity (W/mK)";
  annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)), Documentation(info="<html><img src=\"modelica://SolarTherm/Resources/Properties_Fe2O3H2O.png\"></html>"));
end Fe2O3_H2O;