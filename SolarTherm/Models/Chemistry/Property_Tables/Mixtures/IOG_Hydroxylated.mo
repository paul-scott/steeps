within SolarTherm.Models.Chemistry.Property_Tables.Mixtures;

record IOG_Hydroxylated "51.29wt% Fe2O3.H2O, 41.98wt% Fe2O3, 2.62wt% Al2O3, 4.11wt% SiO2"
  constant Modelica.SIunits.Temperature T_table[122] = {298.15, 300.00, 310.00, 320.00, 330.00, 340.00, 350.00, 360.00, 370.00, 380.00, 390.00, 400.00, 410.00, 420.00, 430.00, 440.00, 450.00, 460.00, 470.00, 480.00, 490.00, 500.00, 510.00, 520.00, 530.00, 540.00, 550.00, 560.00, 570.00, 580.00, 590.00, 600.00, 610.00, 620.00, 630.00, 640.00, 650.00, 660.00, 670.00, 680.00, 690.00, 700.00, 710.00, 720.00, 730.00, 740.00, 750.00, 760.00, 770.00, 780.00, 790.00, 800.00, 810.00, 820.00, 830.00, 840.00, 850.00, 860.00, 870.00, 880.00, 890.00, 900.00, 910.00, 920.00, 930.00, 940.00, 950.00, 960.00, 970.00, 980.00, 990.00, 1000.00, 1010.00, 1020.00, 1030.00, 1040.00, 1050.00, 1060.00, 1070.00, 1080.00, 1090.00, 1100.00, 1110.00, 1120.00, 1130.00, 1140.00, 1150.00, 1160.00, 1170.00, 1180.00, 1190.00, 1200.00, 1210.00, 1220.00, 1230.00, 1240.00, 1250.00, 1260.00, 1270.00, 1280.00, 1290.00, 1300.00, 1310.00, 1320.00, 1330.00, 1340.00, 1350.00, 1360.00, 1370.00, 1380.00, 1390.00, 1400.00, 1410.00, 1420.00, 1430.00, 1440.00, 1450.00, 1460.00, 1470.00, 1480.00, 1490.00, 1500.00} "Absolute temperature (K)";
  constant Modelica.SIunits.SpecificEnthalpy h_table[122] = {0.0, 1548.4, 10024.6, 18674.1, 27492.1, 36474.6, 45616.5, 54915.9, 64368.9, 73972.0, 83724.1, 93622.3, 103611.6, 113636.7, 123695.0, 133786.9, 143910.9, 154064.8, 164248.6, 174461.3, 184701.7, 194969.8, 205263.9, 215585.1, 225932.1, 236305.0, 246703.7, 257127.1, 267576.0, 278050.0, 288548.5, 299072.9, 309622.7, 320198.1, 330798.9, 341426.1, 352079.5, 362758.9, 373466.2, 384201.3, 394964.6, 405756.6, 416578.6, 427430.6, 438313.5, 449228.6, 460176.8, 471159.5, 482177.6, 493232.3, 504325.0, 515457.5, 526630.6, 537847.8, 549110.2, 560419.8, 572105.2, 583435.6, 594817.0, 606252.9, 617746.0, 629301.1, 640922.3, 652613.8, 664380.9, 676230.0, 688166.8, 699780.9, 710835.2, 721840.5, 732804.2, 743731.7, 754627.6, 765496.5, 776341.3, 787164.4, 797968.9, 808756.8, 819529.9, 830289.2, 841036.0, 851771.7, 862497.3, 873212.9, 883920.7, 894620.9, 906740.3, 917432.7, 928118.7, 938799.2, 949474.1, 960144.3, 970809.8, 981471.1, 992128.8, 1002782.1, 1013432.4, 1024079.5, 1034723.1, 1045364.0, 1056002.6, 1066638.8, 1077272.4, 1087904.0, 1098533.6, 1109161.3, 1119786.6, 1130410.9, 1141033.6, 1151654.7, 1162274.7, 1172893.0, 1183510.7, 1194126.2, 1204741.5, 1215355.7, 1225969.4, 1236581.5, 1247193.3, 1257804.6, 1268414.8, 1279024.0} "Specific enthalpy (J/kg)";
  constant Modelica.SIunits.SpecificHeatCapacityAtConstantPressure cp_table[122] = {835.41, 838.74, 856.36, 873.46, 890.08, 906.29, 922.13, 937.65, 952.86, 967.82, 982.54, 997.06, 1000.72, 1004.22, 1007.57, 1010.80, 1013.90, 1016.91, 1019.82, 1022.66, 1025.43, 1028.14, 1030.80, 1033.42, 1036.00, 1038.56, 1041.10, 1043.62, 1046.14, 1048.65, 1051.17, 1053.70, 1056.24, 1058.80, 1061.39, 1064.01, 1066.66, 1069.36, 1072.10, 1074.90, 1077.76, 1080.68, 1083.67, 1086.75, 1089.91, 1093.16, 1096.53, 1100.01, 1103.61, 1107.34, 1111.23, 1115.28, 1119.50, 1123.93, 1128.56, 1133.43, 1130.56, 1135.54, 1140.81, 1146.41, 1152.37, 1158.73, 1165.54, 1172.85, 1180.70, 1189.17, 1198.31, 1108.23, 1102.89, 1098.35, 1094.46, 1091.10, 1088.18, 1085.62, 1083.37, 1081.37, 1079.60, 1078.00, 1076.57, 1075.28, 1074.11, 1073.05, 1072.07, 1071.18, 1070.37, 1069.62, 1069.56, 1068.92, 1068.33, 1067.77, 1067.26, 1066.78, 1066.34, 1065.93, 1065.55, 1065.19, 1064.85, 1064.54, 1064.25, 1063.97, 1063.72, 1063.48, 1063.26, 1063.05, 1062.85, 1062.67, 1062.50, 1062.34, 1062.19, 1062.05, 1061.91, 1061.79, 1061.68, 1061.57, 1061.47, 1061.38, 1061.29, 1061.21, 1061.13, 1061.06, 1060.99, 1060.93} "Specific heat capacity at constant pressure (J/kgK)";
  constant Modelica.SIunits.SpecificEntropy s_table[122] = {615.01, 620.19, 647.98, 675.44, 702.57, 729.39, 755.89, 782.08, 807.98, 833.59, 858.92, 883.98, 908.64, 932.80, 956.47, 979.67, 1002.42, 1024.74, 1046.64, 1068.14, 1089.26, 1110.00, 1130.39, 1150.43, 1170.14, 1189.53, 1208.60, 1227.39, 1245.88, 1264.10, 1282.05, 1299.74, 1317.18, 1334.37, 1351.33, 1368.06, 1384.58, 1400.89, 1416.99, 1432.89, 1448.60, 1464.14, 1479.48, 1494.66, 1509.67, 1524.53, 1539.22, 1553.77, 1568.17, 1582.43, 1596.57, 1610.57, 1624.45, 1638.21, 1651.86, 1665.41, 1679.23, 1692.49, 1705.65, 1718.72, 1731.70, 1744.61, 1757.45, 1770.23, 1782.95, 1795.62, 1808.26, 1820.42, 1831.87, 1843.16, 1854.29, 1865.27, 1876.12, 1886.83, 1897.40, 1907.87, 1918.20, 1928.43, 1938.54, 1948.55, 1958.46, 1968.26, 1977.97, 1987.58, 1997.10, 2006.52, 2017.11, 2026.37, 2035.54, 2044.64, 2053.64, 2062.57, 2071.42, 2080.20, 2088.90, 2097.52, 2106.08, 2114.57, 2122.98, 2131.33, 2139.60, 2147.81, 2155.96, 2164.05, 2172.07, 2180.03, 2187.93, 2195.77, 2203.56, 2211.28, 2218.94, 2226.56, 2234.12, 2241.62, 2249.07, 2256.46, 2263.81, 2271.10, 2278.35, 2285.54, 2292.68, 2299.78} "Specific entropy (J/kgK)";
  constant Modelica.SIunits.Density rho_table[122] = {4488.64, 4488.56, 4488.11, 4487.65, 4487.19, 4486.72, 4486.25, 4485.77, 4485.29, 4484.81, 4484.32, 4483.83, 4483.33, 4482.82, 4482.32, 4481.81, 4481.29, 4480.77, 4480.25, 4479.72, 4479.19, 4478.65, 4478.12, 4477.57, 4477.03, 4476.47, 4475.92, 4475.36, 4474.80, 4474.24, 4473.67, 4473.10, 4472.52, 4471.94, 4471.36, 4470.77, 4470.18, 4469.59, 4468.99, 4468.40, 4467.79, 4467.19, 4466.58, 4465.97, 4465.35, 4464.73, 4464.12, 4463.49, 4462.86, 4462.23, 4461.60, 4460.96, 4460.32, 4459.68, 4459.03, 4458.39, 4457.74, 4457.09, 4456.43, 4455.77, 4455.11, 4454.44, 4453.77, 4453.10, 4452.43, 4451.76, 4451.08, 4450.40, 4449.73, 4449.07, 4448.41, 4447.75, 4447.11, 4446.46, 4445.83, 4445.19, 4444.57, 4443.94, 4443.32, 4442.71, 4442.10, 4441.49, 4440.89, 4440.28, 4439.69, 4439.09, 4438.50, 4437.91, 4437.33, 4436.75, 4436.17, 4435.59, 4435.02, 4434.45, 4433.88, 4433.31, 4432.74, 4432.18, 4431.62, 4431.07, 4430.51, 4429.95, 4429.40, 4428.85, 4428.30, 4427.75, 4427.21, 4426.66, 4426.12, 4425.58, 4425.04, 4424.50, 4423.96, 4423.42, 4422.89, 4422.35, 4421.82, 4421.29, 4420.76, 4420.23, 4419.70, 4419.18} "Density (kg/m3)";
  constant Modelica.SIunits.ThermalConductivity k_table[122] = {12.63, 12.55, 12.15, 11.77, 11.41, 11.08, 10.77, 10.47, 10.19, 9.92, 9.67, 9.44, 9.20, 8.99, 8.78, 8.58, 8.39, 8.21, 8.04, 7.87, 7.71, 7.56, 7.41, 7.26, 7.13, 6.99, 6.86, 6.74, 6.62, 6.50, 6.38, 6.27, 6.16, 6.05, 5.94, 5.84, 5.74, 5.64, 5.54, 5.47, 5.41, 5.35, 5.29, 5.24, 5.18, 5.13, 5.08, 5.03, 4.98, 4.92, 4.84, 4.77, 4.70, 4.63, 4.57, 4.50, 4.44, 4.37, 4.31, 4.26, 4.21, 4.17, 4.12, 4.08, 4.04, 3.99, 3.95, 3.94, 3.93, 3.92, 3.90, 3.89, 3.88, 3.87, 3.86, 3.85, 3.84, 3.82, 3.82, 3.81, 3.80, 3.79, 3.78, 3.77, 3.77, 3.76, 3.75, 3.74, 3.73, 3.73, 3.72, 3.71, 3.70, 3.70, 3.69, 3.68, 3.67, 3.67, 3.66, 3.65, 3.64, 3.63, 3.63, 3.62, 3.61, 3.61, 3.60, 3.59, 3.59, 3.58, 3.57, 3.56, 3.55, 3.55, 3.54, 3.53, 3.53, 3.52, 3.52, 3.51, 3.50, 3.49} "Thermal conductivity (W/mK)";
  annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)), Documentation(info="<html><img src=\"modelica://SolarTherm/Resources/Properties_IOG_Hydroxylated.png\"></html>"));

end IOG_Hydroxylated;