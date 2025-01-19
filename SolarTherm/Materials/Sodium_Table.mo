within SolarTherm.Materials;
package Sodium_Table
  extends SolarTherm.Materials.PartialMaterial(MM = 22.989769e-3, T_melt = 97.79 + 273.15, cost = 3.0);
  import SolarTherm.Utilities.Interpolation.Interpolate1D;
  
  //constant SI.MolarMass MM = 22.989769e-3 "Molar mass (kg/mol)";
  //constant SI.Temperature T_melt = 97.79 + 273.15 "Melting point (K)";
  //constant Real cost = 3.0 "USD/kg";
  constant SI.SpecificEnthalpy h_data[41] = {503167.0, 515982.3, 528775.6, 541547.6, 554299.2, 567031.2, 579744.3, 592439.6, 605117.7, 617779.5, 630425.8, 643057.5, 655675.4, 668280.4, 680873.4, 693455.0, 706026.3, 718588.1, 731141.1, 743686.4, 756224.6, 768756.8, 781283.7, 793806.2, 806325.2, 818841.6, 831356.1, 843869.8, 856383.3, 868897.7, 881413.8, 893932.4, 906454.4, 918980.7, 931512.2, 944049.7, 956594.2, 969146.4, 981707.3, 994277.7, 1006858.6};
  constant SI.Temperature T_data[41] = {673.15, 683.15, 693.15, 703.15, 713.15, 723.15, 733.15, 743.15, 753.15, 763.15, 773.15, 783.15, 793.15, 803.15, 813.15, 823.15, 833.15, 843.15, 853.15, 863.15, 873.15, 883.15, 893.15, 903.15, 913.15, 923.15, 933.15, 943.15, 953.15, 963.15, 973.15, 983.15, 993.15, 1003.15, 1013.15, 1023.15, 1033.15, 1043.15, 1053.15, 1063.15, 1073.15};
  constant SI.SpecificHeatCapacity cp_data[41] = {1282.661, 1280.415, 1278.250, 1276.166, 1274.164, 1272.244, 1270.406, 1268.652, 1266.980, 1265.393, 1263.888, 1262.468, 1261.133, 1259.881, 1258.715, 1257.634, 1256.637, 1255.726, 1254.901, 1254.161, 1253.507, 1252.939, 1252.458, 1252.062, 1251.753, 1251.531, 1251.395, 1251.345, 1251.383, 1251.507, 1251.719, 1252.017, 1252.403, 1252.876, 1253.436, 1254.084, 1254.819, 1255.641, 1256.551, 1257.549, 1258.634};
  constant SI.DynamicViscosity mu_data[41] = {0.0002771840, 0.0002722542, 0.0002675282, 0.0002629939, 0.0002586400, 0.0002544561, 0.0002504326, 0.0002465605, 0.0002428316, 0.0002392381, 0.0002357730, 0.0002324294, 0.0002292011, 0.0002260824, 0.0002230678, 0.0002201521, 0.0002173307, 0.0002145989, 0.0002119527, 0.0002093881, 0.0002069013, 0.0002044888, 0.0002021475, 0.0001998741, 0.0001976657, 0.0001955196, 0.0001934332, 0.0001914040, 0.0001894296, 0.0001875079, 0.0001856368, 0.0001838142, 0.0001820383, 0.0001803073, 0.0001786195, 0.0001769733, 0.0001753671, 0.0001737994, 0.0001722690, 0.0001707744, 0.0001693144};
  constant SI.ThermalConductivity k_data[41] = {69.47130, 68.91885, 68.37259, 67.83245, 67.29835, 66.77024, 66.24803, 65.73166, 65.22106, 64.71615, 64.21686, 63.72312, 63.23487, 62.75202, 62.27451, 61.80227, 61.33523, 60.87331, 60.41645, 59.96457, 59.51760, 59.07547, 58.63811, 58.20546, 57.77743, 57.35396, 56.93497, 56.52040, 56.11017, 55.70422, 55.30247, 54.90485, 54.51128, 54.12171, 53.73606, 53.35425, 52.97621, 52.60188, 52.23119, 51.86406, 51.50041};
  constant SI.Density rho_data[41] = {857.7316, 855.4355, 853.1361, 850.8334, 848.5273, 846.2179, 843.9050, 841.5888, 839.2690, 836.9458, 834.6191, 832.2888, 829.9549, 827.6174, 825.2763, 822.9315, 820.5830, 818.2307, 815.8747, 813.5149, 811.1512, 808.7836, 806.4121, 804.0366, 801.6572, 799.2737, 796.8862, 794.4946, 792.0988, 789.6988, 787.2946, 784.8861, 782.4733, 780.0562, 777.6347, 775.2087, 772.7782, 770.3432, 767.9037, 765.4595, 763.0106};
  
  redeclare model State
    parameter String table_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/MaterialTables/Sodium.txt");
    parameter SI.SpecificEnthalpy h_start=700000;
    SI.Density rho;
    SI.SpecificEnthalpy h (start=h_start);
    SI.Temperature T;
    SI.DynamicViscosity mu;
    SI.SpecificHeatCapacity cp;
    SI.ThermalConductivity k;
    Tables.CombiTable1Ds Tab (tableOnFile=true, tableName="table_1D_1", columns=2:6, fileName=table_file);
  equation
    	Tab.u = h;
	T = Tab.y[1];
	cp = Tab.y[2];
    mu = Tab.y[3];
    k = Tab.y[4];
    rho = Tab.y[5];
  end State;

  redeclare function h_Tf
    input SI.Temperature T;
    input Real f = 0 "Liquid mass melt fraction (No effect on result)";
    output SI.SpecificEnthalpy h;
  algorithm
    h := SolarTherm.Utilities.Interpolation.Interpolate1D(T_data,h_data,T);
  end h_Tf;
  
  redeclare function rho_Tf
    input SI.Temperature T;
    input Real f = 0 "Liquid mass melt fraction (No effect on result)";
    output SI.Density rho;
  algorithm
    rho := SolarTherm.Utilities.Interpolation.Interpolate1D(T_data,rho_data,T);
  end rho_Tf;
  
  function mu_T
    input SI.Temperature T;
    output SI.DynamicViscosity mu;
  algorithm
    mu := SolarTherm.Utilities.Interpolation.Interpolate1D(T_data,mu_data,T);
  end mu_T;
end Sodium_Table;