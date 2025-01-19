within SolarTherm.Materials;
package Dowtherm_Table
  extends SolarTherm.Materials.PartialMaterial(MM = 166.0e-3, T_melt = 12.0 + 273.15, cost = 2.10); //Cost from Heller et al. 2013
  import SolarTherm.Utilities.Interpolation.Interpolate1D;
  
  //constant SI.MolarMass MM = 22.989769e-3 "Molar mass (kg/mol)";
  //constant SI.Temperature T_melt = 97.79 + 273.15 "Melting point (K)";
  //constant Real cost = 3.0 "USD/kg";
  constant SI.SpecificEnthalpy h_data[9] = {0.0, 65180.0, 135480.0, 229680.0, 330855.0, 438955.0, 554055.0, 676555.0, 807855.0};
  constant SI.Temperature T_data[9] = {298.15, 338.15, 378.15, 428.15, 478.15, 528.15, 578.15, 628.15, 678.15};
  constant SI.SpecificHeatCapacity cp_data[9] = {1586.6, 1701.0, 1814.0, 1954.0, 2093.0, 2231.0, 2373.0, 2527.0, 2725.0};
  constant SI.DynamicViscosity mu_data[9] = {0.00500, 0.00158, 0.00091, 0.00056, 0.00038, 0.00027, 0.00020, 0.00016, 0.00012};
  constant SI.ThermalConductivity k_data[9] = {0.1395, 0.1315, 0.1251, 0.1171, 0.1091, 0.1011, 0.0931, 0.0851, 0.0771};
  constant SI.Density rho_data[9] = {1063.5, 1023.7, 990.7, 947.8, 902.5, 854.0, 801.3, 724.3, 672.5};
  
  redeclare model State
    parameter String table_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/MaterialTables/DowthermA.txt");
    parameter SI.SpecificEnthalpy h_start=0.0;
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
end Dowtherm_Table;