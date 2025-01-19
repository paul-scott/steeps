within SolarTherm.Materials;
package Air_Table
  extends SolarTherm.Materials.PartialMaterial(MM = 28.97e-3, T_melt = 58, cost = 0.0);
  import SolarTherm.Utilities.Interpolation.Interpolate1D;
  
  //constant SI.MolarMass MM = 28.97e-3 "Molar mass (kg/mol)";
  //constant SI.Temperature T_melt = 0.0 "Melting point (K)";
  //constant Real cost = 0.0 "USD/kg";
  
  redeclare model State "A model which calculates state and properties"
    parameter SI.SpecificEnthalpy h_start;
    parameter String table_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/MaterialTables/Air_Table.txt");
    SI.Density rho;
    SI.SpecificEnthalpy h (start=h_start);
    SI.Temperature T;
    SI.DynamicViscosity mu;
    SI.SpecificHeatCapacity cp;
    SI.ThermalConductivity k;
    	Tables.CombiTable1Ds Tab (tableOnFile=true, tableName="table_1D_1", columns=2:9, fileName=table_file);
  equation
	Tab.u = h;
	T = Tab.y[1];
	rho = Tab.y[2];
	cp = Tab.y[3];
	mu = Tab.y[4];
	k = Tab.y[6];
  end State;
  
  redeclare function h_Tf "Specific enthalpy of air vs Temperature"
    input SI.Temperature T "Absolute temperature (K)";
    input Real f = 0 "Liquid mass melt fraction (No effect on result)";
    output SI.SpecificEnthalpy h "Specific Enthalpy (J/kg)";
  protected
    Real T_data[35] = {100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,3000};
    Real h_data[35] = {0,51100,101575,151900,202225,252625,303200,354075,405350,457100,509375,562225,615675,669725,724375,779600,835375,891675,962737.5,1063475,1180175,1298375,1418175,1540025,1663925,1789675,1917325,2046975,2179175,2314625,2454075,2598825,2750625,3078275,4009525};
  algorithm
    h := SolarTherm.Utilities.Interpolation.Interpolate1D(T_data,h_data,T);
  end h_Tf;
  
  redeclare function rho_Tf "Density of air vs Temperature"
    input SI.Temperature T;
    input Real f = 0 "Liquid mass melt fraction (No effect on result)";
    output SI.Density rho;
  protected
    Real T_data[35] = {100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,3000};
    Real rho_data[35] = {3.5562,2.3364,1.7458,1.3947,1.1614,0.995,0.8711,0.774,0.6964,0.6329,0.5804,0.5356,0.4975,0.4643,0.4354,0.4097,0.3868,0.3666,0.3482,0.3166,0.2902,0.2679,0.2488,0.2322,0.2177,0.2049,0.1935,0.1833,0.1741,0.1658,0.1582,0.1513,0.1448,0.1389,0.1135};
  algorithm
    rho := SolarTherm.Utilities.Interpolation.Interpolate1D(T_data,rho_data,T);
  end rho_Tf;


  function T_h "Temperature interpolated from enthalpy"
    input SI.SpecificEnthalpy h;
    output SI.Temperature T;
  protected 
    Real T_data[35] = {100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,3000};
    Real h_data[35] = {0,51100,101575,151900,202225,252625,303200,354075,405350,457100,509375,562225,615675,669725,724375,779600,835375,891675,962737.5,1063475,1180175,1298375,1418175,1540025,1663925,1789675,1917325,2046975,2179175,2314625,2454075,2598825,2750625,3078275,4009525};
  algorithm
    T := SolarTherm.Utilities.Interpolation.Interpolate1D(h_data,T_data,h);
  end T_h;
end Air_Table;