within SolarTherm.Materials;
package Air_amb_p_curvefit
  extends SolarTherm.Materials.PartialMaterial(MM = 28.97e-3, T_melt = 58, cost = 0.0);
  import SolarTherm.Utilities.Interpolation.Interpolate1D;
  
  //constant SI.MolarMass MM = 28.97e-3 "Molar mass (kg/mol)";
  //constant SI.Temperature T_melt = 0.0 "Melting point (K)";
  //constant Real cost = 0.0 "USD/kg";
  
  redeclare model State "A model which calculates state and properties"
    parameter SI.SpecificEnthalpy h_start;
    //parameter String table_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/MaterialTables/Air_Table.txt");
    SI.Density rho;
    SI.SpecificEnthalpy h (start=h_start);
    SI.Temperature T;
    SI.DynamicViscosity mu;
    SI.SpecificHeatCapacity cp;
    SI.ThermalConductivity k;
    	//Tables.CombiTable1Ds Tab (tableOnFile=true, tableName="table_1D_1", columns=2:9, fileName=table_file);
  equation
    T = SolarTherm.Media.Air.Air_amb_p_curvefit_utilities.T_h(h);
    rho = SolarTherm.Media.Air.Air_amb_p_curvefit_utilities.rho_T(T);
    cp = SolarTherm.Media.Air.Air_amb_p_curvefit_utilities.cp_T(T);
    mu = SolarTherm.Media.Air.Air_amb_p_curvefit_utilities.mu_T(T);
    k = SolarTherm.Media.Air.Air_amb_p_curvefit_utilities.k_T(T);
	//Tab.u = h;
	//T = Tab.y[1];
	//rho = Tab.y[2];
	//cp = Tab.y[3];
	//mu = Tab.y[4];
	//k = Tab.y[6];
  end State;
  
  redeclare function h_Tf "Specific enthalpy of air vs Temperature"
    input SI.Temperature T "Absolute temperature (K)";
    input Real f = 0 "Liquid mass melt fraction (No effect on result)";
    output SI.SpecificEnthalpy h "Specific Enthalpy (J/kg)";
  algorithm
    h := SolarTherm.Media.Air.Air_amb_p_curvefit_utilities.h_T(T);
  end h_Tf;
  
  redeclare function rho_Tf "Density of air vs Temperature"
    input SI.Temperature T;
    input Real f = 0 "Liquid mass melt fraction (No effect on result)";
    output SI.Density rho;
  algorithm
    rho := SolarTherm.Media.Air.Air_amb_p_curvefit_utilities.rho_T(T);
  end rho_Tf;


  function T_h "Temperature interpolated from enthalpy"
    input SI.SpecificEnthalpy h;
    output SI.Temperature T;
  algorithm
    T := SolarTherm.Media.Air.Air_amb_p_curvefit_utilities.T_h(h);
  end T_h;
end Air_amb_p_curvefit;