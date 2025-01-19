within SolarTherm.Validation.Gan_HBS.Materials.Flue_Gas_5884Pa;

model State "A model which calculates state and properties"
  parameter SI.SpecificEnthalpy h_start;
  parameter String table_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/MaterialTables/Air_Table.txt");
  SI.Density rho;
  SI.SpecificEnthalpy h(start = h_start);
  SI.Temperature T;
  SI.DynamicViscosity mu;
  SI.SpecificHeatCapacity cp;
  SI.ThermalConductivity k;
  Tables.CombiTable1Ds Tab(tableOnFile = true, tableName = "table_1D_1", columns = 2:9, fileName = table_file);
equation
  Tab.u = h;
  T = Tab.y[1];
  rho = Tab.y[2];
  cp = Tab.y[3];
  mu = Tab.y[4];
  k = Tab.y[6];
end State;