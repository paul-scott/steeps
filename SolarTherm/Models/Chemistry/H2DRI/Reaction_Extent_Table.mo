within SolarTherm.Models.Chemistry.H2DRI;

model Reaction_Extent_Table
  parameter String file "File where extent matrix is stored.";
  
  //Inputs
  SI.Temperature T_rxn "Reaction temperature (K)"; //horizontal table axis u2
  Real r_input "Ratio of H2 to Fe2O3 (-)"; //vertical table axis u1
  
  //Raw Outputs
  Real xi_1_raw "Extent of reaction one";
  Real xi_2_raw "Extent of reaction two";
  Real xi_3_raw "Extent of reaction three";
  
  //Real Outputs
  Real xi_1 "Extent of reaction one after checking bounds";
  Real xi_2 "Extent of reaction two after checking bounds";
  Real xi_3 "Extent of reaction three after checking bounds";
  
  Modelica.Blocks.Tables.CombiTable2D xi_1_table(
    tableOnFile=true,
    tableName="xi_1",
    smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
    fileName=file);
  Modelica.Blocks.Tables.CombiTable2D xi_2_table(
    tableOnFile=true,
    tableName="xi_2",
    smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
    fileName=file);
  Modelica.Blocks.Tables.CombiTable2D xi_3_table(
    tableOnFile=true,
    tableName="xi_3",
    smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
    fileName=file);
algorithm    
    xi_1 := min(max(0.0,xi_1_raw),1.0);// 0.0 <= xi_1 <= 1.0
    xi_2 := min(max(0.0,xi_2_raw),xi_1);// 0.0 <= xi_2 <= xi_1
    xi_3 := min(max(0.0,xi_3_raw),xi_2);// 0.0 <= xi_3 <= xi_2
equation
    xi_1_table.u1 = r_input;
    xi_1_table.u2 = T_rxn;
    xi_1_raw = xi_1_table.y;
    
    xi_2_table.u1 = r_input;
    xi_2_table.u2 = T_rxn;
    xi_2_raw = xi_2_table.y;
    
    xi_3_table.u1 = r_input;
    xi_3_table.u2 = T_rxn;
    xi_3_raw = xi_3_table.y;
    
annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)));
end Reaction_Extent_Table;