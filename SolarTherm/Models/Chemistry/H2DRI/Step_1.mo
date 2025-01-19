within SolarTherm.Models.Chemistry.H2DRI;

model Step_1
  extends Icons.ChemEqn;
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  import CV = Modelica.SIunits.Conversions;
  final SI.MolarInternalEnergy R = CN.R "Molar gas constant (J/molK)";
  
  //    A              B                                           C             D
  //0.5Fe2O3(s) + (1/6)H2(g) ---DeltaH_rxn (J/molrxn) ---> (1/3)Fe3O4(s) + (1/6)H2O at different temperatures.
  //H2O is at the most stable state at whatever temperature is defined.
  
  //Chemistry Tables
  Chemical Chem_A=SolarTherm.Models.Chemistry.ChemTable.Fe2O3;
  Chemical Chem_B=SolarTherm.Models.Chemistry.ChemTable.H2;//"Standard thermodynamic constants for B";
  Chemical Chem_C=SolarTherm.Models.Chemistry.ChemTable.Fe3O4;//"Standard thermodynamic constants for C";
  Chemical Chem_D=SolarTherm.Models.Chemistry.ChemTable.H2O;// "Standard thermodynamic constants for D";
  
  //Stochiometric Ratio
  parameter Real n_A = 0.5 "Moles of A per moles of reaction (mol/molrxn)";
  parameter Real n_B = 1.0/6.0 "Moles of B per moles of reaction (mol/molrxn)";
  parameter Real n_C = 1.0/3.0 "Moles of C per moles of reaction (mol/molrxn)";
  parameter Real n_D = 1.0/6.0 "Moles of D per moles of reaction (mol/molrxn)";
  
  //Reference conditions at which molar enthalpy of formation is defined at
  parameter SI.Temperature T_ref = 298.15 "Reference temperature 298.15 K";
  parameter SI.Pressure p_ref = 1.0e5 "Reference pressure 100000 Pa";
  
  //Imposed conditions
  SI.Temperature T "Reaction temperature (K)";
  SI.Pressure p "Reaction pressure (Pa)";

  //Outputs
  SI.MolarInternalEnergy H_rxn "Enthalpy change of reaction at temperature T and pressure p per mole of reaction (J/molrxn)";
  SI.MolarInternalEnergy G_rxn "Gibbs free energy of reaction at temperature T and pressure p per mole of reaction (J/molrxn)";
  SI.MolarEntropy S_rxn "Delta Entropy of reaction at temperature T and pressure p per mole of reaction (J/molrxnK)";
  Real K_rxn "Equilibrium constant (unspecified units)";

protected
  //Chemical Aspect (Thermodynamic properties at standard conditions)
  //Chemicals are declared in lines 15-18. Properties are passed on automatically.
  //A = Fe2O3(s)
  parameter SI.MolarMass M_A = Chem_A.M "Molar mass (kg/mol)"; 
  parameter SI.MolarEnthalpy Hf_A = Chem_A.Hf0 "Molar standard formation enthalpy (J/mol)";
  parameter SI.MolarEntropy Sf_A = Chem_A.S0 "Molar absolute entropy (J/molK)";
  
  //B = H2(g)
  parameter SI.MolarMass M_B = Chem_B.M "Molar mass (kg/mol)"; 
  parameter SI.MolarEnthalpy Hf_B = Chem_B.Hf0 "Molar standard formation enthalpy (J/mol)";
  parameter SI.MolarEntropy Sf_B = Chem_B.S0 "Molar absolute entropy (J/molK)";
  
  //C = Fe3O4(s)
  parameter SI.MolarMass M_C = Chem_C.M "Molar mass (kg/mol)"; 
  parameter SI.MolarEnthalpy Hf_C = Chem_C.Hf0 "Molar standard formation enthalpy (J/mol)";
  parameter SI.MolarEntropy Sf_C = Chem_C.S0 "Molar absolute entropy (J/molK)";
  
  //D = H20(l)
  parameter SI.MolarMass M_D = Chem_D.M "Molar mass (kg/mol)"; 
  parameter SI.MolarEnthalpy Hf_D = Chem_D.Hf0 "Molar standard formation enthalpy (J/mol)";
  parameter SI.MolarEntropy Sf_D = Chem_D.S0 "Molar absolute entropy (J/molK)";
  
  //Thermal Aspect (Deviation from standard temperature and pressure)
  //A = Fe2O3(s)
  SI.SpecificEnthalpy h_refA "Specific enthalpy of substance A at 298.15K (J/kg)";
  SI.SpecificEnthalpy h_A "Specific enthalpy of substance A at temperature T (J/kg)";
  SI.SpecificEntropy s_refA "Reference specific entropy of substance A at 298.15K (J/kgK)";
  SI.SpecificEntropy s_A "Specific entropy of substance A at temperature T (J/kgK)"; 
  
  //B = H2(g)
  SI.SpecificEnthalpy h_refB "Reference specific enthalpy of substance B at 298.15K (J/kg)";
  SI.SpecificEnthalpy h_B "Specific enthalpy of substance B at temperature T (J/kg)";
  SI.SpecificEntropy s_refB "Reference specific entropy of substance B at 298.15K (J/kgK)";
  SI.SpecificEntropy s_B "Specific entropy of substance B at temperature T (J/kgK)";
  
  //C = Fe3O4(s)
  SI.SpecificEnthalpy h_refC "Reference specific enthalpy of substance C at 298.15K (J/kg)";
  SI.SpecificEnthalpy h_C "Specific enthalpy of substance C at temperature T (J/kg)";
  SI.SpecificEntropy s_refC "Reference specific entropy of substance C at 298.15K (J/kgK)";
  SI.SpecificEntropy s_C "Specific entropy of substance C at temperature T (J/kgK)";
  
  //D = H20(l)
  SI.SpecificEnthalpy h_refD "Reference specific enthalpy of substance D at 298.15K (J/kg)";
  SI.SpecificEnthalpy h_D "Specific enthalpy of substance D at temperature T (J/kg)";
  SI.SpecificEntropy s_refD "Reference specific entropy of substance D at 298.15K (J/kgK)";
  SI.SpecificEntropy s_D "Specific entropy of substance D at temperature T (J/kgK)";

  //Molar Enthalpies
  SI.MolarEnthalpy H_chem_reactants "Chemical enthalpy of reactants per reaction (J/molrxn)";
  SI.MolarEnthalpy H_chem_products "Chemical enthalpy of products per reaction (J/molrxn)";  
  SI.MolarEnthalpy H_therm_reactants "Thermal enthalpy of reactants per reaction (J/molrxn)";
  SI.MolarEnthalpy H_therm_products "Thermal enthalpy of products per reaction (J/molrxn)";
  
  //Molar Entropies
  SI.MolarEntropy S_chem_reactants "Standard entropy of reactants per reaction (J/Kmolrxn)";
  SI.MolarEntropy S_chem_products "Standard entropy of products per reaction (J/Kmolrxn)";
  SI.MolarEntropy S_therm_reactants "Delta Entropy of reactants per reaction (J/Kmolrxn)";
  SI.MolarEntropy S_therm_products "Delta Entropy of products per reaction (J/Kmolrxn)";
  
equation
  //Specific enthalpy and entropy differences from (T_ref,p_ref) to (T,p) (J/kg)
  //Unfortunately, due to lack of standardisation between Modelica media models, these functions have to be included manually.
  //A = Fe2O3(s)
  h_refA = SolarTherm.Media.SolidParticles.Fe2O3_utilities.h_T(T_ref);
  h_A = SolarTherm.Media.SolidParticles.Fe2O3_utilities.h_T(T);
  s_refA = SolarTherm.Media.SolidParticles.Fe2O3_utilities.s_T(T_ref);
  s_A = SolarTherm.Media.SolidParticles.Fe2O3_utilities.s_T(T);
  
  //B = H2(g)
  h_refB = Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_ref,T_ref);  
  h_B = Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p,T);
  s_refB = Modelica.Media.IdealGases.SingleGases.H2.specificEntropy_pTX(p_ref,T_ref,{1.0}); //Mass fraction X is 1.0 for a single gas
  s_B = Modelica.Media.IdealGases.SingleGases.H2.specificEntropy_pTX(p,T,{1.0});
  
  //C = Fe3O4(s)
  h_refC = SolarTherm.Media.SolidParticles.Fe3O4_utilities.h_T(T_ref);
  h_C = SolarTherm.Media.SolidParticles.Fe3O4_utilities.h_T(T);
  s_refC = SolarTherm.Media.SolidParticles.Fe3O4_utilities.s_T(T_ref);
  s_C = SolarTherm.Media.SolidParticles.Fe3O4_utilities.s_T(T);
  
  //D = H20(l/g)
  h_refD = Modelica.Media.Water.IF97_Utilities.h_pT(p_ref,T_ref);
  h_D = Modelica.Media.Water.IF97_Utilities.h_pT(p,T);
  s_refD = Modelica.Media.Water.IF97_Utilities.s_pT(p_ref,T_ref);
  s_D = Modelica.Media.Water.IF97_Utilities.s_pT(p,T);
  
  //Sum up the molar enthalpies (J/molrxn)
  H_chem_reactants = n_A*Hf_A + n_B*Hf_B;
  H_chem_products = n_C*Hf_C + n_D*Hf_D;
  H_therm_reactants = n_A*(h_A-h_refA)*M_A + n_B*(h_B-h_refB)*M_B;
  H_therm_products = n_C*(h_C-h_refC)*M_C + n_D*(h_D-h_refD)*M_D;
  
  //Overall Delta_H, molar enthalpy change of the reaction (J/molrxn)
  H_rxn = H_chem_products + H_therm_products - H_chem_reactants - H_therm_reactants;
  
  //Sum up the molar entropies (J/molrxn K)
  S_chem_reactants = n_A*Sf_A + n_B*Sf_B;
  S_chem_products = n_C*Sf_C + n_D*Sf_D;
  
  S_therm_reactants = n_A*(s_A-s_refA)*M_A + n_B*(s_B-s_refB)*M_B;
  S_therm_products = n_C*(s_C-s_refC)*M_C + n_D*(s_D-s_refD)*M_D;
  
  //Overall Delta_S, molar entropy change of the reaction (J/molrxn K)
  S_rxn = S_chem_products - S_chem_reactants + S_therm_products - S_therm_reactants;
  
  //Overall Delta_G, molar Gibbs free energy change of the reaction (J/molrxn)
  G_rxn = H_rxn - T*(S_rxn);

  //Overall K, equilibrium constant of the reaction.
  K_rxn = exp(-1.0*G_rxn/(R*T));
annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)));
end Step_1;