within SolarTherm.Models.Storage;

function CpA_external_insulation "Input T_max (K) and U_loss (W/m2K), outputs cost per unit area (USD/m2)"
  input SI.Temperature T_max "Design maximum temperature (K)";
  input SI.CoefficientOfHeatTransfer U_loss_tank "Design heat loss coefficient (W/m2K)";
  output Real CpA "Cost per unit area of insulation (USD/m2)";
protected
  Real t "Theta, non-dimensional temperature as a ratio of 298.15K";
  Real R "1/U (m2K/W)";
  //LT: Rockwool Only
  Real a_0 = -0.146934;
  Real a_1 = 1.51788;
  Real a_2 = -0.353781;
  Real a_3 = -0.0109734;
  Real a_4 = 0.252587;
  //Regression via StatGraphics 4179 Observations, 240C - 640C at 20C intervals, 0.2 <= R <= 20 at 0.1 intervals
  //R2 = 100.0
  //AdjR2 = 100.0
  //Standard Error = 0.00758971
  //Mean Abs Error = 0.00582226
  
  //HT: Microporous + Rockwool
  Real b_0 = 223.307;
  Real b_1 = 259.08;
  Real b_2 = -790.825;
  Real b_3 = -61.2552;
  Real b_4 = 12.3795;
  //Regression via StatGraphics 1194 Observations, 660C - 760C at 20C intervals, 0.2 <= R <= 20 at 0.1 intervals
  //R2 = 99.9928
  //AdjR2 = 99.9928
  //Standard Error = 1.24006
  //Mean Abs Error = 1.0542
algorithm
  t := T_max/298.15;
  R := 1.0/U_loss_tank;
  
  if T_max > (650.0+273.15) then
    //LT Equation 1, R, R/t, 1/R, t^2*R
    CpA := a_0 + a_1*R + a_2*(R/t) + a_3*(1.0/R) + a_4*((t^2)*R);
  else
    //HT Equation 1, R, R/t, t^2, t^3
    CpA := b_0 + b_1*R + b_2*(R/t) + b_3*(t^2) + b_4*(t^3);
  end if;
annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)));
end CpA_external_insulation;