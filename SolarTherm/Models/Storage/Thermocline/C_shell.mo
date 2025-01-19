within SolarTherm.Models.Storage.Thermocline;

function C_shell "Estimating the cost of a cylindrical shell based on hydrostatic pressure."
import CN = Modelica.Constants;
  input Real rho_bulk "bulk density";
  input Real H "tank height";
  input Real D "tank diameter";
  input Real sigma_yield "yield strength of metal";
  input Real rho_shell "tank metal density";
  input Real pri_shell "tank USD/kg metal";
  input Real f_safety = 4 "safety factor of yield stress/hoop stress";
  output Real C "total cost of tank shell";
protected
  Real t_0 = 6e-3 "minimum shell thickness";
  Real t_top = 6e-3 "roof thickness";
  Real t_bot = 8e-3 "base thickness";
  Real t_1 "thickness of wall at the base";
  Real p "pressure at the base wall";
  Real V_shell "volume of shell material used";
algorithm
  p := rho_bulk*9.81*H;
  t_1 := max(t_0,p*D*f_safety/(2.0*sigma_yield)); //Cannot be smaller than t_0
  V_shell := CN.pi*(H/6)*(-2.0*t_0*t_0 +4.0*t_0*t_1+4.0*t_1*t_1 +3.0*D*(t_0+t_1)) + 0.25*CN.pi*D*D*(t_top+t_bot);
  C := rho_shell*V_shell*pri_shell;
end C_shell;