within SolarTherm.Utilities.Nusselt.Internal_Flow;
//Laminar Liquid: Yang (1962), J. Heat Transfer, Vol 84, pp. 353-362.
//
//Turbulent Liquid: Petukhov & Kirillov (1970), Advances in Heat Transfer, Vol 6, 504-564.
// L/d > 60, 0.08 < mu_s/mu_f < 40, 1e4 < Re < 5e6, 2 < Pr < 140
function Nusselt_Liquid
  input Real Re "Reynolds";
  input Real Pr "Prandtl";
  input Real T_f "Fluid temperature";
  input Real T_s "Wall temperature";
  input Real mu_f "Fluid viscosity";
  input Real mu_s "Fluid viscosity at the wall temperature";
  output Real Nu "Nusselt number";
protected
  Real f "Friction Factor";
algorithm
  if Re <= 2000 then //Laminar
    Nu := 3.66*((max(mu_f/mu_s,1.0e-3))^(0.11));
    f := 64.0/Re;
  elseif Re >= 4000 then //Turbulent
    f := (1.82*log10(Re)-1.64)^(-2.0);
    if T_s > T_f then //Turbulent
      Nu := (((f/8.0)*Re*Pr) / ( 1.07 + 12.7*((f/8.0)^0.5)*((Pr^0.667)-1.0)))*((max(mu_f/mu_s,1.0e-3))^(0.11));
    else
      Nu := (((f/8.0)*Re*Pr) / ( 1.07 + 12.7*((f/8.0)^0.5)*((Pr^0.667)-1.0)))*((max(mu_f/mu_s,1.0e-3))^(0.25));
    end if;
  else //Transition
    f := (1.82*log10(Re)-1.64)^(-2.0);
    if T_s > T_f then //Fluid is being heated
      Nu := 3.66*((max(mu_f/mu_s,1.0e-3))^(0.11)) + ( (((f/8.0)*Re*Pr) / ( (1.07 + 12.7*((f/8.0)^0.5)*((Pr^0.667)-1.0)))*((max(mu_f/mu_s,1.0e-3))^(0.11))) - 3.66*((max(mu_f/mu_s,1.0e-3))^(0.11)) )*SolarTherm.Utilities.Phis((Re - 2000.0) / 2000.0);
    else //Fluid is being cooled
      Nu := 3.66*((max(mu_f/mu_s,1.0e-3))^(0.11)) + ( (((f/8.0)*Re*Pr) / ( (1.07 + 12.7*((f/8.0)^0.5)*((Pr^0.667)-1.0)))*((max(mu_f/mu_s,1.0e-3))^(0.25))) - 3.66*((max(mu_f/mu_s,1.0e-3))^(0.11)) )*SolarTherm.Utilities.Phis((Re - 2000.0) / 2000.0);
    end if;
  end if;

annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)));
end Nusselt_Liquid;