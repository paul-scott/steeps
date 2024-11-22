within SolarTherm.Utilities.Nusselt.Internal_Flow;

function Nusselt_HBS_Rough "Nusselt number based on Haaland approximation"
  input Real Re "Reynolds";
  input Real Pr "Prandtl";
  input Real E_div_d "Ratio of surface roughness to hydraulic diameter of channel (-)";
  //input Real T_f "Fluid temperature";
  //input Real T_s "Wall temperature";
  output Real Nu "Nusselt number";
protected
  Real f "Friction Factor";
algorithm
  if Re <= 2300.0 then //Laminar
    Nu := 3.66;
    f := 64.0/max(Re, 1.0e-6);
    
  elseif Re >= 4000.0 then //Turbulent
    //f := (0.790*log(Re) - 1.64)^(-2.0);
    f := 1.0/((-1.8*log10((6.9/max(Re, 1.0e-6) + (E_div_d/3.7)^1.11)))^2.0);
    Nu := ((f/8.0)*(Re-1000.0)*Pr)/(1.0+1.27*((f/8.0)^0.5)*((Pr^(2.0/3.0))-1.0));
    
  else //Transition
    f := 1.0/((-1.8*log10((6.9/4000.0 + (E_div_d/3.7)^1.11)))^2.0); //f_max
    
    Nu := 3.66 + SolarTherm.Utilities.Phis((Re - 2300.0) / 1700.0)*(((f/8.0)*(3000.0)*Pr)/(1.0+1.27*((f/8.0)^0.5)*((Pr^(2.0/3.0))-1.0)) - 3.66); //Nu = Nu_low + Phis*(Nu_high-Nu_low)
  end if;

annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)));
end Nusselt_HBS_Rough;