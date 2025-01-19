within SolarTherm.Utilities.Nusselt.Internal_Flow;
//Laminar Gas: Worsoe-Schmidt (1966), Int. J. Heat Mass Transfer, Vol. 9 pp. 1291-1295.
//
//Turbulent Gas: Humble, Lowdermilk and Desmon (1951), NACA Report 1020.
//30 < L/d < 120, 0.46 < T_s/T_f < 3.5, 7e3 < Re < 3e5
function Nusselt_Gas
  input Real Re "Reynolds";
  input Real Pr "Prandtl";
  input Real T_f "Fluid temperature";
  input Real T_s "Wall temperature";
  output Real Nu "Nusselt number";
algorithm
  if Re <= 2000 then //Laminar
    Nu := 3.66;
  elseif Re >= 4000 then //Turbulent
    if T_s > T_f then
      Nu := 0.023*(Re^0.8)*(Pr^0.4)*((max(T_s/T_f,1.0e-3))^(-0.55));
    else
      Nu := 0.023*(Re^0.8)*(Pr^0.4);
    end if;
  else //Transition
    if T_s > T_f then //Fluid is being heated
      Nu := 3.66 + ((0.023*(Re^0.8)*(Pr^0.4)*((max(T_s/T_f,1.0e-3))^(-0.55))) - 3.66)*SolarTherm.Utilities.Phis((Re - 2000.0) / 2000.0);
    else //Fluid is being cooled
      Nu := 3.66 + ((0.023*(Re^0.8)*(Pr^0.4)) - 3.66)*SolarTherm.Utilities.Phis((Re - 2000.0) / 2000.0);
    end if;
  end if;

annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)));
end Nusselt_Gas;