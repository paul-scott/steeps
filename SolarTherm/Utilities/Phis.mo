within SolarTherm.Utilities;
//Based on Armando Fontalvo's work
function Phis "Smooth transition function, x capped to a minimum of 1e-3"
  input Real x;
  output Real y;
protected
  Real x1;
algorithm
  x1 := max(1.0e-3, x);
  y := exp(-1.0/x1)/(exp(-1.0/x1) + exp(-1.0/(1.0-x1)));
annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)));
end Phis;