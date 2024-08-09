within SolarTherm.Validation.Gan_HBS.Materials;

package Air_431kPa
  extends SolarTherm.Materials.PartialMaterial(MM = 56.077e-3, T_melt = 3200.0, cost = 0.2);//Composition {N2 = 78.0 mol%, O2 = 22.0 mol%} Pressure = 5884 Pa
  annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)));
end Air_431kPa;