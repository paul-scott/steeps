within SolarTherm.Validation.Gan_HBS.Materials;

package Flue_Gas_107kPa
  //Composition {CO2 = 30.7 mol%, H2O = 1.8 mol%, N2 = 67.2 mol%, O2 = 0.3 mol%} Pressure = 107209 Pa
  extends SolarTherm.Materials.PartialMaterial(MM = 56.077e-3, T_melt = 3200.0, cost = 0.2);
  annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)));
end Flue_Gas_107kPa;