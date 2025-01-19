  function k_Tf "find thermal conductivity from temperature"
    input SI.Temperature T;
    input Real f;
    output SI.ThermalConductivity k;
  protected
    Real T_data[8] = {373.15, 473.15, 673.15, 873.15, 1073.15, 1273.15, 1473.15, 1673.15};
    Real k_data[8] = {4.85344, 4.35136, 3.68192, 3.3472, 3.17984, 3.01248, 3.01248, 3.01248};
  algorithm
    k := SolarTherm.Utilities.Interpolation.Interpolate1D(T_data,k_data,T);
  end k_Tf;