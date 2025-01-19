within SolarTherm.Utilities.Interpolation;

function Interpolate1D "Insert x_data[:], y_data[:],x, obtain a y value as an output based on linear interpolation or extrapolation at the end points"
    input Real x_data[:];
    input Real y_data[:];
    input Real x;
    output Real y;
  protected
    Integer N = size(x_data, 1);
    Integer i = 1;
  algorithm
    while i <= N loop
      if x < x_data[1] then
        y := y_data[1] + (x - x_data[1]) / (x_data[2] - x_data[1]) * (y_data[2] - y_data[1]);
        break;
      elseif x > x_data[N] then
        y := y_data[N - 1] + (x - x_data[N - 1]) / (x_data[N] - x_data[N - 1]) * (y_data[N] - y_data[N - 1]);
        break;
      elseif x >= x_data[i] and x <= x_data[i + 1] then
        y := y_data[i] + (x - x_data[i]) / (x_data[i + 1] - x_data[i]) * (y_data[i + 1] - y_data[i]);
        break;
      else
        i := i + 1;
      end if;
    end while;
//extrapolate
  end Interpolate1D;