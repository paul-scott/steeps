within SolarTherm.Models.Storage.Thermocline;
//Linear interpolation of an array.
function Interpolate_Temperature "Input normalized tank height, output interpolated fluid temperature"
  input Real[:] ZDH_array "Array of tank vertical axis";
  input Real[:] T_array "Array of tank temperatures";
  input Integer N "Number of elements in the array";
  input Real Z "Normalised tank axis value where temperature is wanted";
  output Real T "Temperature output";
protected
  Integer i=1 "Iteration number";
algorithm
  if ZDH_array[1] > Z then //temperature is at the first element
    T := T_array[1];
  elseif ZDH_array[N] <= Z then //temperature is at the last element
    T := T_array[N];
  else //it is somewhere in the middle
    while i < N loop
      if ZDH_array[i] <= Z and ZDH_array[i+1] > Z then
        T := T_array[i] + ((Z-ZDH_array[i])/(ZDH_array[i+1]-ZDH_array[i]))*(T_array[i+1]-T_array[i]); //interpolation
        break;
      else
        i := i + 1;
      end if;
    end while;
  end if;
end Interpolate_Temperature;