within SolarTherm.Validation.Models.PCM;

function K
  input Real f;
  output Real k;
algorithm
  if f <= 0.0 then
    k := 4.0;
  elseif f >= 1.0 then
    k := 1.25;
  else
    k := (1.25*4.0*(((f-1.0)*Rho(1.0))-(f*Rho(0.0))))/(((f-1.0)*1.25*Rho(1.0))-(f*4.0*Rho(0.0)));
    //k := 4.0 + f*(1.25-4.00); //Do not use simple mean
  end if;
end K;