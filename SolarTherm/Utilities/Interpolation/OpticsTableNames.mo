within SolarTherm.Utilities.Interpolation;

function OpticsTableNames "input the prefix to a optical lookup table, outputs an array of filenames for two optical tables"
    input String prefix;
    input Real f_recv;
    output String[2] file_optics_list;
protected
    Real pct;
    Integer count;
    Integer[13] List = {70,75,80,85,90,95,100,105,110,115,120,125,130};
    Boolean success;
    Real weight;
    String file_opticsA;
    String file_opticsB;
    
  algorithm
    success := false;
    count := 1;
    pct := 100.0 * f_recv;
    if pct > 130.0 then
      file_opticsA := prefix + "_125%Arecv_optics.motab";
      file_opticsB := prefix + "_130%Arecv_optics.motab";
      weight := (pct - 125) / (130 - 125);
    elseif pct < 70.0 then
      file_opticsA := prefix + "_70%Arecv_optics.motab";
      file_opticsB := prefix + "_75%Arecv_optics.motab";
      weight := (pct -70) / (75 - 70);
    else
      while count < 13 and success == false loop
        if pct >= List[count] and pct <= List[count + 1] then
          file_opticsA := prefix + "_" + String(List[count]) + "%Arecv_optics.motab";
          file_opticsB := prefix + "_" + String(List[count + 1]) + "%Arecv_optics.motab";
          weight := (pct - List[count]) / (List[count + 1] - List[count]);
          success := true;
        else
          count := count + 1;
        end if;
      end while;
    end if;
    file_optics_list := {file_opticsA, file_opticsB};
  end OpticsTableNames;