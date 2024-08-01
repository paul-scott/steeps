within SolarTherm.Calculators;

model Annual_Solar_Resource
//Calculates the annual solar resource (integral of DNI) based on a given TMY weather data file
//Alice Springs 1996 "example_TMY3.motab" lat=-23.795, lon=133.889, tzone=9.5
//Daggett 2008 "Daggett_Ca_TMY3.motab" lat=34.85, lon=-116.78, tzone=-8.0
//Pilbara "Pilbara_WA_2019.motab"
//Gladstone "Gladstone_QLD_2014.motab"
//Burnie "Burnie_TAS_2015.motab"

  parameter String wea_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Weather/example_TMY3.motab");
  //parameter String wea_file2 = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Weather/Pilbara_WA_2019.motab");
  parameter nSI.Angle_deg lon = -23.795 "Longitude (+ve East)";
  parameter nSI.Angle_deg lat = 133.889 "Latitude (+ve North)";
  parameter nSI.Time_hour t_zone = 9.5 "Local time zone (UCT=0)";
  parameter Integer year = 2008 "Meteorological year";
  Real DNI_Annual(start=0) "Solar resource, continuous smoothing (kWh/m2-year)";
  Real GHI_Annual(start=0) "Solar resource, linear smoothing (kWh/m2-year)";
//SolarTherm.Models.Sources.DataTable.DataTable data_1(lon = lon, lat = lat, t_zone = t_zone, year = year, file = wea_file) annotation(
//    Placement(visible = true, transformation(extent = {{-120, 82}, {-90, 110}}, rotation = 0)));
SolarTherm.Models.Sources.DataTable.DataTable data_2(lon = lon, lat = lat, t_zone = t_zone, year = year, file = wea_file) annotation(
    Placement(visible = true, transformation(extent = {{-120, 82}, {-90, 110}}, rotation = 0)));
SolarTherm.Models.Sources.DataTable.DataTable_Linear data_1(lon = lon, lat = lat, t_zone = t_zone, year = year, file = wea_file) annotation(
    Placement(visible = true, transformation(extent = {{-120, 82}, {-90, 110}}, rotation = 0)));
equation
  der(DNI_Annual) = max(0.0,data_1.DNI)/(3.6e6);
  der(GHI_Annual) = max(0.0,data_1.GHI)/(3.6e6);
annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)),experiment(StopTime = 3.1536e+07, StartTime = 0, Tolerance = 0.001, Interval = 60, maxStepSize = 60, initialStepSize = 60));
end Annual_Solar_Resource;