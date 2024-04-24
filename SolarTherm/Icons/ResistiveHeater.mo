within SolarTherm.Icons;
partial class ResistiveHeater
	extends Modelica.Icons.Package;
	annotation (
		Icon(coordinateSystem(preserveAspectRatio=false), 
			graphics={
			Ellipse(
				origin = {1, 0}, 
				fillColor = {255, 255, 255}, 
				fillPattern = FillPattern.Solid, 
				lineThickness = 0.25, 
				extent = {{-50, 50}, {50, -50}}), 
			Polygon(origin = {1, -3}, 
				fillColor = {252, 233, 79}, 
				fillPattern = FillPattern.Solid, 
				lineThickness = 0.25, 
				points = {{1, 49}, {-31, -5}, {5, -1}, {-1, -43}, {27, 17}, {-7, 11}, {1, 49}}),
			Text(
				origin = {-10, 254}, 
				lineColor = {0, 0, 255}, 
				extent = {{-149, -114}, {151, -154}}, 
				textString = "%name")
				}),
		Diagram(
			coordinateSystem(preserveAspectRatio=false)));
end ResistiveHeater;