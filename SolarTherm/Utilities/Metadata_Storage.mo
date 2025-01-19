within SolarTherm.Utilities;

function Metadata_Storage
    import Modelica.Utilities.Streams;
    import Modelica.Utilities.Strings;
    input String file;
    output Real[21] meta;
  protected
    String line;
    Integer i;
    Integer j;
  algorithm
    line := Streams.readLine(file, 5);
    Streams.close(file);

    i := 10;
    j := 10;
    (meta[1], j) := Strings.scanReal(line, j+1, message="Metadata cannot be read " + file);
    (meta[2], j) := Strings.scanReal(line, j+1, message="Metadata cannot be read " + file);
    (meta[3], j) := Strings.scanReal(line, j+1, message="Metadata cannot be read " + file);
    (meta[4], j) := Strings.scanReal(line, j+1, message="Metadata cannot be read " + file);
    (meta[5], j) := Strings.scanReal(line, j+1, message="Metadata cannot be read " + file);
    (meta[6], j) := Strings.scanReal(line, j+1, message="Metadata cannot be read " + file);
    (meta[7], j) := Strings.scanReal(line, j+1, message="Metadata cannot be read " + file);
    (meta[8], j) := Strings.scanReal(line, j+1, message="Metadata cannot be read " + file);
    (meta[9], j) := Strings.scanReal(line, j+1, message="Metadata cannot be read " + file);
    (meta[10], j) := Strings.scanReal(line, j+1, message="Metadata cannot be read " + file);
    (meta[11], j) := Strings.scanReal(line, j+1, message="Metadata cannot be read " + file);
    (meta[12], j) := Strings.scanReal(line, j+1, message="Metadata cannot be read " + file);
    (meta[13], j) := Strings.scanReal(line, j+1, message="Metadata cannot be read " + file);
    (meta[14], j) := Strings.scanReal(line, j+1, message="Metadata cannot be read " + file);
    (meta[15], j) := Strings.scanReal(line, j+1, message="Metadata cannot be read " + file);
    (meta[16], j) := Strings.scanReal(line, j+1, message="Metadata cannot be read " + file);
    (meta[17], j) := Strings.scanReal(line, j+1, message="Metadata cannot be read " + file);
    (meta[18], j) := Strings.scanReal(line, j+1, message="Metadata cannot be read " + file);
    (meta[19], j) := Strings.scanReal(line, j+1, message="Metadata cannot be read " + file);
    (meta[20], j) := Strings.scanReal(line, j+1, message="Metadata cannot be read " + file);
    (meta[21], j) := Strings.scanReal(line, j+1, message="Metadata cannot be read " + file);
end Metadata_Storage;