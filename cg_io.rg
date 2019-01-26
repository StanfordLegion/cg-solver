import "regent"

require "io"
require "bit"

local c = regentlib.c
local cstring = terralib.includec("string.h")
local cunistd = terralib.includec("unistd.h")

terra read_int(f : &c.FILE, nbytes : int32) : int32
  var vari : uint32[1];
  if c.fread(&vari[0], nbytes, 1, f) ~= 1 then
    c.printf("\n\n\nEMERGENCY\n\n\n");
  end
  return vari[0]
end

terra read_float(f : &c.FILE) : float
  var varf : float[1];
  if c.fread(&varf[0], 4, 1, f) ~= 1 then
    c.printf("\n\n\nEMERGENCY\n\n\n");
  end
  return varf[0]
end

terra read_double(f : &c.FILE) : double
	var varf : double[1];
	if c.fread(&varf[0], 8, 1, f) ~= 1 then
		c.printf("\n\n\nEMERGENCY\n\n");
	end
	--c.printf(" %f", varf[0])
	return varf[0]
end

