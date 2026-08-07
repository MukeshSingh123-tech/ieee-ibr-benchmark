function p = gb_root()
%GB_ROOT  Absolute path to the IEEE/ project root.
%
%   Used so scripts can write to results/ and read data/ regardless of the
%   current working directory when they are invoked.

p = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end
