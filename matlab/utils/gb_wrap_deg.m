function out = gb_wrap_deg(angle_deg)
%GB_WRAP_DEG  Wrap an angle (or array) into (-180, +180] degrees.
%
%   Twin of python/gridbench/metrics.py:wrap_deg.
%
%   Required for every angle comparison in this project. A relay that sees
%   -161.4 deg and 177.0 deg differs by -21.5 deg, not +338.5 deg; reporting the
%   unwrapped value turns a modest, plausible shift into a spectacular but
%   meaningless one.

out = mod(angle_deg + 180, 360) - 180;
out(abs(out + 180) < eps(180)) = 180;    % map the -180 endpoint to +180
end
