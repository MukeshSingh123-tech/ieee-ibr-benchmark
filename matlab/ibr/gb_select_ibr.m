function [gens, penActual] = gb_select_ibr(mpc, penetrationPct, order, preserveSlack)
%GB_SELECT_IBR  Choose which generator ROWS are displaced by IBRs.
%
%   [gens, penActual] = GB_SELECT_IBR(mpc, 60)
%
%   MATLAB twin of python/gridbench/ibr.py:select_ibr_gens.
%
%   Penetration is measured in installed MVA capacity -- how interconnection
%   studies state it -- not dispatched MW. That also handles synchronous
%   condensers correctly: a machine dispatching 0 MW still supplies fault
%   current and voltage support, so replacing it is a real change that a
%   dispatched-MW basis would score as zero.
%
%   Machines with zero weight are never selected. Greedily absorbing them was a
%   real bug on the Python side: it converted EVERY machine in the system at a
%   requested 20% penetration, because a zero-MW unit never advances the target.
%
%   The slack is preserved below 100%: a phasor load flow needs an angle
%   reference. At 100% it is converted too but REMAINS the reference bus, which
%   is the correct representation of a grid-forming inverter.

define_constants;
if nargin < 3 || isempty(order), order = 'largest_first'; end
if nargin < 4, preserveSlack = true; end

on   = mpc.gen(:, GEN_STATUS) > 0;
live = find(on);

refBuses  = mpc.bus(mpc.bus(:, BUS_TYPE) == REF, BUS_I);
isSlack   = ismember(mpc.gen(:, GEN_BUS), refBuses);

% capacity weight: PMAX, with the slack's 1e9 sentinel replaced by its dispatch
w = mpc.gen(:, PMAX);
bad = ~isfinite(w) | w <= 0 | w >= 1e6;
w(bad) = abs(mpc.gen(bad, PG));
w(~isfinite(w) | w < 0) = 0;

total = sum(w(live));
gens = []; penActual = 0;
if total <= 0 || penetrationPct <= 0, return; end

keepSlack = preserveSlack && penetrationPct < 100;
pool = live(w(live) > 0 & ~(keepSlack & isSlack(live)));
if isempty(pool), return; end

switch lower(order)
    case 'largest_first',  [~, s] = sort(w(pool), 'descend');
    case 'smallest_first', [~, s] = sort(w(pool), 'ascend');
    otherwise,             s = randperm(numel(pool));
end
pool = pool(s);

target = penetrationPct / 100 * total;
acc = 0; chosen = [];
for g = pool(:)'
    if acc >= target - 1e-9, break; end
    chosen(end+1) = g; %#ok<AGROW>
    acc = acc + w(g);
end

gens = sort(chosen(:));
penActual = sum(w(gens)) / total * 100;
end
