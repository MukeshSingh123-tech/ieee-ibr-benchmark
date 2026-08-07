function gfmBuses = gb_split_gfm(mpc, ibrGens, gfmSharePct)
%GB_SPLIT_GFM  Which IBR buses are grid-forming, for a given fleet share.
%
%   gfmBuses = GB_SPLIT_GFM(mpc, ibrGens, 25)
%
%   MATLAB twin of python/gridbench/ibr.py:split_gfl_gfm. Share is measured in
%   MVA capacity, matching how an operator would write a grid-forming
%   requirement into a grid code.
%
%   Largest units first: if you are going to mandate grid-forming capability on
%   only part of a fleet, the large plants buy the most system strength per unit
%   of cost. That choice is what makes 25% GFM sufficient to restore a solvable
%   fault problem at 100% IBR penetration -- see wp3_faults_matlab.

define_constants;
gfmBuses = [];
if isempty(ibrGens) || gfmSharePct <= 0, return; end

w = mpc.gen(:, PMAX);
bad = ~isfinite(w) | w <= 0 | w >= 1e6;
w(bad) = abs(mpc.gen(bad, PG));
w(~isfinite(w) | w < 0) = 0;

if gfmSharePct >= 100
    gfmBuses = unique(mpc.gen(ibrGens, GEN_BUS));
    return
end

total = sum(w(ibrGens));
if total <= 0, return; end

[~, s] = sort(w(ibrGens), 'descend');
pool = ibrGens(s);

target = gfmSharePct / 100 * total;
acc = 0; chosen = [];
for g = pool(:)'
    if acc >= target - 1e-9, break; end
    chosen(end+1) = g; %#ok<AGROW>
    acc = acc + w(g);
end
gfmBuses = unique(mpc.gen(chosen, GEN_BUS));
end
