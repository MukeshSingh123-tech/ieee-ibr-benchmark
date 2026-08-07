function Spu = gb_inverter_rating(mpc, genRows)
%GB_INVERTER_RATING  MVA rating of the inverter replacing the given machines, pu.
%
%   MATLAB twin of python/gridbench/ibr.py:inverter_rating_pu.
%
%   Under the default 'match_machine_capability' sizing, S = sqrt(PMAX^2+QMAX^2).
%   Substituting that into Qmax(V,P) = sqrt((V*Ilim*S)^2 - P^2) at V = 1.0,
%   Ilim = 1.0, P = PMAX returns exactly QMAX -- so the replacement inverter
%   starts life with precisely the reactive capability of the machine it
%   displaces, and every later difference is attributable to the VOLTAGE AND
%   POWER DEPENDENCE of the converter limit rather than to a sizing assumption.
%
%   This matters: an independently inferred nameplate once handed the inverter
%   MORE capability than the machine it replaced (33 vs 24 MVAr on IEEE 14-bus),
%   which silently turned the study into a measurement of the sizing choice.

define_constants;
cfg = gb_config();

pmax = mpc.gen(genRows, PMAX);
qmax = abs(mpc.gen(genRows, QMAX));

% the slack's PMAX is a sentinel (1e9); fall back to its dispatch
bad = ~isfinite(pmax) | pmax >= 1e6;
pmax(bad) = abs(mpc.gen(genRows(bad), PG));
qmax(~isfinite(qmax) | qmax >= 1e6) = 0;

if strcmpi(cfg.ibr.sizing, 'match_machine_capability')
    rating = sqrt(sum(pmax.^2 + qmax.^2));
else
    rating = sum(pmax);
end

if ~isfinite(rating) || rating <= 0
    rating = max(sum(abs(mpc.gen(genRows, PG))), 1);
end
Spu = rating / mpc.baseMVA;
end
