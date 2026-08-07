function gb_writetable(rows, filename)
%GB_WRITETABLE  Write a struct array to CSV with the same schema Python uses.
%
%   GB_WRITETABLE(rows, filename)
%
%   Column names must match python/gridbench/metrics.py output so that
%   gb_compare_python and python/studies/compare_matlab.py can join the two
%   sides on identical keys. Logical values are written as the strings
%   "True"/"False" to match Python's csv output rather than MATLAB's 1/0.

if isempty(rows)
    fid = fopen(filename, 'w'); fclose(fid);
    return
end

fields = fieldnames(rows);
T = table();
for k = 1:numel(fields)
    f = fields{k};
    vals = {rows.(f)};

    if all(cellfun(@(v) islogical(v) && isscalar(v), vals))
        col = strings(numel(vals), 1);
        for i = 1:numel(vals)
            if vals{i}, col(i) = "True"; else, col(i) = "False"; end
        end
    elseif all(cellfun(@(v) isnumeric(v) && isscalar(v), vals))
        col = cell2mat(vals(:));
    else
        col = string(cellfun(@(v) toStr(v), vals(:), 'UniformOutput', false));
    end
    T.(f) = col;
end

writetable(T, filename);
end


function s = toStr(v)
if ischar(v) || isstring(v)
    s = char(v);
elseif isnumeric(v) && isscalar(v)
    s = num2str(v);
elseif isnumeric(v)
    s = strjoin(arrayfun(@num2str, v(:)', 'UniformOutput', false), ' ');
elseif islogical(v)
    if v, s = 'True'; else, s = 'False'; end
else
    s = '';
end
end
