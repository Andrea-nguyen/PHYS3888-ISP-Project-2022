function s = getScale(y,h)
if nargin < 2
    h = 1;
end
m = max(y);
s = h / m;
end