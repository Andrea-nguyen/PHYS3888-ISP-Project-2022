function y = gamModel(t,a,b)
y = (b^a / gamma(a)) .* (t .^ (a-1)) .* exp(-b*t);
end