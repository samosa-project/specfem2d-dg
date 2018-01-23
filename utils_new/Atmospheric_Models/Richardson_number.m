function [R] = Richardson_number(WIND, D, N)
  R = (N ./ (D * WIND)) .^ 2.0;
end

