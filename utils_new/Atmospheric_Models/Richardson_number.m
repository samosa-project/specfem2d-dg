function [R] = Richardson_number(WIND, D, N)
  % W (m/s)
  % N (rad/s)
  R = (N ./ (D * WIND)) .^ 2.0;
end

