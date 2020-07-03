
function [theta_2] = snells(v_1, v_2, theta_1)
  theta_2 = asin( (v_2/v_1)*sin(theta_1) );
end