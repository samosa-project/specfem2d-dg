function [theta_out] = snells(v_inc, v_out, theta_inc)
  theta_out = asin( (v_out/v_inc)*sin(theta_inc) );
end