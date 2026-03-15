function [A, B] = linear_spline_natural(x,y,n,h)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Input parametets:  x - data vector of length (n + 1) %
  %                    y - data vector of length (n + 1) %
  %                    n - see previous vectors          %
  %                    h - diff(x)                       %
  %                                                      %
  % Output parameters: M - computed spline coefficients  %
  %                                                      %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  A = diff(y)./h;
  B = y(2:1:end) - A.*x(2:1:end);

endfunction
