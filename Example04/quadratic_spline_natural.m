function M = quadratic_spline_natural(x,y,n,h)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Input parametets:  x - data vector of length (n + 1) %
  %                    y - data vector of length (n + 1) %
  %                    n - see previous vectors          %
  %                    h - diff(x)                       %
  %                                                      %
  % Output parameters: M - computed spline coefficients  %
  %                                                      %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  sigma    = 2*diff(y)./h;
  d_tilde  = diff(sigma);
  A_tilde  = gallery("tridiag", ones(n-2,1), 2*ones(n-1,1), ones(n-2,1));
  M_tilde  = A_tilde\d_tilde;
  M        = [0; M_tilde; 0];

endfunction
