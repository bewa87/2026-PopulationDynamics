function M = cubic_spline_natural(x,y,n,h)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Input parametets:  x - data vector of length (n + 1) %
  %                    y - data vector of length (n + 1) %
  %                    n - see previous vectors          %
  %                    h - diff(x)                       %
  %                                                      %
  % Output parameters: M - computed spline coefficients  %
  %                                                      %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Step 1: Storing needed parameters

  lambda_tilde  = zeros(n-1,1);
  mu_tilde      = zeros(n-1,1);
  sigma         = zeros(n,1);
  d_tilde       = zeros(n-1,1);

  % Step 2: Computation of helping parameters

  for j = 1:1:(n-1)
    lambda_tilde(j)  = h(j+1)/(h(j) + h(j+1));
    mu_tilde(j)      = 1 - lambda_tilde(j);
  endfor

  sigma = diff(y)./h;

  for j = 1:1:(n-1)
    d_tilde(j) = 6*(sigma(j+1) - sigma(j))/(h(j) + h(j + 1));
  endfor

  % Step 3: Final computation of M - spline coefficients

  lambda = [0; lambda_tilde];
  mu     = [mu_tilde; 0];
  d      = [0; d_tilde; 0];

  A      = gallery("tridiag", mu, 2*ones(n+1,1), lambda);
  M      = A\d;

endfunction
