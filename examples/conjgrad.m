function x = conjgrad(A,b,tol)
% CONJGRAD  Conjugate Gradient Method.
%   X = CONJGRAD(A,B) attemps to solve the system of linear equations A*X=B
%   for X. The N-by-N coefficient matrix A must be symmetric and the right
%   hand side column vector B must have length N.
%
%   X = CONJGRAD(A,B,TOL) specifies the tolerance of the method. The
%   default is 1e-10.
%
% Example (highlight lines between %{ and %}, then press F9):
%{
  n = 6000;
  m = 8000;
  A = randn(n,m);
  A = A * A';
  b = randn(n,1);
  tic, x = conjgrad(A,b); toc
  norm(A*x-b)
%}
% By Yi Cao at Cranfield University, 18 December 2008
% Updated on 6 Feb 2014.
%
% Obtained on 2/1/2018 from:
% https://www.mathworks.com/matlabcentral/fileexchange/22494-conjugate-gradient-method?focused=3809602&tab=function
% Some changes added by Richard Payne (denoted with RDP)
    if nargin<3
        tol=1e-10;
    end
    % Error checking from RDP
    if (size(b,1) ~= size(A,1)) || (size(b,2) > 1)
        error('b must be a column vector with length size(A,1)');
    end
        
    
    x = b;
    r = b - A*x;
    if norm(r) < tol
        return
    end
    y = -r;
    z = A*y;
    s = y'*z;
    t = (r'*y)/s;
    x = x + t*y;
 
    %oldnorm = 1; % RDP
    for k = 1:numel(b)
       r = r - t*z;
       thenorm = norm(r); % RDP addition
       if (thenorm < tol) 
            return;
       end
       %oldnorm = thenorm; % RDP
       B = (r'*z)/s;
       y = -r + B*y;
       z = A*y;
       s = y'*z;
       t = (r'*y)/s;
       x = x + t*y;
    end
 end