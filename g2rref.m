% lifted from Github, Thanks to the author @benathon
% modified to find the binary inverse
% Irwan


% This is a modified version of matlab's building rref which calculates
% row-reduced echelon form in gf(2).  Useful for linear codes.
% Tolerance was removed because yolo, and because all values
% should only be 0 or 1.  @benathon

function C = g2rref(A)
%G2RREF   Reduced row echelon form in gf(2).
%   R = RREF(A) produces the reduced row echelon form of A in gf(2).
%
%   Class support for input A:
%      float: with values 0 or 1
%   Copyright 1984-2005 The MathWorks, Inc. 
%   $Revision: 5.9.4.3 $  $Date: 2006/01/18 21:58:54 $

[m,n] = size(A);
B = eye(m); ### added by Irwan

% Loop over the entire matrix.
i = 1;
j = 1;

while (i <= m) && (j <= n)
   % Find value and index of largest element in the remainder of column j.
   k = find(A(i:m,j),1) + i - 1;

   % Swap i-th and k-th rows.
   A([i k],j:n) = A([k i],j:n);
   B([i k],j:n) = B([k i],j:n); ### added by Irwan
   
   % Save the right hand side of the pivot row
   aijn = A(i,j:n);
   bijn = B(i,j:n); ### added by Irwan
   
   % Column we're looking at
   col = A(1:m,j);
   col2 = B(1:m,j); ### added by Irwan
      
   % Never Xor the pivot row against itself
   col(i) = 0;
   
   % This builds an matrix of bits to flip
   flip = col*aijn;
   flip2 = col*bijn;  ### added by Irwan
   
   % Xor the right hand side of the pivot row with all the other rows
   A(1:m,j:n) = xor( A(1:m,j:n), flip );
   B(1:m,j:n) = xor( B(1:m,j:n), flip ); ### added by Irwan
   i = i + 1;
   j = j + 1;
end

C = B; ### binary inverse of A (added by Irwan)
endfunction
