function y = shiftvektor(x);
  n = length(x);
for i= 1:n,
  if i == 1,
    y(i) = x(n);
  else,
  y(i) = x(i-1);
  endif 
endfor
endfunction
