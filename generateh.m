function H = generateh(partisi,x,r);
  p = length(partisi);
  H = x;
  for j = 1:(r-1),
    y = [];
    for i = 1:p,
      if i == 1,
        x1 = shiftvektor(x(1:partisi(i)));
      else,
        ind1 = sum(partisi(1:(i-1)));
        ind2 = ind1+partisi(i);
        x1 = shiftvektor(x((ind1+1):ind2));
      endif
      y = [y x1];
    endfor
    x = y;
    H=[H;y];
  endfor
endfunction
