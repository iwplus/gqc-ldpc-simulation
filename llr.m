function b = llr(yi,t,n);
      if yi == 1,
        b = log(t/(n-t)); ### LLR value in Binary symmetric channel
      else,
        b = log((n-t)/t); ############################################
      endif
endfunction
