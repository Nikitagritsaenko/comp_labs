function [line] = find_line_eq(x, y)
  a = x(2) - y(2);
  b = y(1) - x(1);
  c = x(1) * y(2) - y(1) * x(2);
  mu = sqrt(a ^ 2 + b ^ 2);
  a = a / mu;
  b = b / mu;
  c = c / mu;
  
  if (c > 0)
    a = a * -1;
    b = b * -1;
  end
  c = abs(c);
  
  line = [a, b, c];
end
