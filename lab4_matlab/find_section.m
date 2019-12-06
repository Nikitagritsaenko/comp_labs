% -----------------------------
%   DATE 2018.08.02
%   getYZSection.m
%
%   Функция находит проекцию фигуры заданной в XZ на плоскость x = h
%
%   [in] inSet - фигура, заданная в XZ
%   [in] h - задает плоскость YZ
%   [out] outSet - проекция фигуры на плоскость x = h
function outSet = find_section(inSet, h)
outSet = [];
inSet = squeeze(inSet);
k = 1;
for i = 1:length(inSet)
    r = inSet(i,1);
    z = inSet(i,2);
    a = 1;
    b = 0;
    c = h;
    s_len = sqrt(a * a + b * b);
    x0 = a * c / s_len;
    y0 = b * c / s_len;
    d = r * r - c * c / s_len + 0.01;
    if d < -1e-9
      %  outSet = [];
    elseif abs(d) < 1e-9
        outSet(k, :) = [y0 z];
        k = k + 1;
    else
        delta = sqrt(d / s_len);
        x1 = x0 + b * delta;
        x2 = x0 - b * delta;
        y1 = y0 - a * delta;
        y2 = y0 + a * delta;
        outSet(k, :) = [y2 z];
      %  outSet(4+k, :) = [y2 z];
        k = k + 1;
    end 
    if length(outSet) == 3
       outSet(4, :) = outSet(1, :);
    end
end

