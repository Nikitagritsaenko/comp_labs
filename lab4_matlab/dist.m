function [result] = dist(A, B)

if(nargin == 1)
    result = sqrt((A.x - 0)^2 + (A.y - 0)^2);
else
    result = sqrt((A.x - B.x)^2 + (A.y - B.y)^2);
end

end