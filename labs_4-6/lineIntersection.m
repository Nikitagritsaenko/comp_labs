function [intersection] = lineIntersection(points, line, closure)
a = line(1);
b = line(2);
c = line(3);

k = - a / b;
b = c / b;

N = length(points);

intersection=[];
if (closure ~= 1 && N <= 2)
    tmp = points;
else
    tmp = points;
    tmp(N+1, :) = points(1, :);
end

count_inters = 1;
if (N == 2)
    N = 1;
end
if (N == 3)
    N = 3;
end
for i = 1:N
    A = tmp(i, :);
    B = tmp(i+1, :);
       
    if((k * A(1) + b - A(2)) * (k * B(1) + b - B(2)) > 0)
%         hord_dist = 0;
        continue;
    end
    
    if( (k * A(1) + b - A(2)) * (k * B(1) + b - B(2)) == 0)
        if(k * A(1) + b - A(2) == 0)
            intersection(count_inters, :) = A;
            count_inters = count_inters + 1;
        end

        continue
    end
    koef = abs((k * A(1) + b - A(2))/(k * B(1) + b - B(2)));
    P = [(A(1) + koef * B(1)) / (1 + koef) (A(2) + koef * B(2)) / (1 + koef)];
    intersection(count_inters, :) = P;
    count_inters = count_inters + 1;
end

if(length(intersection) <= 1)

%     hord_dist = 0;

    return 
end
    
% hord_dist = dist(intersection(1), intersection(2));
% if(length(intersection) > 2)
%     ELEMENT.draw(elem, "b", -1, true)
%     disp(length(intersection));
% end
end