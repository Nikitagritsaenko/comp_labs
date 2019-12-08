function [points] = findIntersection(curve, line)
eps = 1e-11;
points = [];
pointsInd = 1;
%     length(curve(:, 1))
% print()
if ~isempty(curve)
    for i = 1:length(curve(:, 1)) - 1
        curveLine = find_line_eq([curve(i, 1), curve(i, 2)], [curve(i + 1, 1), curve(i + 1, 2)]);
        p = linsolve([line(1) line(2); curveLine(1) curveLine(2)], [line(3); curveLine(3)]);
        
        minx = min(curve(i, 1), curve(i + 1, 1));
        maxx = max(curve(i, 1), curve(i + 1, 1));
        miny = min(curve(i, 2), curve(i + 1, 2));
        maxy = max(curve(i, 2), curve(i + 1, 2));
        
        if ((minx - eps <= p(1, 1)) && (p(1, 1) <= maxx + eps)  && (miny - eps <= p(2, 1)) && (p(2, 1) <= maxy + eps))
            points(pointsInd, 1:4) = [p(1, 1), p(2, 1), i, i + 1];
            pointsInd = pointsInd + 1;
       end
    end
    
    if (size(points, 1) > 1)
        points = uniquetol(points, eps, 'ByRows', true);
        
        [~,idx] = sort(points(:,1));
        points = points(idx, :);
    end
end
end