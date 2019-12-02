function [points] = findIntersection(curve, line)
    eps = 1e-3;
    points = [];
    pointsInd = 1;

    for i = 1:length(curve(:, 1)) - 1
        eq = find_line_eq([curve(i, 1), curve(i, 2)], [curve(i + 1, 1), curve(i + 1, 2)]);
        a = eq(1); b = eq(2); c = eq(3);
        tx = -100:0.01:100;
        ty = (c - a * tx) / b;
        curveLine = [tx ty]';
        
        for j = 1:size(line, 1)
            p1 = line(j, :);

            for k = 1:size(curveLine, 1)
                p2 = curveLine(k, :);
                dist = norm(p1 - p2);
        
                if dist < eps
                    points(pointsInd, 1:4) = [p2(1, 1), p2(2, 1), i, i + 1];
                    pointsInd = pointsInd + 1;
                end
            end  
        end   
    end

    if (size(points, 1) > 1)
        points = uniquetol(points, eps, 'ByRows', true);
        
        [~,idx] = sort(points(:,1));
        points = points(idx, :);
    end
end

