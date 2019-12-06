%%%%% LAB4
clear all
clc
%%
% угол между направлением камеры-обскуры и направлением на центр (между 8 и 9 лучами)
ang = acos((700^2 + 720^2 - 31^2) / (2 * 700 * 720));
% положение края детектора (1-го столбца)
spd_start = [0, -0.708];
% положение 16-го столбца
spd_end = [0.72 * sin(ang), 0.72 * -cos(ang)];
% вектор направления камеры-обскуры в экваториальной плоскости
spd_vect = (spd_end - spd_start) / norm(spd_end - spd_start);
% шаг между столбцами в плоскости детектора, 2 числа
spd_xy_step = [2.3375 - 0.88 , 3.81 - 2.3375 + 0.88 ] * 1e-03;
%  центр детектора
pp = spd_start + spd_vect * ((spd_xy_step(1) + spd_xy_step(2)) * 8 + 0.52 * 1e-03) / 2;% + spd_vect * 0.35 / 2 * 1e-03;
% отступ апертуры от центра детектора
aperture_xy_offset = 0.0395;
% координата апертуры
aperture_xy = [pp(1) - spd_vect(2) * aperture_xy_offset, pp(2) + spd_vect(1) * aperture_xy_offset];
%spd xz – устройство детектора в меридиональной плоскости
spd_z_start = (27.52 - 0.49) / 2 * 1e-03;
spd_z_step = -1.72 * 1e-03;


spd_xy = spd_start + spd_vect * (spd_xy_step(2) / 2 + 0.26 * 1e-03);


[flux,RBDRY,ZBDRY,NBDRY,R,Z,time,rdim,zdim] = gfile_extractor_1t(035685,00150,65);
c_r = 0
c_z = 0
r_size = size(R, 2);
z_size = size(Z, 2);

% plot custom grid
r_min = min(R)-1
r_max = max(R)+1
z_min = min(Z)-1
z_max = max(Z)+1
r_step = rdim / r_size
z_step = zdim / z_size
min_flux = 1000000
for i = 1:r_size
    for j = 1:z_size
        if (min_flux > flux(i, j))
            min_flux = flux(i, j)
            c_r = R(j);
            c_z = Z(i);
            %       hold on
        end
    end
end

plot_b = 1;
K = zeros(256, 208);
points = zeros(208, 4, 2);%208 sectors of 4x[360 points with equal (r, z)] => 4 points with (r,z)


N = 13; % число хорд в секторе
M = 8; % число разбиений средними поверхностями
EXTR_N = 2; % число экстремумов
sector_ids = [];
% ========= if j >= 6 => 2 центра else (0;0)
k = 1;

step = NBDRY/26;
start_index = 1;
index = start_index;
end_index = NBDRY;
while index < end_index
    sector_ids(k) = round(index);
    index = index + step;
    k = k + 1;
end

%%
% sector_ids(25) = sector_ids(25) + 2;
% sector_ids(26) = sector_ids(1);
%sector_ids = [1 3 6 9 12 15 18 21 24 27 30 32 34 36 39 42 45 48 51 54 58 60 63 66 69 72];
num_sec = 0;

step = 1/8
for k = 1:1:8
    for i = 1 : EXTR_N * N - 1
        num_sec = num_sec + 1;
        num_point = 0;
        alpha = (k-1) * step;
        i1 = sector_ids(i);
        i2 = sector_ids(i+1);
        points(num_sec, 1, :) = (1 - alpha) * [RBDRY(i1) ZBDRY(i1)] + alpha * [c_r c_z];
        points(num_sec, 2, :) = (1 - alpha) * [RBDRY(i2) ZBDRY(i2)] + alpha * [c_r c_z];
        alpha = (k) * step;
        points(num_sec, 3, :) = (1 - alpha) * [RBDRY(i2) ZBDRY(i2)] + alpha * [c_r c_z];
        points(num_sec, 4, :) = (1 - alpha) * [RBDRY(i1) ZBDRY(i1)] + alpha * [c_r c_z];
        if (k == 8)
            points(num_sec, 3, :) = [c_r c_z];
            points(num_sec, 4, :) = [c_r c_z];
        end
    end
    
    alpha = (k-1) * step;
    num_sec = num_sec + 1;
    points(num_sec, 1, :) = (1 - alpha) * [RBDRY(sector_ids(26)) ZBDRY(sector_ids(26))] + alpha * [c_r c_z];
    points(num_sec, 2, :) = (1 - alpha) * [RBDRY(1) ZBDRY(1)] + alpha * [c_r c_z];
    alpha = (k) * step;
    points(num_sec, 3, :) = (1 - alpha) * [RBDRY(1) ZBDRY(1)] + alpha * [c_r c_z];
    points(num_sec, 4, :) = (1 - alpha) * [RBDRY(sector_ids(26)) ZBDRY(sector_ids(26))] + alpha * [c_r c_z];
    if (k == 8)
        points(num_sec, 3, :) = [c_r c_z];
        points(num_sec, 4, :) = [c_r c_z];
    end
end

LINE = zeros(16, 3);
SPD_R = zeros(16, 1);
SPD_XY = zeros(17, 2);
SPD_XY(1, :) = spd_xy;


for j = 0:15
    line = find_line_eq(SPD_XY(j + 1, :), aperture_xy);
    LINE(j + 1, :) = line;
    SPD_R(j + 1, :) = -sqrt(norm(SPD_XY(j + 1, :))^2 - line(3)^2);
    SPD_XY(j + 2, :) = SPD_XY(j + 1, :) + spd_vect * spd_xy_step(mod(j, 2) + 1);
end
%%
for j = 0:15
    line = LINE(j + 1, :);
    spd_r = SPD_R(j + 1, :);
    spd_xy = SPD_XY(j + 1, :);
    
    figure(); hold on; grid on;
    axis( [ -0.8, 0.7, -0.6, 0.6 ] );
    title(line(3));
    
    aperture_xz_offset = pdist([spd_xy(1), spd_xy(2); aperture_xy(1), aperture_xy(2)]);
    aperture_xz = [1, 0] * aperture_xz_offset;
    
    % цикл по сегментам разбиения
    segments = points;
    
    
    for i = 1:length(segments)
        set = find_section(segments(i,:,:), line(3));
        set_addit = [];
        index_addit = 1;
        if (~isempty(set))
            if (plot_b == 1)
                plot(set(:, 1), set(:, 2), 'b', 'linewidth', 0.8);
                plot(-set(:, 1), set(:, 2), 'b', 'linewidth', 0.8);
                if (size(set, 1) == 3)
                    x_max = max(set(:,1));
                    for q = 1:size(set(:, 1), 1)
                        x = [set(q, 1), -set(q, 1)];
                        y = [set(q, 2), set(q, 2)];
                        if x(1) ~= x_max
                            set_addit(index_addit, :) = [x(1) x(2)];
                            set_addit(index_addit+1, :) = [y(1) y(2)];
                            index_addit = index_addit + 2;
                            plot(x, y, 'b', 'linewidth', 0.8);
                        end
                    end
                elseif (size(set, 1) < 3)
                    for q = 1:size(set(:, 1), 1)
                        x = [set(q, 1), -set(q, 1)];
                        y = [set(q, 2), set(q, 2)];
                        set_addit(index_addit, :) = [x(1) x(2)];
                        set_addit(index_addit+1, :) = [y(1) y(2)];
                        index_addit = index_addit + 2;
                        plot(x, y, 'b', 'linewidth', 0.8);
                    end
                    
                end
            end
            
            % ----- цикл по лучам
            for spd_z = 1:16
                line2 = find_line_eq([spd_r, spd_z_start + spd_z_step * (spd_z - 1)], [spd_r, 0] + aperture_xz);
                
                % лучи в 4-х плоскостях попадают в центральный столб токамака
                if (j < 13)
                    points_int = lineIntersection(set, line2, 1);
                    
                    % пересечение прямой и кривой
                    if (size(points_int, 1) > 1)
                        if (plot_b == 1)
                            plot(points_int(:, 1), points_int(:, 2), 'ko', 'MarkerSize', 3);
                        end
                        for point_ind = 1:1
                            K(j * 16 + spd_z, i) = K(j * 16 + spd_z, i) + pdist([points_int(point_ind, 1) points_int(point_ind, 2); points_int(point_ind+1, 1) points_int(point_ind+1, 2)]);
                        end
                    end
                    if (~isempty(set_addit))
                        for s_a_i = 1:2:length(set_addit(:, 1)) - 1
                            set_a(:,1) = set_addit(s_a_i, :);
                            set_a(:,2) = set_addit(s_a_i + 1, :);
                            points_int = lineIntersection(set_a, line2, 0);
                            
                            % пересечение прямой и кривой
                            if (size(points_int, 1) > 0)
                                if (plot_b == 1)
                                    plot(points_int(:, 1), points_int(:, 2), 'ko', 'MarkerSize', 3);
                                end
                                for point_ind = 1:1
%                                     K(j * 16 + spd_z, i) = K(j * 16 + spd_z, i) + pdist([points_int(point_ind, 1) points_int(point_ind, 2); points_int(point_ind+1, 1) points_int(point_ind+1, 2)]);
                                end
                            end
                        end
                    end
                end
                % if (j >= 6) % мб и не сработает
                
                points_int = lineIntersection([-set(:, 1), set(:, 2)], line2, 1);
                
                if (size(points_int, 1) > 1)
                    if (plot_b == 1)
                        plot(points_int(:, 1), points_int(:, 2), 'ko', 'MarkerSize', 3);
                    end
                    for point_ind = 1:1
                        K(j * 16 + spd_z, i) = K(j * 16 + spd_z, i) + pdist([points_int(point_ind, 1) points_int(point_ind, 2); points_int(point_ind+1, 1) points_int(point_ind+1, 2)]);
                    end
                end
                
                
                if (plot_b == 1)
                    plot([spd_r spd_r + (spd_r+aperture_xz(1)-spd_r)*100], [(spd_z_start+spd_z_step*(spd_z-1)),...
                        (spd_z_start+spd_z_step*(spd_z-1))+(0-(spd_z_start+spd_z_step*(spd_z-1)))*100 ], 'r', 'linewidth', 0.8);
                end
            end
            
        end
    end
    
end
