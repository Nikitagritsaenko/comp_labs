%%%%% LAB4
clear all

plot_b = 1;
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


[flux,RBDRY,ZBDRY,NBDRY,R,Z,time,rdim,zdim] = gfile_extractor_1t(037000,00156,65);
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
K = zeros(256, 208);
%for j = 0:15
j = 10
    spd_xy = spd_xy - (j-1)*2.5*spd_vect * spd_xy_step(mod(j, 2) + 1);
    line = find_line_eq(spd_xy, aperture_xy);
    %       расстояние до плоскости сечения
    spd_r = -sqrt(norm(spd_xy)^2 - line(3)^2);
    aperture_xz_offset = pdist([spd_xy(1), spd_xy(2); aperture_xy(1), aperture_xy(2)]);
    aperture_xz = [1, 0] * aperture_xz_offset;
    H = spd_r;
    [separ_yj, separ_zj, ind_separj, c_rj, c_zj] = get_cut(RBDRY, ZBDRY, NBDRY, c_r, c_z, H);
    if (plot_b == 1)
        figure();
        hold on;
        grid on;
        axis( [ -0.8, 0.7, -0.6, 0.6 ] );
        title(line(3));
        plot(separ_yj, separ_zj, 'bo');
    end
    %spd_xy = spd_xy - 2.5*spd_vect * spd_xy_step(mod(j, 2) + 1);

    N = 13; % число хорд в секторе
    M = 8; % число разбиений средними поверхностями
    EXTR_N = 2; % число экстремумов
    sector_ids = [];
    sector_idxs = zeros(EXTR_N, N);
    % ========= if j >= 6 => 2 центра else (0;0)
    k = 1;
    if (j < 6)
        step = size(separ_yj, 2)/26;
        start_index = 1;
        index = start_index;
        end_index = size(separ_yj, 2);
        while index < end_index
            sector_ids(k) = round(index);
            index = index + step;
            k = k + 1;
        end  
        cr_j = 0;
        cz_j = 0;
    else
       step = ind_separj/26;
       start_index = ind_separj;
       index = start_index;
       end_index = size(separ_yj, 2);
        while index < end_index
            sector_ids(k) = round(index);
            index = index + step;
            k = k + 1;
        end     
        cr_j = c_rj(1);
        cz_j = c_zj(1);
    end
   % sector_ids = [1 3 6 9 12 15 18 21 24 27 30 32 34 36 39 42 45 48 51 54 58 60 63 66 69 72];
    c_idx = 1;
    num_sec = 0;
    segments_tmp = zeros(M, (EXTR_N) * N, 2);
    for k = 1:1:8
        cur_start_sec = num_sec + 1;
        for i = 1 : EXTR_N * N - 1
            num_sec = num_sec + 1;
            num_point = 0;
            alpha = (k-1) * 0.1;
            mul = 0.033;
            line_points = 2;
            if (k == 8)
                line_points = 8;
            end
            for ii = sector_ids(i):sector_ids(i+1)
                num_point = num_point + 1;
                segments(num_sec, num_point, :) = (1 - alpha) * [separ_yj(ii) separ_zj(ii)] + alpha * [cr_j cz_j];
            end
            for l = 1 : line_points
                alpha = (k-1) * 0.1 + l * mul;
                num_point = num_point + 1;
                segments(num_sec, num_point, :) = (1 - alpha) * [separ_yj(ii) separ_zj(ii)] + alpha * [cr_j cz_j];
            end
            if (k ~= 8)
                alpha = k * 0.1;
                for ii = sector_ids(i+1): -1 :sector_ids(i)
                    num_point = num_point + 1;
                    segments(num_sec, num_point, :) = (1 - alpha) * [separ_yj(ii) separ_zj(ii)] + alpha * [cr_j cz_j];
                end
            else
                num_point = num_point + 1;
                segments(num_sec, num_point, :) = [cr_j cz_j];
            end
            for l = 1 : line_points
                alpha = k * 0.1 - l * mul;
                num_point = num_point + 1;
                segments(num_sec, num_point, :) = (1 - alpha) * [separ_yj(ii) separ_zj(ii)] + alpha * [cr_j cz_j];
            end
        end
        num_sec = num_sec + 1;
        num_point = 0;
        num_point = num_point + 1;
        alpha = (k-1) * 0.1;
        segments(num_sec, num_point, :) = (1 - alpha) * [separ_yj(start_index) separ_zj(start_index)] + alpha * [cr_j cz_j];
        num_point = num_point + 1;
        segments(num_sec, num_point, :) = (1 - alpha) * [separ_yj(end_index) separ_zj(end_index)] + alpha * [cr_j cz_j];
        alpha = k * 0.1;
        num_point = num_point + 1;
        segments(num_sec, num_point, :) = (1 - alpha) * [separ_yj(end_index) separ_zj(end_index)] + alpha * [cr_j cz_j];
        num_point = num_point + 1;
        segments(num_sec, num_point, :) = (1 - alpha) * [separ_yj(start_index) separ_zj(start_index)] + alpha * [cr_j cz_j];
    end
    
   
    % цикл по сегментам разбиения
    
    for i = 1:length(segments)
        set = find_section(segments(i,:,:), line(3));
        
        
        if (~isempty(set))
            if (plot_b == 1)
                plot(set(:, 1), set(:, 2), 'g');
                text(min(set(:, 1)) + (max(set(:, 1)) - min(set(:, 1)))/ 2,min(set(:, 2)) + (max(set(:, 2)) - min(set(:, 2)))/ 2, num2str(i));
                plot(-set(:, 1), set(:, 2), 'g');
            end
            % ----- цикл по лучам
            for spd_z = 1:16
                line2 = find_line_eq([spd_r, spd_z_start + spd_z_step * (spd_z - 1)], [spd_r, 0] + aperture_xz);
                %pdist([spd_r, spd_z_start + spd_z_step * (spd_z - 1); [spd_r, 0] + aperture_xz])
                
                % лучи в 4-х плоскостях попадают в центральный столб токамака
                if (j < 13)
                    points = findIntersection(set, line2);
                    % пересечение прямой и кривой
                    if (size(points, 1) > 1)
                        if (plot_b == 1)
                   %         plot(points(:, 1), points(:, 2), 'ro');
                        end
                        for point_ind = 1:1
                            K(j * 16 + spd_z, i) = K(j * 16 + spd_z, i) + pdist([points(point_ind, 1) points(point_ind, 2); points(point_ind+1, 1) points(point_ind+1, 2)]);
                        end
                    end
                end
                if (j >= 6) % мб и не сработает
                    points = findIntersection([-set(:, 1), set(:, 2)], line2);
                    if (size(points, 1) > 1)
                        if (plot_b == 1)
                     %       plot(points(:, 1), points(:, 2), 'ro');
                        end
                        for point_ind = 1:1
                            K(j * 16 + spd_z, i) = K(j * 16 + spd_z, i) + pdist([points(point_ind, 1) points(point_ind, 2); points(point_ind+1, 1) points(point_ind+1, 2)]);
                        end
                    end
                end
                if (plot_b == 1)
                    plot([spd_r spd_r + (spd_r+aperture_xz(1)-spd_r)*100], [(spd_z_start+spd_z_step*(spd_z-1)),...
                        (spd_z_start+spd_z_step*(spd_z-1))+(0-(spd_z_start+spd_z_step*(spd_z-1)))*100 ], 'k');
                end
            end
            
        end
    end
    
    spd_xy = spd_xy - 2.5*spd_vect * spd_xy_step(mod(j, 2) + 1);
%end
K
