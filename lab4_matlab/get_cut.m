function [result_y, result_z, N_split, magnet_axis_y, magnet_axis_z]=get_cut(RBDRY, ZBDRY, NBDRY, magnet_axis_R, magnet_axis_Z, H)
result_y_left=[];
result_z_left=[];
result_y_right=[];
result_z_right=[];
N_right = 0;
N_left = 0;
%можно передавать длинну массивов, а можно вычислять "по ходу"
%for i = 1:NBDRY
for i = 1:length(RBDRY)
    z = ZBDRY(i);
    r = RBDRY(i);
    y1 = sqrt(r*r-H*H);
    y2 = -sqrt(r*r-H*H);
    
    if(isreal(y1))
        N_right = N_right+1;
        result_y_right(N_right) = y1;
        result_z_right(N_right) = z;
    end
    if(isreal(y2))
        N_left = N_left+1;
        result_y_left(N_left) = y2;
        result_z_left(N_left) = z;

    end     
end

result_y_left = rot90(result_y_left, 2);
result_z_left = rot90(result_z_left, 2);
result_y = [result_y_left, result_y_right];
result_z = [result_z_left, result_z_right];
if( H < min(RBDRY))
    N_split = N_left + 1;
else
    N_split = -1;
end
    
z = magnet_axis_Z;
r = magnet_axis_R;
magnet_axis_y = [];
magnet_axis_z = [];
y1 = sqrt(r*r-H*H);
y2 = -sqrt(r*r-H*H); 
if(isreal(y1))
    magnet_axis_y(1) = y1;
    magnet_axis_z(1) = z;
end
if(isreal(y2))
    magnet_axis_y(2) = y2;
    magnet_axis_z(2) = z;
end     

end