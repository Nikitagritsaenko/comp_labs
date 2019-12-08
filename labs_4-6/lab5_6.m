%% LAB 5
clearvars -except K
clc
%%
% ��������� ��� ������ 4 � �������� ������� K
A = K;
cond(A)
figure('color', 'white');
title('L-matrix values');
hold on;
imagesc(A)
colorbar;
%% ������������ ����������� �����
t = eig(A * A');
figure('color', 'white');
title('Singular values');
hold on;
hist(t);
xlabel('lamda');
ylabel('values');
%% ������������ ���������������
k = 0.1;
infA = A * (1-k);
supA = A * (1+k);

for i = 10:10:100
    i
    HeurMinCond(infA, supA, i)
end

%%

% ������� ������ � ��������� (������ b)
time = 150;
frame = 30401;
file_name = '35685_SPD16x16.mat';

data = load(file_name);

S = data.sign_bb(:,:,:);

% 3 ������ ���� ����� �� ������ ����� ����� ����� � ������� S ��� ���������
% �������, �� ������ ������ ����� �� ����� ����, ������� ���� � ���������� ����������
tp = cell2mat(data.Data(1,2)) * 1e-3; 	% time resolution, ms
tz=cell2mat(data.Data(2,2)); 			% delay, ms
t0 = tz;       						% ��� ��� � ��
%% ������������ ���������
figure('color', 'white');
title('shot = 35685; time = 150');
hold on;
B = S(:, :, frame);
imagesc(B);
colorbar;
%% ����� 0 - ���
figure('color', 'white');
title('Solution histogram (method 0)');
hold on;
B = S(:, :, frame);
b = reshape(B, [256, 1]);
x = A \ b;
histogram(x);
xlabel('x');
%% ����� 1; x = inv(A' * A) * A' * b;
% ������-�� ������ ������ ����������
figure('color', 'white');
title('Solution histogram (method 1)');
hold on;
B = S(:, :, frame);
b = reshape(B, [256, 1]);
%x = (A' * b)' / (A' * A);
x = inv(A' * A) * A' * b;

histogram(x);
xlabel('x');

%% �����

%% �������� ����, ��� tolmax = 0
k = 0.1
infA = A * (1-k);
supA = A * (1+k);

B = S(:, :, frame);
bb = reshape(B, [256, 1]);

infb = bb * 0.99;
supb = bb * 1.01;

B = S(:, :, frame);
bb = reshape(B, [256, 1]);

[tolmax, argmax, envs, ccode, flag] = tolsolvty(infA, supA, infb, supb);

tolmax
infb = infb - abs(tolmax) - 0.000002;
supb = supb + abs(tolmax) + 0.000002;

[tolmax, argmax, envs, ccode, flag] = tolsolvty(infA, supA, infb, supb);
if flag == 1
    figure('color', 'white');
    histogram(argmax);
    title("Solution histogram (method 2), tolmax = " + tolmax);
    xlabel('x');
end
    

%% ����� 3 (�������� ����������������)
N = 256;
M = 208;

rad = 1;

B = S(:, :, frame);
b = reshape(B, [N, 1]);
e = ones(N + M, 1);
for i = 1:M
    e(i) = 0;
end
diag_rad_b = zeros(N, N);
for i = 1:N
    for j = 1:N
        if (i == j)
            diag_rad_b(i, j) = rad;
        end
    end
end

C = [A -diag_rad_b
    -A -diag_rad_b];
b = double(b);
d = [b
    -b];
lb = zeros(N+M, 1);

xw = linprog(e,C,d,[],[],lb);
x = xw(1:M, :);

figure('color', 'white');
histogram(x);
title("Solution histogram (method 3)");
xlabel('x');

figure('color', 'white');
title("x-values");
hold on;
imagesc(x);
colorbar;
        
        