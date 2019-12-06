%% LAB 5
clearvars -except K
clc
%%
% ��������� ��� ������ 4 � �������� ������� K
A = K;
cond(A);
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
x = b \ A;
histogram(x);
xlabel('x');
%% ����� 1; x = inv(A' * A) * A' * b;
% ������-�� ������ ������ ����������
figure('color', 'white');
title('Solution histogram (method 1)');
hold on;
B = S(:, :, frame);
b = reshape(B, [256, 1]);
x = (A' * b)' / (A' * A);
histogram(x);
xlabel('x');

%% �����
%% ����� 2-1 - ������� ������ ������� ���� 3*3, 5*5, ..., � ��� ������ ��� � ����
ws = 1; % ������ ���� (��������)
flag = 0; % ������� �� ���������� ���������
while (1)
    d = floor(ws / 2);
    B = S(:, :, frame);
    bb = reshape(B, [256, 1]);
    b = [bb bb];
    N = size(B, 1); M = size(B, 2);
    
    for i = 1:N
        
        for j = 1:M
            
            max = b((i-1) * N + j); min = b((i-1) * N + j);
            x1 = d; x2 = d; y1 = d; y2 = d;
            while i - x1 < 1
                x1 = x1 - 1;
            end
            while j - y1 < 1
                y1 = y1 - 1;
            end
            while i + x2 > N
                x2 = x2 - 1;
            end
            while j + y2 > M
                y2 = y2 - 1;
            end
            for k = i-x1:1:i+x2
                for m = j-y1:1:j+y2
                    if B(k, m) > max
                        max = B(k, m);
                    end
                    if B(k, m) < min
                        min = B(k, m);
                    end
                end
            end
            
            b((i-1) * N + j, :) = [min max];
        end
    end
    
    k = 0; % ���������� � ������� A (����������� ��������� ���� ���� ��� ��)
    infA = A * (1-k);
    supA = A * (1+k);
    infb = b(:,1);
    supb =  b(:,2);
    sprintf('����������� ��� ������� ���� = %i; ������ ���� = %i', k, ws)
    ws = ws + 2;
    [tolmax, argmax, envs, ccode, flag] = tolsolvty(infA, supA, infb, supb);
    
    if ws > 5
        if flag == 0
            disp('������� ������� ����, � ������� ��� � �� �������')
        end
        break;
    end
    
    if flag == 1
        figure('color', 'white');
        histogram(argmax);
        title("Solution histogram (method 2-1), window (" + ws + "*" + ws + "), "+ "frame = " + frame);
        xlabel('x');
    end
end

%% ����� 2-2 - ������� ��������� ����, � ��� ������ ��� � ���� �������������.
df = 1; % ������ �������� ���������� ����
F = frame;
flag = 0;
while flag == 0
    B = S(:, :, frame);
    bb = reshape(B, [256, 1]);
    b = [bb bb];
    for frame_idx = F-df:1:F+df
        B = S(:, :, frame_idx);
        bb1 = reshape(B, [256, 1]);
        b1 = [bb1 bb1];
        
        % ������� ������ b
        for q = 1:size(b,1)
            if (b1(q, 1) < b(q, 1))
                b(q, 1) = b1(q, 1);
            end
            if (b1(q, 2) > b(q, 2))
                b(q, 2) = b1(q, 2);
            end
        end
    end
    k = 0.1; % ���������� � ������� A (����������� ��������� ���� ���� ��� ��)
    infA = A * (1-k);
    supA = A * (1+k);
    infb = b(:,1);
    supb =  b(:,2);
    
    sprintf('����������� ��� ������� ���� = %i; ����� [%i ... %i]', k, F - df, F + df)
    
    [tolmax, argmax, envs, ccode, flag] = tolsolvty(infA, supA, infb, supb);
    df = df + 1;
    
    if df > 10 && flag == 0
        disp('������� ������� ��������� ����������, � ������� ��� � �� �������')
        break;
    end
    
    if flag == 1
        figure('color', 'white');
        histogram(argmax);
        l = (F - df); r = (F + df);
        title("Solution histogram (method 2-2), window (" + ws + "*" + ws + "frames = [" + l + "; " + r + "]");
        xlabel('x');
    end
end

%% ����� 2-3 (���������������) - ������� ��������� ����, � ��� ������ ��� � ���� �������������. ���� �� ������� 3*3 
df = 1;
ws = 3;
d = floor(ws / 2);
F = frame;
flag = 0;
while (1)
    B = S(:, :, frame);
    bb = reshape(B, [256, 1]);
    b = [bb bb];
    for frame_idx = F-df:1:F+df
        B = S(:, :, frame_idx);
        bb1 = reshape(B, [256, 1]);
        b1 = [bb1 bb1];
        
        N = size(B, 1); M = size(B, 2);
        
        for i = 1:N
            
            for j = 1:M               
                max = b1((i-1) * N + j); min = b1((i-1) * N + j);
                
                x1 = d; x2 = d; y1 = d; y2 = d;
                if i - x1 < 1
                    x1 = 0;
                end
                if j - y1 < 1
                    y1 = 0;
                end
                if i + x2 > N
                    x2 = 0;
                end
                if j + y2 > M
                    y2 = 0;
                end
                for k = i-x1:1:i+x2
                    for m = j-y1:1:j+y2
                        if B(k, m) > max
                            max = B(k, m);
                        end
                        if B(k, m) < min
                            min = B(k, m);
                        end
                    end
                end
                
                b1((i-1) * N + j, :) = [min max];
            end
        end
        
        % ������� ������ b
        for q = 1:size(b,1)
            if (b1(q, 1) < b(q, 1))
                b(q, 1) = b1(q, 1);
            end
            if (b1(q, 2) > b(q, 2))
                b(q, 2) = b1(q, 2);
            end
        end
    end
    k = 0; % ���������� � ������� A (����������� ��������� ���� ���� ��� ��)
    infA = A * (1-k);
    supA = A * (1+k);
    infb = b(:,1);
    supb =  b(:,2);
    
    sprintf('����������� ��� ������� ���� = %i; ���� �� ������� 3 * 3 ;����� [%i ... %i]', k, F - df, F + df)
    
    [tolmax, argmax, envs, ccode, flag] = tolsolvty(infA, supA, infb, supb);
    df = df + 1;
    
    if df > 10
        break;
    end  
    
    if flag == 1
        figure('color', 'white');
        histogram(argmax);
        l = (F - df); r = (F + df);
        title("Solution histogram (method 2-3), window 3*3, frames = [" + l + "; " + r + "]");
        xlabel('x');
    end
end

