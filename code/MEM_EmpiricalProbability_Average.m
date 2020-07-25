clear all; clc;

%%
N = 300; M = 9;
L = 10000;
THRESHOLD = -50;

DATA = zeros(N, N, L);
fid = fopen('Result_CaL70.dat', 'r');
for i = 1:1000
    fread(fid, N*N, 'single');
end
for i = 1:L
    tmp = fread(fid, N*N, 'single');
    DATA(:,:,i) = reshape(tmp, N, N);
end
fclose(fid);

mean(mean(mean(DATA)))

%%
% c = [1 76 151 226 301];
c = [1 101 201 301];
for i = 1:3
    for j = 1:3
        data(i,j,:) = squeeze(mean(mean(DATA(c(i):c(i+1)-1, c(j):c(j+1)-1, :)))) > THRESHOLD;
    end
end
clear DATA;
data = reshape(data, M, L);

%% Calculate empirical state probability
for i = 1:M
    b(i) = 2^(i-1);
end

P = zeros(1,2^M);
for t = 1:L
    V = 1;
    for i = 1:M
        if data(i,t)
            V = V + b(i);
        end
    end
    P(V) = P(V) + 1;
end
P = P / sum(P);

%%
fid = fopen('PN.txt', 'w');
fprintf(fid, '%.10f\n', P);
fclose(fid);