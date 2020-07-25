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

%%
s = zeros(M,1);
for i = 1:M
    s(i) = mean(data(i,:));
end

ss = zeros(M,M);
for i = 1:M
    for j = i+1:M
        ss(i,j) = mean(data(i,:).*data(j,:));
    end
end

%%
fid = fopen('s.txt', 'w');
fprintf(fid, '%.10f\n', s);
fclose(fid);

fid = fopen('ss.txt', 'w');
for i = 1:M
    for j = 1:M
        if i<j
            fprintf(fid, '%.10f ', ss(i,j));
        elseif i>j
            fprintf(fid, '%.10f ', ss(j,i));
        else
            fprintf(fid, '0 ');
        end
    end
    fprintf(fid, '\n');
end
fclose(fid);
