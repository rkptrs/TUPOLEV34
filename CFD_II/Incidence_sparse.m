clear all
close all
clc

N = 50;
tic;
for iter = 1:20
tE21 = sparse([]);
for i = 1:N^2
    tE21(i, i+floor((i-1)/N)) = -1;
    tE21(i, i+1+floor((i-1)/N)) = 1;
    tE21(i, (N+1)*N+i) = -1;
    tE21(i, (N+1)*N+i+N) = 1;
end 
% sE21 = sparse(tE21);
% tE21;


% tic;
% d = 0:N;
% values = [-ones(N,1); zeros(N^2 - N, 1)];
% for i = 1:N-1
%     value = [zeros((i-1)*N,1);ones(N,1); -ones(N,1);zeros(N^2-(i+1)*N,1)];
%     values = [values, value];
% end
% values = [values, [zeros(N^2 - N, 1);ones(N,1)]];
% sparsetE21_1 = spdiags(values, d, N^2, N*(N+1));
% sparsetE21_2 = spdiags([-ones(N^2,1), ones(N^2,1)], [0 N], N^2, N*(N+1));
% sparsetE21 = [sparsetE21_1, sparsetE21_2];
% toc;
% full(sparsetE21)-tE21;



% transpose(tE21);
% size(tE21);
% tE21 = sparse(tE21);
% 

for i = 1:N
    LEFT(i) = 1 + (i-1)*(N+1);
    RIGHT(i) = LEFT(i) + N;
    BOTTOM(i) = N*(N+1) + i;
    TOP(i) = BOTTOM(i) + N^2;
end

% wall velocities
U_wall_top = -1;
U_wall_bot = 0;
U_wall_left = 0;
U_wall_right = 0;
V_wall_top = 0;
V_wall_bot = 0;
V_wall_left = 0;
V_wall_right = 0;

u_norm_E = sparse([]);
U_norm = sparse([]);
extra = 0;
index = [];
for i = 1:N
    u_norm_E(:, i + extra) = tE21(:, LEFT(i));
    u_norm_E(:,i+1 + extra) = tE21(:, RIGHT(i));
    U_norm(i + extra,1) = U_wall_left;
    U_norm(i+1+extra,1) = U_wall_right;
    extra = extra + 1;
end
u_norm_E(:, [2*N+1:3*N]) = tE21(:, BOTTOM);
u_norm_E(:, [3*N+1:4*N]) = tE21(:, TOP);
U_norm([2*N+1:3*N],1) = V_wall_bot;
U_norm([3*N+1:4*N],1) = V_wall_top;

tE21(:,[LEFT, RIGHT, TOP, BOTTOM]) = [];
u_norm = u_norm_E*U_norm;

E10 = sparse([]);
for i = 1:(N-1)*N
    E10(i, i+floor((i-1)/(N-1))) = -1;
    E10(i, i+1+floor((i-1)/(N-1))) = 1;
end 
for i = (N-1)*N+1:(N-1)*2*N
    E10(i, i-(N-1)*N) = -1;
    E10(i, i-(N-1)*N+N) = 1;
end
end
toc;
% transpose(tE21)+E10





