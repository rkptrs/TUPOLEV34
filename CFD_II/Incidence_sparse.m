clear all
close all
clc

N = 4;
tE21 = zeros(N^2, 2*N*(N+1));
for i = 1:N^2
    tE21(i, i+floor((i-1)/N)) = -1;
    tE21(i, i+1+floor((i-1)/N)) = 1;
    tE21(i, (N+1)*N+i) = -1;
    tE21(i, (N+1)*N+i+N) = 1;
end 
tE21;

n1 = N^2;
n2 = 2*N*(N+1);
d = 0:N;
values = [-ones(N,1); zeros(N^2 - N, 1)];
for i = 1:N-1
    value = [zeros((i-1)*N,1);ones(N,1); -ones(N,1);zeros(N^2-(i+1)*N,1)];
    values = [values, value];
end
values = [values, [zeros(N^2 - N, 1);ones(N,1)]];
sparsetE21_1 = spdiags(values, d, n1, n2/2);
sparsetE21_2 = spdiags([-ones(n1,1), ones(n1,1)], [0 N], n1, n2/2);
sparsetE21 = [sparsetE21_1, sparsetE21_2];
% full(sparsetE21)-tE21;

% transpose(tE21);
% size(tE21);
% tE21 = sparse(tE21);
% 
% for i = 1:N
%     LEFT(i) = 1 + (i-1)*(N+1);
%     RIGHT(i) = LEFT(i) + N;
%     BOTTOM(i) = N*(N+1) + i;
%     TOP(i) = BOTTOM(i) + N^2;
% end
% 
% % wall velocities
% U_wall_top = -1;
% U_wall_bot = 0;
% U_wall_left = 0;
% U_wall_right = 0;
% V_wall_top = 0;
% V_wall_bot = 0;
% V_wall_left = 0;
% V_wall_right = 0;
% 
% u_norm_E = zeros(N^2, 4*N);
% extra = 0;
% index = [];
% for i = 1:N
%     u_norm_E(:, i + extra) = tE21(:, LEFT(i));
%     u_norm_E(:,i+1 + extra) = tE21(:, RIGHT(i));
%     U_norm(i + extra,1) = U_wall_left;
%     U_norm(i+1+extra,1) = U_wall_right;
%     extra = extra + 1;
% end
% u_norm_E(:, [2*N+1:3*N]) = tE21(:, BOTTOM);
% u_norm_E(:, [3*N+1:4*N]) = tE21(:, TOP);
% U_norm([2*N+1:3*N],1) = V_wall_bot;
% U_norm([3*N+1:4*N],1) = V_wall_top;
% 
% tE21(:,[LEFT, RIGHT, TOP, BOTTOM]) = [];
% size(transpose(tE21))
% u_norm = u_norm_E*U_norm;
% 
% E10 = zeros((N-1)*(2*N), N^2);
% for i = 1:(N-1)*N
%     E10(i, i+floor((i-1)/(N-1))) = -1;
%     E10(i, i+1+floor((i-1)/(N-1))) = 1;
% end 
% for i = (N-1)*N+1:(N-1)*2*N
%     E10(i, i-(N-1)*N) = -1;
%     E10(i, i-(N-1)*N+N) = 1;
% end
% 
% % transpose(tE21)+E10





