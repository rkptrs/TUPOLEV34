function [Aur,Alr,Aut,Alt] = GeomtoCST()
M = 12  %Number of CST-coefficients in design-vector x

%Define optimization parameters
x0 = 0*ones(M,1);     %initial value of design vector x(starting vector for search process)
lb = -1*ones(M,1);    %upper bound vector of x
ub = 1*ones(M,1);     %lower bound vector of x

%options=optimset('Display','Iter');
fid= fopen('withcomb135.dat','r'); % Filename can be changed as required
Coor = fscanf(fid,'%g %g',[2 Inf]); 
fclose(fid);
Coor = Coor';
global data;
data.coords = Coor;


%perform optimization

[x,fval,exitflag] = fmincon(@CST_objective,x0,[],[],[],[],lb,ub,[],[]);
M_break=M/2;
Aur=x(1:M_break);
Alr=x(1+M_break:end);

global data;
data.coords = [Coor(:,1), Coor(:,2)*0.6];

[x,fval,exitflag] = fmincon(@CST_objective,x0,[],[],[],[],lb,ub,[],[]);
M_break=M/2;
Aut=x(1:M_break);
Alt=x(1+M_break:end);

end

