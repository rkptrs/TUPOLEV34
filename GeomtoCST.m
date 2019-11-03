function [Au,Al] = GeomtoCST()
M = 12  %Number of CST-coefficients in design-vector x

%Define optimization parameters
x0 = 0*ones(M,1);     %initial value of design vector x(starting vector for search process)
lb = -1*ones(M,1);    %upper bound vector of x
ub = 1*ones(M,1);     %lower bound vector of x

options=optimset('Display','Iter');


%perform optimization
global data;
Coor = data.coords;
[x,fval,exitflag] = fmincon(@CST_objective,x0,[],[],[],[],lb,ub,[],options);
M_break=M/2
Au=x(1:M_break);
Al=x(1+M_break:end);


end

