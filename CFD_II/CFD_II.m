clear all
close all
clc

% This syetm that you need to solve will be singular. Matlab gives you a
% warning at each time step. To switch of this warning, remove the comment
% in the next line

warning on

% 00D#MMXIX#

% This file contains the skeleton of the program which will solve the lid
% driven cavity problem on a unit square. The parts that have to be
% supplemented are described in the assignment.
%
% The pieces that need to be filled in are indicated
%

%
% When running the code, determine a suitable time step. A too small time
% step will make the calculation very long, while for a too large time step
% the solution will blow up due to numerical instability.
%

Re = 1000;              % Reynolds number
N = 3;                 % Number of volumes in the x- and y-direction
Delta = 1/N;            % uniform spacing to be used in the mapping to compute tx


% Determine a suitable time step and stopping criterion, tol

dt = ;             % time step
tol =;             % tol determines when steady state is reached and the program terminates

% wall velocities
U_wall_top = -1;
U_wall_bot = 0;
U_wall_left = 0;
U_wall_right = 0;
V_wall_top = 0;
V_wall_bot = 0;
V_wall_left = 0;
V_wall_right = 0;

%
%   Generation of a non-uniform mesh
%

%
%   tx are the coordinates of the nodal points on the outer-oriented grid
%
tx = zeros(1,N+1);
for i=1:N+1
    xi = (i-1)*Delta;
    tx(i) = 0.5*(1. - cos(pi*xi));
end

% Local mesh size on outer oriented grid
th = zeros(N,1);
th = tx(2:N+1) - tx(1:N);

%
%  x are the coordinates of the nodal points on the inner-oriented grid (including
%  endpoints 0 and 1)
%  h contains the edge lengths on the inner-oriented grid
%
x = 0.5*(tx(1:N) + tx(2:N+1));
x = [0 x 1];

h = zeros(N+1,1);
h = x(2:N+2) - x(1:N+1);

%
%   Initial condition u=v=0
%
%   Both u and v will be stored in one big vector called 'u'
%
%   The vector u only contains the true unknowns, not the velocities that
%   are prescribed by the boundary conditions
%
%   The vector u contains the *inner-oriented* circulations as unknowns

u = zeros(2*N*(N-1),1);

% Set up the Incindence matrix 'tE21' which connects the fluxes to the
% volumes. Use the orientation described in the assignment.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
tE21 = zeros(N^2, 2*N*(N+1));
for i = 1:N^2
    tE21(i, i+floor((i-1)/N)) = -1;
    tE21(i, i+1+floor((i-1)/N)) = 1;
    tE21(i, (N+1)*N+i) = -1;
    tE21(i, (N+1)*N+i+N) = 1;
end 
tE21 = sparse(tE21);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%  Inserting boundary conditions for normal velocity components
%  Multiplication by h components is done in order to obtain integral flux
%  values 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Here you need to modify the incidence matrix tE21 to include the
%    boundary conditions and store the boundary terms un u_norm (see
%    assignment
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setting up simple Hodge matrix which converts the fluxes on the
% outer-oriented grid to circulation on the inner-oriented grid. Assume
% that the fluxes and circulation are constant over each 1-cell. This will
% give a diagonal Hodge matrix. Call this Hodge matrix 'H1t1'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Here you need to construct the incidence matrix H1t1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hu_norm is the vector which will contain the Hodge of the prescribed
% normal fluxes. Calculate the vector 'Hu_norm'

%Vector ubc_norm containing fluxes of prescribed velocities is reused

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Multiply H1t1 with u_norm to het Hu_norm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

%  Remove corresponding row and columns from the Hodge matrix and also
%  remove the corresp0onding 'rows' from Hu_norm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Remove rows and colums in H1t1 and Hu_norm for prescribed values
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up the incidence E^{2,1} between 1-cochain circulation and 2-cochain vorticity on
% the inner-oriented (extended) grid
%
% This incidence matrix will be called 'E21' in the program
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Set up the incidence matrix E21 on the dual grid
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Inserting prescribed tangential bundary conditions

%Vector ubc of size 2(N+2)(N+1) is constructed. It contains all the
%prescribed velocities (normal and tangential) as represented in the inner
%oriented grid.

% Remove columns from the incidence matrix E21 corresponding to both the
% prescribed tangental velocities and normal velocities


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Remove prescribed values from the matrix (just as was done for tE21)
%    and store the known values in a vector u_pres
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

% Set up the Hodge matrix which maps inner-oriented 2-cochains to
% outer-oriented 0-cochains. Assume that the vorticity is constant in the
% inner-oriented 2-cells. This will give a diagonal Hodge matrix. Call this
% Hodhe matrix 'Ht02'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Set up the Hodge matrix which converts integrated values on the dual
%    grid to point values on the primal grid assuming that the physical
%    quantity is constant over the dual surfaces.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up the Hodge matrix which maps inner-oriented 1-cochains to
% outer-oriented 1-cochains. Call this Hodge matrix 'Ht11'. Assume again
% that everything is constant over the 1-cells, which will then yield a
% diagonal Hdoge matrix.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Set up the Hodge matrix which converts integrated values along dual
%    edge to integral values over primal edges
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%
% The prescribed velocties will play a role in the momentum equation
%
u_pres_tvort=Ht02*u_pres; %U_pres to outer oriented 0 form representing contribution of boundary conditions to point wise vorticity
u_pres = H1t1*E21'*Ht02*u_pres; %U_pres to inner oriented 1 forms


% Now all matrices are set up and the time stepping cam start. 'iter' will
% record the number of time steps. This allows you to give output after a
% preselected number of time steps.
%
% 'diff' will be the maximal du/dt or dv/dt. If 'diff' is sufficiently
% small, steady state has been reached. Determine a suitable value for
% 'tol'
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    From here the code enters the time-stepping loop and no new parts in
%    the code need to be inserted. If done correctly, everything should
%    work now.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diff = 1;
iter = 1;


% Set up the matrix for the Poisson equation    

A = -tE21*Ht11*tE21';

% Perform an LU-decomposition for the pressure matrix A

[L,U] = lu(A);

% Abbreviation for some matrix products which are constant in the time loop

VLaplace = H1t1*E21'*Ht02*E21;
DIV = tE21*Ht11;

while diff > tol
        
    %Vector chi is obtained. It corresponds with the point-wise vorticity
    %at each cell
    chi=Ht02*E21*u+u_pres_tvort;  
    
    %Vectors uxchi and uychi correspond with the multiplications of
    %chi with the horizontal and vertical velocity components at each cell.
    %Only the cells required for vector convective are calculated. The
    %ordering of each vector with respect to the ordering of cells in the
    %grid is different (left to right for uxchi and bottom to top for
    %uychi)
    
    uxchi=zeros((N+1)*(N-1),1);
    uychi=zeros((N+1)*(N-1),1);
    
    for i=1:N-1
        for j=1:N+1
            k=j+(i-1)*(N+1); %Number of vector component
            if j==1                
                uxchi(k)=U_wall_left*chi(j+i*(N+1));                
                uychi(k)=V_wall_bot*chi((i+1)+(j-1)*(N+1));
            elseif j==N+1
                uxchi(k)=U_wall_right*chi(j+i*(N+1));
                uychi(k)=V_wall_top*chi((i+1)+(j-1)*(N+1));
            else
                uxchi(k)=(u(j-1+(i-1)*(N-1))+u(j-1+i*(N-1)))/(2*h(j))*chi(j+i*(N+1));
                uychi(k)=(u(N*(N-1)+i+(j-2)*N)+u(N*(N-1)+i+1+(j-2)*N))/(2*h(j))*chi((i+1)+(j-1)*(N+1));
            end
        end
    end
    
    %Given vectors uxchi and uychi, vector convective can be constructed
    convective=zeros(2*N*(N-1),1);
    for i=1:N
        for j=1:N-1
            convective(j+(i-1)*(N-1))=-h(j+1)/2*(uychi(i+(j-1)*(N+1))+uychi(i+1+(j-1)*(N+1))); %Components along horizontal cell lines
            convective(N*(N-1)+i+(j-1)*N)=h(j+1)/2*(uxchi(i+(j-1)*(N+1))+uxchi(i+1+(j-1)*(N+1))); %Components along vertical cell lines
        end
    end
    
    % Set up the right hand side for the Poisson equation for the pressure
    
    rhs_Poisson  =   DIV*(u/dt  - convective - VLaplace*u/Re - u_pres/Re) + u_norm/dt; 
    
    % Solve for the new pressure
    
    temp = L\rhs_Poisson;
    p = U\temp;
    
    % Store the velocity from the previous time step in the vector u_old
    
    uold = u;
    
    % Udate the velocity field
    
    u = u - dt* (convective - tE21'*p + VLaplace*u/Re + u_pres/Re); 
    
    %
    %  Every other 1000 iterations check whether you approach steady state
    %  and check whether yopu satisfy conservation of mass. The largest
    %  rate at which mass is destroyed or created is denoted by 'maxdiv'.
    %  This number should be very small, in the order of machine precision.
    
    if mod(iter,1000) == 0
    
        maxdiv = max(DIV*u + u_norm) 
        
        diff = max(abs(u-uold))/dt
        
    end
    iter = iter + 1;
end
