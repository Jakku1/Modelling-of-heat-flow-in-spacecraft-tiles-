function [x, t, u] = shuttleNoPlot(tmax, nt, xmax, nx, method,tileChoice)
% Function for modelling temperature in a space shuttle tile
% D N Johnston  30/1/19
%
% Input arguments:
% tmax   - maximum time
% nt     - number of timesteps
% xmax   - total thickness
% nx     - number of spatial steps
% method - solution method ('forward', 'backward' etc)
% doplot - true to plot graph; false to suppress graph.
%
% Return arguments:
% x      - distance vector
% t      - time vector
% u      - temperature matrix
%
% For example, to perform a  simulation with 501 time steps
%   [x, t, u] = shuttle(4000, 501, 0.05, 21, 'forward', true);
%

% Set tile properties
thermcon = 0.141; % W/(m K)
density  = 351;   % 22 lb/ft^3
specheat = 1259;  % ~0.3 Btu/lb/F at 500F
alpha = thermcon / (density * specheat);

%load data for selected tile from file
load(strcat('temp',tileChoice,'.mat'));

%convert temperature values into kelvin
tempdata=((tempdata-32)*5/9)+273;


% Initialise everything.
dt = tmax / (nt-1);
t = (0:nt-1) * dt;
dx = xmax / (nx-1);
x = (0:nx-1) * dx;
u = zeros(nt, nx);
p = alpha * dt / dx^2;

% set initial conditions to 16C throughout.
% Do this for first two timesteps.
u([1 2], :) = 298;
ivec = 2:nx-1; % set up index vector

%backwards boundary conditions
%L = 0;

% Main timestepping loop.
for n = 2:nt - 1
    
    % RHS boundary condition: outer surface. 
    % Use interpolation to get temperature R at time t(n+1).
    R = interp1(timedata, tempdata, t(n+1), 'linear', 'extrap');
    
    % Select method.
    switch method
        case 'forward'
            
            %u(n+1,1)=L;
            u(n+1,1)=(1-2*p)*u(n,1)+2*p*u(n,2);
            u(n+1,nx)=R;
               
            for i=2:nx-1
                u(n+1,i) =(1-2*p)*u(n,i)+p*(u(n,i-1)+u(n,i+1));
                
            end
        case 'dufort-frankel'
            u(n+1,1)=((1-2*p)*u(n-1,1)+4*p*u(n,2))/(1+2*p);
            u(n+1,nx)=R;
            
            for i=2:nx-1
                u(n+1,i)=((1-2*p)*u(n-1,i)+ 2*p*(u(n,i-1)+u(n,i+1)))/(1+2*p);
            end
        case 'backward'
            %calculate internal values using backward differencing
            
            u(n+1,1)=(u(n,1)+2*p*u(n+1,2))/(1+2*p);
            u(n+1,nx)=R;
                    
            %neumann boundary conditions
            
            b(1)      = 1+2*p;
            c(1)      = -2*p;
            d(1)      = u(n,1);        
                        
            a(ivec)   = -p;
            b(ivec)   = 1 + 2*p;
            c(ivec)   = -p;
            d(ivec)   = u(n, ivec);
            
            a(nx)     = 0;
            b(nx)     = 1;
            d(nx)     = R;
    
            u(n+1,:) = tdm(a,b,c,d);
            
        case 'crank-nicolson'
            u(n+1,1)=((1-p)*u(n,1)+p*u(n,2)+p*u(n+1,2))/(1+p);
            u(n+1,nx)=R;
            
            %neumann boundaries
            
            b(1)      = 1+p;
            c(1)      = -p;
            d(1)      = (1-p)*u(n,1)+p*u(n,2);
            
            a(ivec)   = -p/2;
            b(ivec)   = 1 + p;
            c(ivec)   = -p/2;
            d(ivec)   = (p/2)*u(n,1:nx-2)+(1-p)*u(n,2:nx-1)+(p/2)*u(n,3:nx);
            
            a(nx)     = 0;
            b(nx)     = 1;
            d(nx)     = R;
            
                        
            u(n+1,:) = tdm(a,b,c,d);
         otherwise
            error (['Undefined method: ' method])
            return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tri-diagonal matrix solution 
function x = tdm(a,b,c,d)
n = length(b);

% Eliminate a terms
for i = 2:n
    factor = a(i) / b(i-1);
    b(i) = b(i) - factor * c(i-1);
    d(i) = d(i) - factor * d(i-1);
end

x(n) = d(n) / b(n);

% Loop backwards to find other x values by back-substitution
for i = n-1:-1:1
    x(i) = (d(i) - c(i) * x(i+1)) / b(i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% End of shuttle function



    