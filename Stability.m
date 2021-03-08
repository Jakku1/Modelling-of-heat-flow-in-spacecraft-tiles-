i=0; %counter
nx = 21; %number of spatial steps 
thick = 0.05; %thickness approximation
tmax = 4000; %maximum time
for nt = 41:20:1001 % run each method for timesteps in increments of 20
    i=i+1; %increase count
    dt(i) = tmax/(nt-1); %finds the size of the timestep in seconds
    disp (['nt = ' num2str(nt) ', dt = ' num2str(dt(i)) ' s']) 
    [~, ~, u] = shuttle(tmax, nt, thick, nx, 'forward', false); % run method with timestep
    uf(i) = u(end, 1); %save the final temperature calculated to the vector
    [~, ~, u] = shuttle(tmax, nt, thick, nx, 'backward', false); 
    ub(i) = u(end, 1); 
    [~, ~, u] = shuttle(tmax, nt, thick, nx, 'dufort-frankel', false); 
    ud(i) = u(end, 1); 
    [~, ~, u] = shuttle(tmax, nt, thick, nx, 'crank-nicolson', false); 
    uc(i) = u(end, 1); 
end 
plot(dt, [uf; ub; ud; uc]) %plot all 
ylim([350 550]) 
xlabel('time step /s')
ylabel('final temperature /K')
legend ('Forward', 'Backward','dufort-frankel','crank-nicolson')

