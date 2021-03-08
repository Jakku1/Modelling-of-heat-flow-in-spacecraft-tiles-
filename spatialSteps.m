i=0; %set count to 0
nt = 501; % default fumber of timesteps
thick = 0.05; 
tmax = 4000; 
for nx = 2:51 %test 2 to 51 spatial steps
    i=i+1; 
    dx(i) = thick/(nx-1); %finds the spatial step
    disp (['nx = ' num2str(nx) ', dx = ' num2str(dx(i)) ' m']) %outputs number of and size of spatial steps
    %run shuttle code for each method
    % save first temperature value at the inner surface
    [~, ~, u] = shuttle(tmax, nt, thick, nx, 'forward', false);
    uf(i) = u(end, 1); 
    [~, ~, u] = shuttle(tmax, nt, thick, nx, 'backward', false); 
    ub(i) = u(end, 1); 
    [~, ~, u] = shuttle(tmax, nt, thick, nx, 'dufort-frankel', false); 
    ud(i) = u(end, 1); 
    [~, ~, u] = shuttle(tmax, nt, thick, nx, 'crank-nicolson', false); 
    uc(i) = u(end, 1); 
end 
plot(dx, [uf; ub; ud; uc]) %plot size of spatial steps against each temperature value
ylim([350 550]) 
xlabel('spatial step /m')
ylabel('final temperature /K')
%xlim([0 0.03])
legend ('Forward', 'Backward','dufort-frankel','crank-nicolson')
