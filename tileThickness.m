function [x, t, u] = tileThickness(tmax, nt, xmax, nx, method,tileChoice)
    %run shuttle funcion to find the initial temperature approximation
    [x, t, u] = shuttleNoPlot(tmax, nt, xmax, nx, method, tileChoice);
    
    while max(u(:,1))>449.81
        xmax = xmax + 0.001; %increase tile thickness and re attempt thickness calculation
        %new temperature caluclation for new thickness
        [x, t, u] = shuttleNoPlot(tmax, nt, xmax, nx, method, tileChoice);
    % contour plot
    surf(x,t,u)
    % comment out the next line to change the surface appearance
    shading interp 
    
    % Rotate the view
    view(150,30)
    
    %label the axes
    xlabel('\itx\rm - m')
    ylabel('\itt\rm - s')
    zlabel('\itu\rm - K')
    switch method
        case 'forward'
            title('Forward Differencing')
        case'dufort-frankel'
            title('Dufort-Frankel')
        case 'backward' 
            title('Backwards Differencing')
        case 'crank-nicolson'
            title('Crank-Nicolson')
            
    pause(0.1)
    end
    %output required thickness
    xmax
end
