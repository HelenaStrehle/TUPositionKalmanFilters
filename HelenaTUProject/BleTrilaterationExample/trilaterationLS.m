function position = trilaterationLS(beacons, distances)
    % Input
    %   beacons: Nx2 matrix of beacon positions
    %   rssi distances
    % output
    %   position
    
    A = [];
    B = [];
    
    x0 = beacons(1,1); 
    y0 = beacons(1,2); 
    d0 = distances(1);
    
    for i = 2:size(beacons,1)
        xi = beacons(i,1);
        yi = beacons(i,2);
        di = distances(i);
        
        A = [A; 2*(xi - x0), 2*(yi - y0)];
        B = [B; d0^2 - di^2 - x0^2 + xi^2 - y0^2 + yi^2];
    end
    
    position = A \ B; % Solve Ax = B using least squares
end
