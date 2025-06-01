function position = trilateration_explicit(B1, B2, B3, r1, r2, r3)
    % trilateration_explicit - Computes (x, y) position using explicit trilateration formula
    %
    % INPUTS:
    %   Beacons pos
    %   rssi distance
    %
    % OUTPUT:
    %   position = [x, y]

    x1 = B1(1); y1 = B1(2);
    x2 = B2(1); y2 = B2(2);
    x3 = B3(1); y3 = B3(2);

    % Compute distances
    U = x2 - x1; 
    Vx = x3 - x1; 
    Vy = y3 - y1; 

    V2 = Vx^2 + Vy^2;

    x = (r1^2 - r2^2 + U^2) / (2 * U);

    if Vy == 0
        error('Vy is zero, division error in trilateration.');
    end

    % Compute y 
    y = (r1^2 - r3^2 + V2 - 2 * Vx * x) / (2 * Vy);

    %global coordinates
    position = [x + x1, y + y1];
end
