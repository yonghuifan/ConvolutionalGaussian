function drawCircles(X,r,c)
    center = X;
    th = linspace(0, 2*pi);
    x = r.*cos(angle)+ center(1);
    y = r.*sin(angle)+ center(2);
    z = center(3);
    plot(x, y, z,c)
end

