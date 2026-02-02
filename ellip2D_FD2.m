function [A, F] = ellip2D_FD2(f, g, a, b, c, d, Nx, Ny)
    % Grid setup
    x = linspace(a, b, Nx);
    y = linspace(c, d, Ny);
    hx = x(2) - x(1);
    hy = y(2) - y(1);
    [X, Y] = meshgrid(x, y);

    % Total number of unknowns
    N = Nx * Ny;
    A = sparse(N, N);
    F = zeros(N, 1);

    % Interior points
    for j = 2:Ny-1
        for i = 2:Nx-1
            idx = sub2ind([Ny, Nx], j, i);
            A(idx, idx) = -2/hx^2 - 2/hy^2 + Y(j,i);  % center
            A(idx, sub2ind([Ny, Nx], j, i+1)) = 1/hx^2;  % right
            A(idx, sub2ind([Ny, Nx], j, i-1)) = 1/hx^2;  % left
            A(idx, sub2ind([Ny, Nx], j+1, i)) = 1/hy^2;  % up
            A(idx, sub2ind([Ny, Nx], j-1, i)) = 1/hy^2;  % down
            F(idx) = f(X(j,i), Y(j,i));
        end
    end

    % Dirichlet boundary conditions
    for j = 1:Ny
        for i = [1, Nx]
            idx = sub2ind([Ny, Nx], j, i);
            A(idx, :) = 0;
            A(idx, idx) = 1;
            F(idx) = g(X(j,i), Y(j,i));
        end
    end

    for i = 1:Nx
        for j = [1, Ny]
            idx = sub2ind([Ny, Nx], j, i);
            A(idx, :) = 0;
            A(idx, idx) = 1;
            F(idx) = g(X(j,i), Y(j,i));
        end
    end
end
