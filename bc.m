%
%  bc.m
%  code3
%
%  Created by Yuan Leiqi (袁磊祺) on 2021/6/11.
%

function [W] = bc(W)
    % set boundary conditions
    %
    global Nx Ny
    d = (Nx * 0.2);
    h = (Ny * 0.2);

    % Entry condition
    for j = 1:Ny
        W(1, j, :) = [1.4, 3, 0, 1];
        W(2, j, :) = [1.4, 3, 0, 1];
        W(3, j, :) = [1.4, 3, 0, 1];
    end

    % Export condition
    for j = 1:Ny
        W(Nx, j, :) = W(Nx - 2, j, :);
        W(Nx - 1, j, :) = W(Nx - 2, j, :);
        W(Nx - 2, j, :) = W(Nx - 2, j, :);
    end

    % Top condition
    for i = 2:Nx - 1
        W(i, Ny, :) = W(i, Ny - 5, :);
        W(i, Ny, 3) = -W(i, Ny - 5, 3);
        W(i, Ny - 1, :) = W(i, Ny - 4, :);
        W(i, Ny - 1, 3) = -W(i, Ny - 4, 3);
        W(i, Ny - 2, :) = W(i, Ny - 3, :);
        W(i, Ny - 2, 3) = -W(i, Ny - 3, 3);
    end

    % Bottom condition
    for i = 2:d - 1
        W(i, 1, :) = W(i, 6, :);
        W(i, 1, 3) = -W(i, 6, 3);
        W(i, 2, :) = W(i, 5, :);
        W(i, 2, 3) = -W(i, 5, 3);
        W(i, 3, :) = W(i, 4, :);
        W(i, 3, 3) = -W(i, 4, 3);
    end

    for i = d + 1:Nx - 1
        W(i, h - 2, :) = W(i, h + 3, :);
        W(i, h - 2, 3) = -W(i, h + 3, 3);
        W(i, h - 1, :) = W(i, h + 2, :);
        W(i, h - 1, 3) = -W(i, h + 2, 3);
        W(i, h, :) = W(i, h + 1, :);
        W(i, h, 3) = -W(i, h + 1, 3);
    end

    for j = 2:h - 1
        W(d + 2, j, :) = W(d - 3, j, :);
        W(d + 2, j, 2) = -W(d - 3, j, 2);
        W(d + 1, j, :) = W(d - 2, j, :);
        W(d + 1, j, 2) = -W(d - 2, j, 2);
        W(d, j, :) = W(d - 1, j, :);
        W(d, j, 2) = -W(d - 1, j, 2);
    end

end
