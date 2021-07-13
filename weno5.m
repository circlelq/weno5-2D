%
%  weno5.m
%  weno5-2D
%
%  Created by Yuan Leiqi (袁磊祺) on 2021/7/8.
%

function [flux, dt] = weno5(U)
    % WENO5 scheme
    global Nx Ny CFL dx

    flux = zeros(Nx, Ny, 4); % flux vector

    W = U2W(U);
    f_1 = zeros(Nx, Ny, 4); % 正通量
    f_2 = zeros(Nx, Ny, 4); % 负通量
    g_1 = zeros(Nx, Ny, 4); % 正通量
    g_2 = zeros(Nx, Ny, 4); % 负通量
    h_1 = zeros(Nx, Ny, 4); % 插值后的f_1
    h_2 = zeros(Nx, Ny, 4);
    f_1x = zeros(Nx, Ny, 4); % f1x为正通量导数 f2x为负通量导数
    f_2x = zeros(Nx, Ny, 4);
    g_1y = zeros(Nx, Ny, 4); % g1x为正通量导数 f2x为负通量导数
    g_2y = zeros(Nx, Ny, 4);

    D_1 = zeros(4, 4); % 正特征值
    D_2 = zeros(4, 4); % 负特征值

    big = zeros(1, 2); % 记录通量导数最大值
    Ry = ones(4, 4);

    dt = CFL * dx / max(max(max(abs(W(:, :, 2)), abs(W(:, :, 3))) + sqrt(1.4 * W(:, :, 4) ./ W(:, :, 1))));

    ep = 1e-6;
    gam = zeros(1, 3);
    gam(1) = 1/16;
    gam(2) = 5/8;
    gam(3) = 5/16;
    beta = zeros(4, 4);
    omg = zeros(3, 4);
    omgm = zeros(4, 1);
    q = zeros(3, 4);

    for j = 1:Ny - 1

        for i = 1:Nx - 1
            p = W(i, j, 4); % 压强
            a = sqrt(1.4 * p / U(i, j, 1)); % 声速
            u = W(i, j, 2);
            v = W(i, j, 3);
            h = (U(i, j, 4) + p) / W(i, j, 1); % 单位质量流体热焓 《计算》P21

            % x 方向
            % 右特征向量矩阵
            Rx = ones(4, 4);
            Rx(2, 1) = u;
            Rx(3, 1) = 0;
            Rx(4, 1) = 0.5 * (u^2 - v^2);
            Rx(1, 2) = 0;
            Rx(2, 2) = 0;
            Rx(4, 2) = v;
            Rx(2, 3) = u - a;
            Rx(3, 3) = v;
            Rx(4, 3) = h - u * a;
            Rx(2, 4) = u + a;
            Rx(3, 4) = v;
            Rx(4, 4) = h + u * a;

            % 左特征向量矩阵
            Lx = ones(4, 4);
            b2 = 0.4 / a^2;
            b1 = b2 * (u^2 + v^2) / 2;
            Lx(1, 1) = 1 - b1;
            Lx(2, 1) = -b1 * v;
            Lx(3, 1) = 0.5 * (b1 + u / a);
            Lx(4, 1) = 0.5 * (b1 - u / a);
            Lx(1, 2) = b2 * u;
            Lx(2, 2) = b2 * u * v;
            Lx(3, 2) = -0.5 * (b2 * u + 1 / a);
            Lx(4, 2) = -0.5 * (b2 * u - 1 / a);
            Lx(1, 3) = b2 * v;
            Lx(2, 3) = 1 + b2 * v^2;
            Lx(3, 3) = -0.5 * b2 * v;
            Lx(4, 3) = -0.5 * b2 * v;
            Lx(1, 4) = -b2;
            Lx(2, 4) = -b2 * v;
            Lx(3, 4) = 0.5 * b2;
            Lx(4, 4) = 0.5 * b2;

            % Lx * Rx

            % 对角矩阵
            Gamma_x = diag([u, u, u - a, u + a]);
            tempU = reshape(U(i, j, :), 4, 1);
            f_1(i, j, :) = Rx * ((Gamma_x + abs(Gamma_x)) / 2) * Lx * tempU;
            f_2(i, j, :) = Rx * ((Gamma_x - abs(Gamma_x)) / 2) * Lx * tempU;

            % y 方向
            % 右特征向量矩阵
            Ry = ones(4, 4);
            Ry(1, 1) = 0;
            Ry(3, 1) = 0;
            Ry(4, 1) = u;
            Ry(2, 2) = 0;
            Ry(3, 2) = v;
            Ry(4, 2) = 0.5 * (v^2 - u^2);
            Ry(2, 3) = u;
            Ry(3, 3) = v - a;
            Ry(4, 3) = h - v * a;
            Ry(2, 4) = u;
            Ry(3, 4) = v + a;
            Ry(4, 4) = h + v * a;

            % 左特征向量矩阵
            Ly = ones(4, 4);
            Ly(1, 1) = -b1 * u;
            Ly(2, 1) = 1 - b1;
            Ly(3, 1) = 0.5 * (b1 + v / a);
            Ly(4, 1) = 0.5 * (b1 - v / a);
            Ly(1, 2) = 1 + b2 * u^2;
            Ly(2, 2) = b2 * u;
            Ly(3, 2) = -0.5 * b2 * u;
            Ly(4, 2) = -0.5 * b2 * u;
            Ly(1, 3) = b2 * v * u;
            Ly(2, 3) = b2 * v;
            Ly(3, 3) = -0.5 * (b2 * v + 1 / a);
            Ly(4, 3) = -0.5 * (b2 * v - 1 / a);
            Ly(1, 4) = -b2 * u;
            Ly(2, 4) = -b2;
            Ly(3, 4) = 0.5 * b2;
            Ly(4, 4) = 0.5 * b2;

            % Ly * Ry
            % 对角矩阵
            Gamma_y = diag([v, v, v - a, v + a]);
            g_1(i, j, :) = Ry * ((Gamma_y + abs(Gamma_y)) / 2) * Ly * tempU;
            g_2(i, j, :) = Ry * ((Gamma_y - abs(Gamma_y)) / 2) * Ly * tempU;

        end

    end

    % x 方向
    for i = 3:Nx - 2

        for j = 3:Ny - 2

            for k = 1:4
                beta(1, k) = (f_1(i - 2, j, k) * f_1(i - 2, j, k) * 4 - f_1(i - 2, j, k) * f_1(i - 1, j, k) * 19 + f_1(i - 1, j, k) * f_1(i - 1, j, k) * 25 + f_1(i - 2, j, k) * f_1(i, j, k) * 11 - f_1(i - 1, j, k) * f_1(i, j, k) * 31 + f_1(i, j, k) * f_1(i, j, k) * 10) / 3;

                beta(2, k) = (f_1(i - 1, j, k) * f_1(i - 1, j, k) * 4 - f_1(i - 1, j, k) * f_1(i, j, k) * 13 + f_1(i, j, k) * f_1(i, j, k) * 13 + f_1(i - 1, j, k) * f_1(i + 1, j, k) * 5 - f_1(i + 1, j, k) * f_1(i, j, k) * 13 + f_1(i + 1, j, k) * f_1(i + 1, j, k) * 4) / 3;

                beta(3, k) = (f_1(i + 2, j, k) * f_1(i + 2, j, k) * 4 - f_1(i + 2, j, k) * f_1(i + 1, j, k) * 19 + f_1(i + 1, j, k) * f_1(i + 1, j, k) * 25 + f_1(i + 2, j, k) * f_1(i, j, k) * 11 - f_1(i + 1, j, k) * f_1(i, j, k) * 31 + f_1(i, j, k) * f_1(i, j, k) * 10) / 3;
            end

            for k = 1:4

                for l = 1:3
                    omg(l, k) = gam(l) / ((ep + beta(l, k))^2);
                end

                omgm(k) = omg(1, k) + omg(2, k) + omg(3, k);
            end

            % 三种长度为3的模版的差值结果
            for k = 1:4
                q(1, k) = 0.375 * f_1(i - 2, j, k) - 1.25 * f_1(i - 1, j, k) + 1.875 * f_1(i, j, k);
                q(2, k) = -0.125 * f_1(i - 1, j, k) + 0.75 * f_1(i, j, k) + 0.375 * f_1(i + 1, j, k);
                q(3, k) = 0.375 * f_1(i, j, k) + 0.75 * f_1(i + 1, j, k) - 0.125 * f_1(i + 2, j, k);
            end

            % 模版组合
            for k = 1:4
                h_1(i, j, k) = (q(1, k) * omg(1, k) + q(2, k) * omg(2, k) + q(3, k) * omg(3, k)) / omgm(k);
            end

        end

    end

    for i = 3:Nx - 2

        for j = 3:Ny - 2

            for k = 1:4
                beta(1, k) = (f_2(i + 2, j, k)^2 * 4 - f_2(i + 2, j, k) * f_2(i + 1, j, k) * 19 + f_2(i + 1, j, k) * f_2(i + 1, j, k) * 25 + f_2(i + 2, j, k) * f_2(i, j, k) * 11 - f_2(i + 1, j, k) * f_2(i, j, k) * 31 + f_2(i, j, k) * f_2(i, j, k) * 10) / 3;

                beta(2, k) = (f_2(i + 1, j, k) * f_2(i + 1, j, k) * 4 - f_2(i + 1, j, k) * f_2(i, j, k) * 13 + f_2(i, j, k) * f_2(i, j, k) * 13 + f_2(i + 1, j, k) * f_2(i - 1, j, k) * 5 - f_2(i - 1, j, k) * f_2(i, j, k) * 13 + f_2(i - 1, j, k) * f_2(i - 1, j, k) * 4) / 3;

                beta(3, k) = (f_2(i - 2, j, k) * f_2(i - 2, j, k) * 4 - f_2(i - 2, j, k) * f_2(i - 1, j, k) * 19 + f_2(i - 1, j, k) * f_2(i - 1, j, k) * 25 + f_2(i - 2, j, k) * f_2(i, j, k) * 11 - f_2(i - 1, j, k) * f_2(i, j, k) * 31 + f_2(i, j, k) * f_2(i, j, k) * 10) / 3;
            end

            for k = 1:4

                for l = 1:3
                    omg(l, k) = gam(l) / ((ep + beta(l, k))^2);
                end

                omgm(k) = omg(1, k) + omg(2, k) + omg(3, k);
            end

            % 三种长度为3的模版的差值结果
            for k = 1:4
                q(1, k) = 0.375 * f_2(i + 2, j, k) - 1.25 * f_2(i + 1, j, k) + 1.875 * f_2(i, j, k);
                q(2, k) = -0.125 * f_2(i + 1, j, k) + 0.75 * f_2(i, j, k) + 0.375 * f_2(i - 1, j, k);
                q(3, k) = 0.375 * f_2(i, j, k) + 0.75 * f_2(i - 1, j, k) - 0.125 * f_2(i - 2, j, k);
            end

            % 模版组合
            for k = 1:4
                h_2(i - 1, j, k) = (q(1, k) * omg(1, k) + q(2, k) * omg(2, k) + q(3, k) * omg(3, k)) / omgm(k);
            end

        end

    end

    for i = 4:Nx - 2

        for j = 4:Ny - 2
            f_1x(i, j, :) = (h_1(i, j, :) - h_1(i - 1, j, :)) / dx;
        end

    end

    for i = 3:Nx - 3

        for j = 3:Ny - 3
            f_2x(i, j, :) = (h_2(i, j, :) - h_2(i - 1, j, :)) / dx;
        end

    end

    % y 方向

    for i = 3:Nx - 2

        for j = 3:Ny - 2

            for k = 1:4
                beta(1, k) = (g_1(i, j - 2, k) * g_1(i, j - 2, k) * 4 - g_1(i, j - 2, k) * g_1(i, j - 1, k) * 19 + g_1(i, j - 1, k) * g_1(i, j - 1, k) * 25 + g_1(i, j - 2, k) * g_1(i, j, k) * 11 - g_1(i, j - 1, k) * g_1(i, j, k) * 31 + g_1(i, j, k) * g_1(i, j, k) * 10) / 3;

                beta(2, k) = (g_1(i, j - 1, k) * g_1(i, j - 1, k) * 4 - g_1(i, j - 1, k) * g_1(i, j, k) * 13 + g_1(i, j, k) * g_1(i, j, k) * 13 + g_1(i, j - 1, k) * g_1(i, j + 1, k) * 5 - g_1(i, j + 1, k) * g_1(i, j, k) * 13 + g_1(i, j + 1, k) * g_1(i, j + 1, k) * 4) / 3;

                beta(3, k) = (g_1(i, j + 2, k) * g_1(i, j + 2, k) * 4 - g_1(i, j + 2, k) * g_1(i, j + 1, k) * 19 + g_1(i, j + 1, k) * g_1(i, j + 1, k) * 25 + g_1(i, j + 2, k) * g_1(i, j, k) * 11 - g_1(i, j + 1, k) * g_1(i, j, k) * 31 + g_1(i, j, k) * g_1(i, j, k) * 10) / 3;
            end

            for k = 1:4

                for l = 1:3
                    omg(l, k) = gam(l) / ((ep + beta(l, k))^2);
                end

                omgm(k) = omg(1, k) + omg(2, k) + omg(3, k);
            end

            % 三种长度为3的模版的差值结果
            for k = 1:4
                q(1, k) = 0.375 * g_1(i, j - 2, k) - 1.25 * g_1(i, j - 1, k) + 1.875 * g_1(i, j, k);
                q(2, k) = -0.125 * g_1(i, j - 1, k) + 0.75 * g_1(i, j, k) + 0.375 * g_1(i, j + 1, k);
                q(3, k) = 0.375 * g_1(i, j, k) + 0.75 * g_1(i, j + 1, k) - 0.125 * g_1(i, j + 2, k);
            end

            % 模版组合
            for k = 1:4
                h_1(i, j, k) = (q(1, k) * omg(1, k) + q(2, k) * omg(2, k) + q(3, k) * omg(3, k)) / omgm(k);
            end

        end

    end

    for i = 3:Nx - 2

        for j = 3:Ny - 2

            for k = 1:4
                beta(1, k) = (g_2(i, j + 2, k)^2 * 4 - g_2(i, j + 2, k) * g_2(i, j + 1, k) * 19 + g_2(i, j + 1, k) * g_2(i, j + 1, k) * 25 + g_2(i, j + 2, k) * g_2(i, j, k) * 11 - g_2(i, j + 1, k) * g_2(i, j, k) * 31 + g_2(i, j, k) * g_2(i, j, k) * 10) / 3;

                beta(2, k) = (g_2(i, j + 1, k) * g_2(i, j + 1, k) * 4 - g_2(i, j + 1, k) * g_2(i, j, k) * 13 + g_2(i, j, k) * g_2(i, j, k) * 13 + g_2(i, j + 1, k) * g_2(i, j - 1, k) * 5 - g_2(i, j - 1, k) * g_2(i, j, k) * 13 + g_2(i, j - 1, k) * g_2(i, j - 1, k) * 4) / 3;

                beta(3, k) = (g_2(i, j - 2, k) * g_2(i, j - 2, k) * 4 - g_2(i, j - 2, k) * g_2(i, j - 1, k) * 19 + g_2(i, j - 1, k) * g_2(i, j - 1, k) * 25 + g_2(i, j - 2, k) * g_2(i, j, k) * 11 - g_2(i, j - 1, k) * g_2(i, j, k) * 31 + g_2(i, j, k) * g_2(i, j, k) * 10) / 3;
            end

            for k = 1:4

                for l = 1:3
                    omg(l, k) = gam(l) / ((ep + beta(l, k))^2);
                end

                omgm(k) = omg(1, k) + omg(2, k) + omg(3, k);
            end

            % 三种长度为3的模版的差值结果
            for k = 1:4
                q(1, k) = 0.375 * g_2(i, j + 2, k) - 1.25 * g_2(i, j + 1, k) + 1.875 * g_2(i, j, k);
                q(2, k) = -0.125 * g_2(i, j + 1, k) + 0.75 * g_2(i, j, k) + 0.375 * g_2(i, j - 1, k);
                q(3, k) = 0.375 * g_2(i, j, k) + 0.75 * g_2(i, j - 1, k) - 0.125 * g_2(i, j - 2, k);
            end

            % 模版组合
            for k = 1:4
                h_2(i, j - 1, k) = (q(1, k) * omg(1, k) + q(2, k) * omg(2, k) + q(3, k) * omg(3, k)) / omgm(k);
            end

        end

    end

    for i = 4:Nx - 2

        for j = 4:Ny - 2
            g_1y(i, j, :) = (h_1(i, j, :) - h_1(i, j - 1, :)) / dx;
        end

    end

    for i = 3:Nx - 3

        for j = 3:Ny - 3
            g_2y(i, j, :) = (h_2(i, j, :) - h_2(i, j - 1, :)) / dx;
        end

    end

    for i = 4:Nx - 3

        for j = 4:Ny - 3
            flux(i, j, :) = f_1x(i, j, :) + f_2x(i, j, :) + g_1y(i, j, :) + g_2y(i, j, :);
        end

    end

end
