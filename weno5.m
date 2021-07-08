%
%  weno5.m
%  weno5-2D
%
%  Created by Yuan Leiqi (袁磊祺) on 2021/7/8.
%

function [flux, dt] = weno5(U)
    % WENO5 scheme
    global Nx Ny

    flux = zeros(4, N); % flux vector

    W = U2W(U);
    f_1 = zeros(4, N); % 正通量
    f_2 = zeros(4, N); % 负通量
    h_1 = zeros(4, N); % 插值后的f_1
    h_2 = zeros(4, N);
    f_1x = zeros(4, N); %f1x为正通量导数 f2x为负通量导数
    f_2x = zeros(4, N);


    D_1 = zeros(4, 4); % 正特征值
    D_2 = zeros(4, 4); % 负特征值

    big = zeros(1, 2); % 记录通量导数最大值
	Ry = ones(4, 4);

    for i = 1:Nx-1
    	for j = 1:Ny-1
	        p = W(i, j, 4); % 压强
	        a = sqrt(1.4 * p / U(i, j, 1)); % 声速
	        H = (U(i, j, 4)+W(i, j, 4)) / U(i, j, 1); % 焓
	        u = W(i, j, 2);
	        v = W(i, j, 3);

	        % x 方向
	        % 右特征向量矩阵
	    	Rx = ones(4, 4);
	    	R_x(2, 1) = u;
	    	R_x(3, 1) = 0;
	    	R_x(4, 1) = 0.5*(u^2-v^2);
	    	R_x(1, 2) = 0;
	    	R_x(2, 2) = 0;
	    	R_x(4, 2) = v;
	    	R_x(2, 3) = u - a;
	    	R_x(3, 3) = v;
	    	R_x(4, 3) = H - u * a;
	    	R_x(2, 4) = u + a;
	    	R_x(3, 4) = v;
	    	R_x(2, 4) = H + u * a;

	        % 左特征向量矩阵
			Lx = ones(4, 4);
			b2 = 0.4 / a^2;
			b1 = b2 * (u^2 + v^2)/2;
	        Lx(1, 1) = 1 - b1;
	        Lx(2, 1) = - b1 * v;
	        Lx(3, 1) = 0.5 * (b1 + u/a);
	        Lx(4, 1) = 0.5 * (b1 - u/a);
	        Lx(1, 2) = b2 * u;
	        Lx(2, 2) = b2 * u * v;
	        Lx(3, 2) = -0.5 * (b2 * u + 1/a);
	        Lx(4, 2) = -0.5 * (b2 * u - 1/a);
	        Lx(1, 3) = b2 * v;
	        Lx(2, 3) = 1 + b2 * v^2;
	        Lx(3, 3) = -0.5 * b2 * v;
	        Lx(4, 3) = -0.5 * b2 * v;
	        Lx(1, 4) = -b2;
	        Lx(2, 4) = -b2 * v;
	        Lx(3, 4) = 0.5 * b2;
	        Lx(4, 4) = 0.5 * b2;

	        Lx * Rx;
	        % 对角矩阵
	        Gammax = diag([u, u, u-a, u+a]);
	        f_1(:, i) = Rx * ((Gamma+abs(Gamma)) / 2) * Lx * U(i, j, :)';
	        f_2(:, i) = Rx * ((Gamma-abs(Gamma)) / 2) * Lx * U(i, j, :)';

	        % Jacobi 矩阵
	        % A = [0, 1, 0;
	        %     -0.8 * u^2, 1.6 * u, 0.4;
	        %     -0.3 * u^3-a^2 * u / 0.4, 0.1 * u^2+a^2/0.4, 1.4 * u];
	    end
    end

    dt = CFL * dx / max(abs(W(:, 2))+sqrt(1.4 * W(:, 3) ./ W(:, 1)));

    ep = 1e-6;
    gam = zeros(1, 3);
    gam(1) = 1/16;
    gam(2) = 5/8;
    gam(3) = 5/16;
    beta = zeros(4, 4);
    omg = zeros(4, 4);
    omgm = zeros(3, 1);
    q = zeros(4, 4);

    for i = 3:N-2

        for j = 1:3
            beta(j, 1) = (f_1(j, i-2) * f_1(j, i-2) * 4-f_1(j, i-2) * f_1(j, i-1) * 19+f_1(j, i-1) * f_1(j, i-1) * 25+f_1(j, i-2) * f_1(j, i) * 11-f_1(j, i-1) .* f_1(j, i) * 31+f_1(j, i) * f_1(j, i) * 10) / 3;

            beta(j, 2) = (f_1(j, i-1) * f_1(j, i-1) * 4-f_1(j, i-1) * f_1(j, i) * 13+f_1(j, i) * f_1(j, i) * 13+f_1(j, i-1) * f_1(j, i+1) * 5-f_1(j, i+1) .* f_1(j, i) * 13+f_1(j, i+1) * f_1(j, i+1) * 4) / 3;

            beta(j, 3) = (f_1(j, i+2) * f_1(j, i+2) * 4-f_1(j, i+2) * f_1(j, i+1) * 19+f_1(j, i+1) * f_1(j, i+1) * 25+f_1(j, i+2) * f_1(j, i) * 11-f_1(j, i+1) * f_1(j, i) * 31+f_1(j, i) * f_1(j, i) * 10) / 3;
        end

        for j = 1:3

            for k = 1:3
                omg(j, k) = gam(k) / ((ep+beta(j, k))^2);
            end

            omgm(j) = omg(j, 1)+omg(j, 2)+omg(j, 3);
        end

        for j = 1:3
            q(j, 1) = 0.375 * f_1(j, i-2)-1.25 * f_1(j, i-1)+1.875 * f_1(j, i);
            q(j, 2) = -0.125 * f_1(j, i-1)+0.75 * f_1(j, i)+0.375 * f_1(j, i+1);
            q(j, 3) = 0.375 * f_1(j, i)+0.75 * f_1(j, i+1)-0.125 * f_1(j, i+2);
        end

        for j = 1:3
            h_1(j, i) = (q(j, 1) * omg(j, 1)+q(j, 2) * omg(j, 2)+q(j, 3) * omg(j, 3)) / omgm(j);
        end

    end

    for i = 3:N-2

        for j = 1:3
            beta(j, 1) = (f_2(j, i+2)^2 * 4-f_2(j, i+2) * f_2(j, i+1) * 19+f_2(j, i+1) * f_2(j, i+1) * 25+f_2(j, i+2) * f_2(j, i) * 11-f_2(j, i+1) * f_2(j, i) * 31+f_2(j, i) * f_2(j, i) * 10) / 3;

            beta(j, 2) = (f_2(j, i+1) * f_2(j, i+1) * 4-f_2(j, i+1) * f_2(j, i) * 13+f_2(j, i) * f_2(j, i) * 13+f_2(j, i+1) * f_2(j, i-1) * 5-f_2(j, i-1) * f_2(j, i) * 13+f_2(j, i-1) * f_2(j, i-1) * 4) / 3;

            beta(j, 3) = (f_2(j, i-2) * f_2(j, i-2) * 4-f_2(j, i-2) * f_2(j, i-1) * 19+f_2(j, i-1) * f_2(j, i-1) * 25+f_2(j, i-2) * f_2(j, i) * 11-f_2(j, i-1) * f_2(j, i) * 31+f_2(j, i) * f_2(j, i) * 10) / 3;
        end

        for j = 1:3

            for k = 1:3
                omg(j, k) = gam(k) / ((ep+beta(j, k))^2);
            end

            omgm(j) = omg(j, 1)+omg(j, 2)+omg(j, 3);
        end

        for j = 1:3
            q(j, 1) = 0.375 * f_2(j, i+2)-1.25 * f_2(j, i+1)+1.875 * f_2(j, i);
            q(j, 2) = -0.125 * f_2(j, i+1)+0.75 * f_2(j, i)+0.375 * f_2(j, i-1);
            q(j, 3) = 0.375 * f_2(j, i)+0.75 * f_2(j, i-1)-0.125 * f_2(j, i-2);
        end

        for j = 1:3
            h_2(j, i-1) = (q(j, 1) * omg(j, 1)+q(j, 2) * omg(j, 2)+q(j, 3) * omg(j, 3)) / omgm(j);
        end

    end

    for i = 4:N-2
        f_1x(:, i) = (h_1(:, i)-h_1(:, i-1)) / dx;
    end

    for i = 3:N-3
        f_2x(:, i) = (h_2(:, i)-h_2(:, i-1)) / dx;
    end

    for i = 4:N-3
        flux(:, i) = f_1x(:, i)+f_2x(:, i);
    end

end
