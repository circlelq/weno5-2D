%
%  main1.m
%  weno5-2D
%
%  Created by Yuan Leiqi (袁磊祺) on 2021/7/8.
%

clear;clc;close all;
set(0, 'defaultlinelinewidth', 3)
set(0, 'defaultaxeslinewidth', 2);
set(0, 'defaultaxesfontsize', 28);
set(0, 'defaulttextfontsize', 28);
set(0, 'DefaultLineMarkerSize', 2);
set(0, 'Defaultaxesfontname', 'Times New Roman');
% addpath('./circleData');
% set(gcf,'unit','centimeters','position',[20 20 20 20])
% figure(); box on;
% plot(x/D,W0_axi,'-');
% xlabel('$x/D$','interpreter','latex');
% ylabel('$U_0$','interpreter','latex');

% saveas(gcf,'Re_80_U0_axi','epsc')

global gamma
global dx dy
global Nx Ny
global CFL
gamma = 1.4;

Ny = 50; % y grid number
Nx = Ny * 3; % x grid number
CFL = 0.9; % CFL number
x = linspace(0, 3, Nx)'; % x grid
dx = (x(end) - x(1)) / (Nx - 1); % x grid spacing
y = linspace(0, 1, Ny)'; % y grid
dy = (y(end) - y(1)) / (Ny - 1); % y grid spacing

W0 = [1.4, 3, 0, 1];
W = zeros(Nx, Ny, 4);

for i = 1:Nx

    for j = 1:Ny
        W(i, j, :) = W0;
    end

end

dt = 1e-2;
dT = 0; % every 0.1 time to record
U = W2U(W);
steps = 0;
flag = 1;
current_time = 0;
t_max = 4;
% maxSteps     = 3e4;
h = (Ny * 0.2);
d = (Nx * 0.2);

aviobj = VideoWriter('weno5', 'MPEG-4');
open(aviobj);
fig = figure;
set(gcf, 'unit', 'centimeters', 'position', [20 20 90 30])

while (flag)
    steps = steps + 1;

    % timestep
    disp(current_time);

    % boundary conditon
    W = U2W(U);
    W = bc(W);
    U = W2U(W);

    if current_time + dt >= t_max % stoping criteria
        dt = t_max - current_time;
        disp('Finished!');
        flag = 0;
    end

    current_time = current_time + dt;

    % weno5 scheme
    [tempflux, dt] = weno5(U);
    tempU1 = U - dt * tempflux;
    W = U2W(tempU1);
    W = bc(W);
    tempU1 = W2U(W);
    [tempflux, ~] = weno5(tempU1);
    tempU2 = 0.75 * U + 0.25 * tempU1 - 0.25 * dt * tempflux;
    W = U2W(tempU2);
    W = bc(W);
    tempU2 = W2U(W);
    [tempflux, ~] = weno5(tempU2);
    U = 1/3 * U + 2/3 * tempU2 - 2/3 * dt * tempflux;

    % visualization

    W = U2W(U);

    for i = 1:d - 1

        for j = 3:Ny - 3
            W1(i, j, :) = W(i, j, :);
        end

    end

    for i = d:Nx - 3

        for j = h + 1:Ny - 3
            W1(i, j, :) = W(i, j, :);
        end

    end

    contourf(W1(:, :, 1))
    xlabel('$y$', 'interpreter', 'latex');
    ylabel('$x$', 'interpreter', 'latex');
    title(sprintf('weno5 scheme, time = %.3f', current_time), 'interpreter', 'latex')

    view(90, -90)
    drawnow
    currFrame = getframe(fig);
    % if abs(current_time - dT) < 1e-2
    writeVideo(aviobj, currFrame);
    % save(strcat(num2str(current_time),'W.mat'),'W')
    % dT = dT + 0.1;
    % end
end

close(aviobj); %关闭

% save W.mat W
