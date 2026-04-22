function solve_report_example3_viz()
% =========================================================================
% ВІЗУАЛІЗАЦІЯ ЯК У ЗВІТІ (Покращена)
% =========================================================================
    clc; close all;
    fprintf('Розрахунок та побудова графіків (High Quality)...\n');

    % --- 1. ПАРАМЕТРИЗАЦІЯ ---
    geom.x1 = @(t) cos(t); geom.y1 = @(t) sin(t);
    geom.dx1 = @(t) -sin(t); geom.dy1 = @(t) cos(t);
    geom.ddx1 = @(t) -cos(t); geom.ddy1 = @(t) -sin(t);

    geom.x2 = @(t) 3*cos(t); geom.y2 = @(t) 3*sin(t);
    geom.dx2 = @(t) -3*sin(t); geom.dy2 = @(t) 3*cos(t);
    geom.ddx2 = @(t) -3*cos(t); geom.ddy2 = @(t) -3*sin(t);

    y_star = [5, 5];
    u_exact = @(x) (1/(2*pi)) * log(1 / norm(x - y_star));
    du_dn_exact = @(x, nu) -1/(2*pi) * dot(x - y_star, nu) / sum((x - y_star).^2);

    % --- 2. РОЗРАХУНОК (N=64) ---
    N = 64; h = 2*pi/N; t = (0:N-1)'*h; w = 1/N;

    [P1, dP1, ddP1, nu1, jac1] = get_geom(t, geom.x1, geom.y1, geom.dx1, geom.dy1, geom.ddx1, geom.ddy1);
    nu1 = -nu1;
    [P2, dP2, ddP2, nu2, jac2] = get_geom(t, geom.x2, geom.y2, geom.dx2, geom.dy2, geom.ddx2, geom.ddy2);

    Dim = 2*N; A = zeros(Dim); RHS = zeros(Dim,1);

    for i=1:N
        for j=1:N
            if i==j, val=dot(ddP1(i,:),nu1(i,:))/(2*jac1(i)); else, r=P1(i,:)-P1(j,:); val=dot(r,nu1(j,:))/sum(r.^2)*jac1(j); end
            A(i,j) = val*w; if i==j, A(i,j)-=0.5; end
            r=P1(i,:)-P2(j,:); A(i,j+N) = log(1/norm(r))*jac2(j)*w;
            r=P2(i,:)-P1(j,:); d2=sum(r.^2); t1=dot(nu2(i,:),nu1(j,:)); t2=2*dot(r,nu2(i,:))*dot(r,nu1(j,:))/d2;
            A(i+N,j) = -(t1-t2)/d2*jac1(j)*w;
            if i==j, val=dot(ddP2(i,:),nu2(i,:))/(2*jac2(i)); else, r=P2(j,:)-P2(i,:); val=dot(r,nu2(i,:))/sum(r.^2)*jac2(j); end
            A(i+N,j+N) = val*w; if i==j, A(i+N,j+N)+=0.5; end
        end
        RHS(i) = u_exact(P1(i,:)); RHS(i+N) = du_dn_exact(P2(i,:), nu2(i,:));
    end
    Psi = A \ RHS; psi1 = Psi(1:N); psi2 = Psi(N+1:end);

    % --- 3. ВІЗУАЛІЗАЦІЯ ---

    % Налаштування сітки (Густіша для плавності)
    n_grid = 60; % Більше точок
    [X, Y] = meshgrid(linspace(-3.5, 3.5, n_grid));
    Q = nan(size(X));

    fprintf('Рендеринг...\n');
    for i = 1:numel(X)
        pt = [X(i), Y(i)]; r = norm(pt);
        % Умова області: 1 < r < 3 (з малим запасом, щоб краї були рівні)
        if r >= 1.0 && r <= 3.0
            Q(i) = get_u(pt, psi1, psi2, P1, nu1, jac1, P2, jac2);
        end
    end

    % ФІГУРА 1: MESH (Як Рис. 4)
    figure(1); clf;
    mesh(X, Y, Q);
    view(-45, 40); % Підбираємо кут
    axis([-4 4 -4 4]);
    colormap jet;
    title('Рис. 4 Графік наближеного розв''язку');
    xlabel('x'); ylabel('y');

    % ФІГУРА 2: CONTOUR (Як Рис. 5)
    figure(2); clf; hold on; axis equal;

    % Малюємо границі (щоб було видно область)
    tt = linspace(0, 2*pi, 200);
    plot(cos(tt), sin(tt), 'k-');
    plot(3*cos(tt), 3*sin(tt), 'k-');

    % Лінії рівня (більше ліній для краси)
    [C, h] = contour(X, Y, Q, 25);
    set(h, 'LineWidth', 1.2); % Трохи товстіші лінії
    colorbar; colormap jet;
    axis([-3.5 3.5 -3.5 3.5]);
    title('Рис. 5 Лінії рівня');

    fprintf('Готово.\n');
end

function [P, dP, ddP, nu, jac] = get_geom(t, fx, fy, fdx, fdy, fddx, fddy)
    N = length(t); P=[fx(t), fy(t)]; dP=[fdx(t), fdy(t)]; ddP=[fddx(t), fddy(t)];
    jac = sqrt(sum(dP.^2, 2)); nu = [dP(:,2), -dP(:,1)] ./ jac;
end
function u = get_u(pt, psi1, psi2, P1, nu1, jac1, P2, jac2)
    N = length(psi1);
    R1 = pt - P1; term1 = sum(R1 .* nu1, 2) ./ sum(R1.^2, 2);
    int1 = sum(psi1 .* term1 .* jac1);
    R2 = pt - P2; term2 = log(1 ./ sqrt(sum(R2.^2, 2)));
    int2 = sum(psi2 .* term2 .* jac2);
    u = (1/N) * (int1 + int2);
end
