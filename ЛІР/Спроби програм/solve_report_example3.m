function solve_report_example3()
% =========================================================================
% РЕАЛІЗАЦІЯ ПРИКЛАДУ №3 (ВИПРАВЛЕНО N -> 2N)
% =========================================================================
    clc; close all;
    fprintf('Запуск Прикладу №3 (R1=1, R2=3)...\n\n');

    % --- 1. ПАРАМЕТРИЗАЦІЯ ---
    geom.x1   = @(t) [cos(t), sin(t)];
    geom.dx1  = @(t) [-sin(t), cos(t)];
    geom.ddx1 = @(t) [-cos(t), -sin(t)];

    geom.x2   = @(t) [3.*cos(t), 3.*sin(t)];
    geom.dx2  = @(t) [-3.*sin(t), 3.*cos(t)];
    geom.ddx2 = @(t) [-3.*cos(t), -3.*sin(t)];

    % --- 2. ТОЧНИЙ РОЗВ'ЯЗОК ---
    y_star = [5, 5];
    u_exact = @(x) (1/(2*pi)) * log(1 / norm(x - y_star));
    du_dn_exact = @(x, nu) -1/(2*pi) * dot(x - y_star, nu) / sum((x - y_star).^2);

    % --- 3. ТЕСТОВІ ТОЧКИ ---
    test_points = [2.0, 0.0; 1.0, 1.0; 0.0, -2.0];

    % --- 4. ЦИКЛ ПО N ---
    % N - це параметр. Кількість точок буде 2*N.
    N_values = [4, 8, 16, 32, 64];

    fprintf('------------------------------------------------------------\n');
    fprintf('  N | 2N (Pts)|  Err(2,0)   |  Err(1,1)   |  Err(0,-2) \n');
    fprintf('------------------------------------------------------------\n');

    res = struct();

    for N = N_values
        % ВУЗЛІВ ТЕПЕР 2*N
        num_nodes = 2 * N;

        h = 2*pi / num_nodes;
        t = (0 : num_nodes-1)' * h;

        % ВАГА КВАДРАТУРИ: 1 / (2N)
        w = 1 / num_nodes;

        [P1, dP1, ddP1, nu1, jac1] = get_geom(t, geom.x1, geom.dx1, geom.ddx1);
        nu1 = -nu1; % Нормаль в дірку
        [P2, dP2, ddP2, nu2, jac2] = get_geom(t, geom.x2, geom.dx2, geom.ddx2);

        % РОЗМІР СИСТЕМИ: 2 * (2N) = 4N
        Dim = 2 * num_nodes;
        A = zeros(Dim, Dim);
        RHS = zeros(Dim, 1);

        for i = 1:num_nodes
            for j = 1:num_nodes
                % L11 (Г1->Г1)
                if i == j
                    val = dot(ddP1(i,:), nu1(i,:)) / (2 * jac1(i));
                else
                    r = P1(i,:) - P1(j,:);
                    val = dot(r, nu1(j,:)) / sum(r.^2) * jac1(j);
                end
                A(i, j) = val * w;
                if i == j, A(i, j) = A(i, j) - 0.5; end

                % L12 (Г2->Г1)
                r = P1(i,:) - P2(j,:);
                val = log(1 / norm(r)) * jac2(j);
                % Зсув стовпчика на num_nodes
                A(i, j+num_nodes) = val * w;

                % L21 (Г1->Г2)
                r_vec = P2(i,:) - P1(j,:); d2 = sum(r_vec.^2);
                t1 = dot(nu2(i,:), nu1(j,:));
                t2 = 2 * dot(r_vec, nu2(i,:)) * dot(r_vec, nu1(j,:)) / d2;
                val = -(t1 - t2) / d2 * jac1(j);
                % Зсув рядка на num_nodes
                A(i+num_nodes, j) = val * w;

                % L22 (Г2->Г2)
                if i == j
                    val = dot(ddP2(i,:), nu2(i,:)) / (2 * jac2(i));
                else
                    r = P2(j,:) - P2(i,:);
                    val = dot(r, nu2(i,:)) / sum(r.^2) * jac2(j);
                end
                % Зсув і рядка, і стовпчика
                A(i+num_nodes, j+num_nodes) = val * w;
                if i == j, A(i+num_nodes, j+num_nodes) = A(i+num_nodes, j+num_nodes) + 0.5; end
            end

            RHS(i) = u_exact(P1(i,:));
            RHS(i+num_nodes) = du_dn_exact(P2(i,:), nu2(i,:));
        end

        if N < 4, Psi = pinv(A)*RHS; else, Psi = A\RHS; end

        psi1 = Psi(1:num_nodes);
        psi2 = Psi(num_nodes+1:end);

        if N == 64
            res.psi1 = psi1; res.psi2 = psi2; res.P1 = P1; res.P2 = P2; res.nu1 = nu1; res.jac1 = jac1; res.jac2 = jac2;
        end

        errs = zeros(1, 3);
        for k = 1:3
            u_val = get_u(test_points(k,:), psi1, psi2, P1, nu1, jac1, P2, jac2);
            errs(k) = abs(u_val - u_exact(test_points(k,:)));
        end

        fprintf(' %-3d | %-7d  | %.3e  | %.3e  | %.3e \n', N, num_nodes, errs);
    end

    % --- ГРАФІКИ (Без змін) ---
    plot_results(geom, test_points, y_star, res);
end

% Допоміжні функції (ті самі)
function [P, dP, ddP, nu, jac] = get_geom(t, fx, fdx, fddx)
    N = length(t); P=[fx(t)]; dP=[fdx(t)]; ddP=[fddx(t)];
    jac = sqrt(sum(dP.^2, 2)); nu = [dP(:,2), -dP(:,1)] ./ jac;
end

function u = get_u(pt, psi1, psi2, P1, nu1, jac1, P2, jac2)
    N = length(psi1); % Це вже 2*N
    w = 1 / N;        % Вага 1/(2N)

    R1 = pt - P1; term1 = sum(R1 .* nu1, 2) ./ sum(R1.^2, 2);
    int1 = sum(psi1 .* term1 .* jac1);
    R2 = pt - P2; term2 = log(1 ./ sqrt(sum(R2.^2, 2)));
    int2 = sum(psi2 .* term2 .* jac2);

    u = w * (int1 + int2);
end

function plot_results(geom, pts, y_star, res)
    % (Код малювання такий самий, як був, з surf та contour)
    figure(1); clf; hold on; axis equal; grid on;
    tt = linspace(0, 2*pi, 200);
    plot(geom.x1(tt), 'b'); plot(geom.x2(tt), 'r');
    plot(pts(:,1), pts(:,2), 'ko', 'MarkerFaceColor', 'g');

    figure(2); clf;
    [X, Y] = meshgrid(linspace(-3.5, 3.5, 80));
    U = nan(size(X));
    for k=1:numel(X)
        pt = [X(k), Y(k)];
        if norm(pt) > 1.05 && norm(pt) < 2.95
            U(k) = get_u(pt, res.psi1, res.psi2, res.P1, res.nu1, res.jac1, res.P2, res.jac2);
        end
    end
    surf(X, Y, real(U)); shading interp; view(2); axis equal; colorbar;
    hold on; plot(geom.x1(tt), 'w'); plot(geom.x2(tt), 'k');
end
