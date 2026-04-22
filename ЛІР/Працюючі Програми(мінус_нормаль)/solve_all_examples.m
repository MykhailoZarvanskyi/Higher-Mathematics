function solve_all_examples()
    clc; close all;
    fprintf('Оберіть приклад для розрахунку:\n');
    fprintf('1 - Коло (R=1.1) і Коло (R=2.5)\n');
    fprintf('2 - Еліпс (1.2, 0.75) і Коло (R=3.5)\n');
    fprintf('3 - Повітряний змій всередині Еліпса (3.5, 2.5)\n');
    choice = input('Ваш вибір (1-3): ');
    if isempty(choice), choice = 1; end

    if choice == 1
        fprintf('\nЗапуск Прикладу 1: Два кола (R=1.1, R=2.5) \n');
        R1 = 1.1;
        geom.x1 = @(t) R1*cos(t);   geom.y1 = @(t) R1*sin(t);
        geom.dx1 = @(t) -R1*sin(t); geom.dy1 = @(t) R1*cos(t);
        geom.ddx1 = @(t) -R1*cos(t); geom.ddy1 = @(t) -R1*sin(t);

        R2 = 2.5;
        geom.x2 = @(t) R2*cos(t);  geom.y2 = @(t) R2*sin(t);
        geom.dx2 = @(t) -R2*sin(t); geom.dy2 = @(t) R2*cos(t);
        geom.ddx2 = @(t) -R2*cos(t); geom.ddy2 = @(t) -R2*sin(t);

    elseif choice == 2
        fprintf('\nЗапуск Прикладу 2: Еліпс (1.2, 0.75) у Колі (R=3.5) \n');
        a=1.2; b=0.75;
        geom.x1 = @(t) a*cos(t);   geom.y1 = @(t) b*sin(t);
        geom.dx1 = @(t) -a*sin(t); geom.dy1 = @(t) b*cos(t);
        geom.ddx1 = @(t) -a*cos(t); geom.ddy1 = @(t) -b*sin(t);

        R_out = 3.5;
        geom.x2 = @(t) R_out*cos(t);  geom.y2 = @(t) R_out*sin(t);
        geom.dx2 = @(t) -R_out*sin(t); geom.dy2 = @(t) R_out*cos(t);
        geom.ddx2 = @(t) -R_out*cos(t); geom.ddy2 = @(t) -R_out*sin(t);

    elseif choice == 3
        fprintf('\nЗапуск Прикладу 3: Змій у Еліпсі \n');
        geom.x1 = @(t) cos(t) + 0.65*cos(2*t) - 0.65;
        geom.y1 = @(t) 1.5*sin(t);
        geom.dx1 = @(t) -sin(t) - 1.3*sin(2*t);
        geom.dy1 = @(t) 1.5*cos(t);
        geom.ddx1 = @(t) -cos(t) - 2.6*cos(2*t);
        geom.ddy1 = @(t) -1.5*sin(t);

        a_out = 3.5; b_out = 2.5;
        geom.x2 = @(t) a_out*cos(t);  geom.y2 = @(t) b_out*sin(t);
        geom.dx2 = @(t) -a_out*sin(t); geom.dy2 = @(t) b_out*cos(t);
        geom.ddx2 = @(t) -a_out*cos(t); geom.ddy2 = @(t) -b_out*sin(t);
    end

    y_star = [5, 5];
    u_exact = @(x) (1/(2*pi)) * log(1 / norm(x - y_star));
    du_dn_exact = @(x, nu) -1/(2*pi) * dot(x - y_star, nu) / sum((x - y_star).^2);

    test_points = [1.8, 0.1;  0.0, 2.0;  -1.5, -0.1];

    N_values = [4, 8, 16, 32, 64];

    fprintf('\nТАБЛИЦЯ ПОХИБОК:\n');
    fprintf('------------------------------------------------------------\n');
    fprintf(' N  | 2N (Вузлів) |  Err(Pt1)   |  Err(Pt2)   |  Err(Pt3) \n');
    fprintf('------------------------------------------------------------\n');

    res = struct();

    for N = N_values
        M = 2 * N;
        h = 2*pi / M;
        t = (0 : M-1)' * h;
        w = 1 / M;

        [P1, dP1, ddP1, nu1, jac1] = get_geom(t, geom.x1, geom.y1, geom.dx1, geom.dy1, geom.ddx1, geom.ddy1);
        nu1 = -nu1;
        [P2, dP2, ddP2, nu2, jac2] = get_geom(t, geom.x2, geom.y2, geom.dx2, geom.dy2, geom.ddx2, geom.ddy2);

        Dim = 2 * M;
        A = zeros(Dim, Dim);
        RHS = zeros(Dim, 1);

        for i = 1:M
            for j = 1:M
                if i==j, val = dot(ddP1(i,:), nu1(i,:)) / (2 * jac1(i));
                else, r = P1(i,:) - P1(j,:); val = dot(r, nu1(j,:)) / sum(r.^2) * jac1(j); end
                A(i, j) = val * w;
                if i == j, A(i, j) -= 0.5; end

                r = P1(i,:) - P2(j,:);
                A(i, j+M) = log(1 / norm(r)) * jac2(j) * w;

                r_vec = P2(i,:) - P1(j,:); d2 = sum(r_vec.^2);
                t1 = dot(nu2(i,:), nu1(j,:));
                t2 = 2 * dot(r_vec, nu2(i,:)) * dot(r_vec, nu1(j,:)) / d2;
                val = -(t1 - t2) / d2 * jac1(j);
                A(i+M, j) = val * w;

                if i==j, val = dot(ddP2(i,:), nu2(i,:)) / (2 * jac2(i));
                else, r = P2(j,:) - P2(i,:); val = dot(r, nu2(i,:)) / sum(r.^2) * jac2(j); end
                A(i+M, j+M) = val * w;
                if i == j, A(i+M, j+M) += 0.5; end
            end

            RHS(i) = u_exact(P1(i,:));
            RHS(i+M) = du_dn_exact(P2(i,:), nu2(i,:));
        end

        if M < 8, Psi = pinv(A)*RHS; else, Psi = A\RHS; end

        psi1 = Psi(1:M);
        psi2 = Psi(M+1:end);

        if N == 64
            res.psi1 = psi1; res.psi2 = psi2;
            res.P1 = P1; res.P2 = P2; res.nu1 = nu1; res.jac1 = jac1; res.jac2 = jac2;
        end

        errs = zeros(1, 3);
        for k = 1:3
            u_val = get_u(test_points(k,:), psi1, psi2, P1, nu1, jac1, P2, jac2);
            errs(k) = abs(u_val - u_exact(test_points(k,:)));
        end
        fprintf(' %-2d | %-10d  | %.3e  | %.3e  | %.3e \n', N, M, errs);
    end
    fprintf('------------------------------------------------------------\n');

    figure(1); clf; hold on; axis equal; grid on; box on; set(gca, 'FontSize', 12);
    tt = linspace(0, 2*pi, 200);

    if choice == 1
        p1_dr = [1.1*cos(tt); 1.1*sin(tt)]'; p2_dr = [2.5*cos(tt); 2.5*sin(tt)]';
    elseif choice == 2
        p1_dr = [1.2*cos(tt); 0.75*sin(tt)]'; p2_dr = [3.5*cos(tt); 3.5*sin(tt)]';
    elseif choice == 3
        p1_dr = [cos(tt)+0.65*cos(2*tt)-0.65; 1.5*sin(tt)]'; p2_dr = [3.5*cos(tt); 2.5*sin(tt)]';
    end

    plot(p1_dr(:,1), p1_dr(:,2), 'b-', 'LineWidth', 2);
    plot(p2_dr(:,1), p2_dr(:,2), 'r-', 'LineWidth', 2);
    plot(test_points(:,1), test_points(:,2), 'ko', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
    title('Геометрія');

    bound = 4;
    n_grid = 10;
    [X, Y] = meshgrid(linspace(-bound, bound, n_grid));
    Q = nan(size(X));

    fprintf('Побудова графіків...\n');
    for k=1:numel(X)
        pt = [X(k), Y(k)];

        in_outer = inpolygon(pt(1), pt(2), p2_dr(:,1), p2_dr(:,2));
        in_inner = inpolygon(pt(1), pt(2), p1_dr(:,1), p1_dr(:,2));

        d1 = min(sqrt(sum((p1_dr - pt).^2, 2)));
        d2 = min(sqrt(sum((p2_dr - pt).^2, 2)));

        if in_outer && ~in_inner && d1 > 0.1 && d2 > 0.1
            Q(k) = get_u(pt, res.psi1, res.psi2, res.P1, res.nu1, res.jac1, res.P2, res.jac2);
        end
    end

    figure(2); clf;
    mesh(X, Y, real(Q));
    view(-30, 60); axis([-bound bound -bound bound]); colormap jet;
    title('Графік наближеного розв''язку (Mesh)');
    xlabel('x'); ylabel('y'); zlabel('u(x)');

    figure(3); clf; hold on; axis equal; box on;
    contour(X, Y, real(Q), 25, 'LineWidth', 1.5);
    plot(p1_dr(:,1), p1_dr(:,2), 'k-', 'LineWidth', 2);
    plot(p2_dr(:,1), p2_dr(:,2), 'k-', 'LineWidth', 2);
    colorbar; colormap jet;
    title('Лінії рівня');
end

function [P, dP, ddP, nu, jac] = get_geom(t, fx, fy, fdx, fdy, fddx, fddy)
    N = length(t);
    P = [fx(t), fy(t)];
    dP = [fdx(t), fdy(t)];
    ddP = [fddx(t), fddy(t)];
    jac = sqrt(sum(dP.^2, 2));
    nu = [dP(:,2), -dP(:,1)] ./ jac;
end

function u = get_u(pt, psi1, psi2, P1, nu1, jac1, P2, jac2)
    M = length(psi1);
    w = 1 / M;

    R1 = pt - P1; term1 = sum(R1 .* nu1, 2) ./ sum(R1.^2, 2);
    int1 = sum(psi1 .* term1 .* jac1);

    R2 = pt - P2; term2 = log(1 ./ sqrt(sum(R2.^2, 2)));
    int2 = sum(psi2 .* term2 .* jac2);

    u = w * (int1 + int2);
end
