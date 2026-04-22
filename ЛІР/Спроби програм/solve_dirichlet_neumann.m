function solve_report_example3()
% =========================================================================
% ПРОГРАМА: Чисельне розв'язування мішаної задачі (Приклад №3)
% Метод: Нистрьом для системи інтегральних рівнянь 2-го роду
% =========================================================================
    clc; close all;
    fprintf('=== ЗАПУСК РОЗРАХУНКУ (ПРИКЛАД №3) ===\n');
    fprintf('Геометрія: Кільце (R1=1, R2=3)\n');
    fprintf('Точний розв''язок: від джерела y*=(5,5)\n\n');

    % --- 1. ВХІДНІ ДАНІ ---

    % Параметризація Г1 (Внутрішня, Діріхле): Коло R=1
    geom.x1   = @(t) cos(t);
    geom.y1   = @(t) sin(t);
    % Похідні (для точності задаємо аналітично)
    geom.dx1  = @(t) -sin(t); geom.dy1  = @(t) cos(t);
    geom.ddx1 = @(t) -cos(t); geom.ddy1 = @(t) -sin(t);

    % Параметризація Г2 (Зовнішня, Нейман): Коло R=3
    geom.x2   = @(t) 3*cos(t);
    geom.y2   = @(t) 3*sin(t);
    geom.dx2  = @(t) -3*sin(t); geom.dy2  = @(t) 3*cos(t);
    geom.ddx2 = @(t) -3*cos(t); geom.ddy2 = @(t) -3*sin(t);

    % Точний розв'язок та його похідна
    y_star = [5, 5];
    u_ex = @(x) (1/(2*pi)) * log(1 / norm(x - y_star));
    du_dn_ex = @(x, n) -1/(2*pi) * dot(x - y_star, n) / sum((x - y_star).^2);

    % Тестові точки
    pts = [2.0, 0.0; 1.0, 1.0; 0.0, -2.0];

    % --- 2. ДОСЛІДЖЕННЯ ЗБІЖНОСТІ ---
    N_vals = [4, 8, 16, 32, 64];

    fprintf('------------------------------------------------------------\n');
    fprintf(' N  |  Err(2,0)   |  Err(1,1)   |  Err(0,-2) \n');
    fprintf('------------------------------------------------------------\n');

    % Змінні для зберігання результатів останньої ітерації
    res = struct();

    for N = N_vals
        % 2.1. Дискретизація
        h = 2*pi / N;
        t = (0 : N-1)' * h;

        % 2.2. Геометрія у вузлах
        % is_inner=1 для Г1 (інвертує нормаль), is_inner=0 для Г2
        [P1, n1, dP1, ddP1, J1] = compute_geom(t, geom.x1, geom.y1, geom.dx1, geom.dy1, geom.ddx1, geom.ddy1, 1);
        [P2, n2, dP2, ddP2, J2] = compute_geom(t, geom.x2, geom.y2, geom.dx2, geom.dy2, geom.ddx2, geom.ddy2, 0);

        % 2.3. Формування СЛАР
        [A, RHS] = assemble_system(N, P1, n1, dP1, ddP1, J1, P2, n2, dP2, ddP2, J2, u_ex, du_dn_ex);

        % 2.4. Розв'язок
        if N < 8, psi = pinv(A)*RHS; else, psi = A\RHS; end
        p1 = psi(1:N);
        p2 = psi(N+1:end);

        % Збереження
        if N == N_vals(end)
            res.p1=p1; res.p2=p2;
            res.P1=P1; res.n1=n1; res.J1=J1;
            res.P2=P2; res.n2=n2; res.J2=J2;
        end

        % 2.5. Похибка
        errs = zeros(1,3);
        for k=1:3
            u_calc = compute_u(pts(k,:), p1, p2, P1, n1, J1, P2, J2);
            errs(k) = abs(u_calc - u_ex(pts(k,:)));
        end
        fprintf(' %-2d | %.3e  | %.3e  | %.3e \n', N, errs);
    end
    fprintf('------------------------------------------------------------\n');

    % --- 3. ВІЗУАЛІЗАЦІЯ ---
    plot_results(geom, pts, y_star, res);
end

% =========================================================================
% БЛОК 1: ГЕОМЕТРІЯ
% =========================================================================
function [P, n, dP, ddP, J] = compute_geom(t, fx, fy, fdx, fdy, fddx, fddy, is_inner)
    % Обчислює координати, похідні, якобіан та нормаль
    x = fx(t); y = fy(t);
    dx = fdx(t); dy = fdy(t);
    ddx = fddx(t); ddy = fddy(t);

    P = [x, y];
    dP = [dx, dy];
    ddP = [ddx, ddy];
    J = sqrt(dx.^2 + dy.^2);

    % Нормаль (y', -x') / J
    n = [dy, -dx] ./ J;

    % Корекція для внутрішньої границі (має дивитися в дірку)
    if is_inner
        n = -n;
    end
end

% =========================================================================
% БЛОК 2: ЗБІРКА СИСТЕМИ (ЯДРА)
% =========================================================================
function [A, RHS] = assemble_system(N, P1, n1, dP1, ddP1, J1, P2, n2, dP2, ddP2, J2, u_ex, du_dn_ex)
    Dim = 2*N;
    A = zeros(Dim);
    RHS = zeros(Dim, 1);
    w = 1/N; % Вага квадратури

    for i = 1:N
        for j = 1:N
            % --- Блок 1,1: L11 (Г1 -> Г1) ---
            % Потенціал подвійного шару. Особливість -> Кривина
            if i == j
                curv = dot(ddP1(i,:), n1(i,:));
                val = curv / (2 * J1(i));
            else
                r = P1(i,:) - P1(j,:);
                val = dot(r, n1(j,:)) / sum(r.^2) * J1(j);
            end
            A(i, j) = val * w;
            if i == j, A(i, j) = A(i, j) - 0.5; end

            % --- Блок 1,2: L12 (Г2 -> Г1) ---
            % Потенціал простого шару. Регулярне.
            r = P1(i,:) - P2(j,:);
            val = log(1 / norm(r)) * J2(j);
            A(i, j+N) = val * w;

            % --- Блок 2,1: L21 (Г1 -> Г2) ---
            % Мішана похідна (Гіперсингулярна формула). Регулярне.
            vec_r = P2(i,:) - P1(j,:);
            d2 = sum(vec_r.^2);
            term1 = dot(n2(i,:), n1(j,:));
            term2 = 2 * dot(vec_r, n2(i,:)) * dot(vec_r, n1(j,:)) / d2;
            val = -(term1 - term2) / d2 * J1(j);
            A(i+N, j) = val * w;

            % --- Блок 2,2: L22 (Г2 -> Г2) ---
            % Похідна простого шару. Особливість -> Кривина
            if i == j
                curv = dot(ddP2(i,:), n2(i,:));
                val = curv / (2 * J2(i));
            else
                r = P2(j,:) - P2(i,:); % Увага: y-x для похідної по x
                val = dot(r, n2(i,:)) / sum(r.^2) * J2(j);
            end
            A(i+N, j+N) = val * w;
            if i == j, A(i+N, j+N) = A(i+N, j+N) + 0.5; end
        end

        % Праві частини
        RHS(i) = u_ex(P1(i,:));
        RHS(i+N) = du_dn_ex(P2(i,:), n2(i,:));
    end
end

% =========================================================================
% БЛОК 3: ВІДНОВЛЕННЯ РОЗВ'ЯЗКУ
% =========================================================================
function u = compute_u(pt, p1, p2, P1, n1, J1, P2, J2)
    N = length(p1);

    % Внесок від Г1 (Подвійний шар)
    R1 = pt - P1;
    term1 = sum(R1 .* n1, 2) ./ sum(R1.^2, 2);
    int1 = sum(p1 .* term1 .* J1);

    % Внесок від Г2 (Простий шар)
    R2 = pt - P2;
    term2 = log(1 ./ sqrt(sum(R2.^2, 2)));
    int2 = sum(p2 .* term2 .* J2);

    % Сума
    u = (1/N) * (int1 + int2);
end

% =========================================================================
% БЛОК 4: ГРАФІКА
% =========================================================================
function plot_results(geom, pts, y_star, res)
    % 1. ГЕОМЕТРІЯ (ЛІВИЙ ГРАФІК)
    figure(1); clf;
    hold on; axis equal; grid on; box on; set(gca, 'FontSize', 12);

    tt = linspace(0, 2*pi, 200);
    % Малюємо межі
    plot(geom.x1(tt), geom.y1(tt), 'b-', 'LineWidth', 2);
    plot(geom.x2(tt), geom.y2(tt), 'r-', 'LineWidth', 2);

    % Малюємо точки
    plot(pts(:,1), pts(:,2), 'ko', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
    plot(y_star(1), y_star(2), 'm*', 'MarkerSize', 10);

    legend('\Gamma_1 (Внутр.)', '\Gamma_2 (Зовн.)', 'Точки', 'Джерело');
    title('Геометрія задачі'); xlabel('x'); ylabel('y');

    % 2. ЛІНІЇ РІВНЯ / РОЗПОДІЛ ПОЛЯ (ПРАВИЙ ГРАФІК)
    figure(2); clf;

    % Створюємо густу сітку для красивої картинки
    [X, Y] = meshgrid(linspace(-3.5, 3.5, 120));
    U = nan(size(X));

    fprintf('Побудова графіку ліній рівня... ');
    for k=1:numel(X)
        pt = [X(k), Y(k)];
        r = norm(pt);
        % Перевіряємо, чи точка всередині кільця (1 < r < 3)
        % Даємо маленький відступ від меж, щоб не було артефактів
        if r > 1.02 && r < 2.98
            U(k) = compute_u(pt, res.p1, res.p2, res.P1, res.n1, res.J1, res.P2, res.J2);
        end
    end
    fprintf('Готово.\n');

    % Використовуємо surf для стабільності + view(2) для 2D вигляду
    surf(X, Y, real(U));
    view(2);            % Вигляд зверху
    shading interp;     % Плавні переходи кольорів (як на фото)
    axis equal; axis tight;
    colorbar;           % Шкала кольорів
    colormap jet;       % Колірна схема "jet" (синій-червоний)

    % Домальовуємо границі поверх (білим або чорним кольором)
    hold on;
    plot(geom.x1(tt), geom.y1(tt), 'k-', 'LineWidth', 1.5);
    plot(geom.x2(tt), geom.y2(tt), 'k-', 'LineWidth', 1.5);

    title('Розподіл поля u(x) (Наближений розв''язок)');
    xlabel('x'); ylabel('y');
end
