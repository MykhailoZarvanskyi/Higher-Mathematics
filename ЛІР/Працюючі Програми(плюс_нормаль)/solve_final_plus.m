function solve_final_plus()

    clc; close all;
    fprintf('Запуск програми \n');

    % 1. Параметризація кривих


    % Г1 (Внутрішня межа, R=1, умова Діріхле)
    % Параметричне рівняння кола: x(t) = (cos t, sin t)
    geom.x1 = @(t) cos(t);   geom.y1 = @(t) sin(t);

    % Перші похідні (потрібні для дотичної, нормалі та якобіана):
    % x'(t) = -sin t, y'(t) = cos t
    geom.dx1 = @(t) -sin(t); geom.dy1 = @(t) cos(t);

    % Другі похідні (потрібні ТІЛЬКИ для обчислення кривини на діагоналі):
    % x''(t) = -cos t, y''(t) = -sin t
    geom.ddx1 = @(t) -cos(t); geom.ddy1 = @(t) -sin(t);

    % Г2 (Зовнішня межа, R=3, умова Неймана)
    % Параметричне рівняння кола: x(t) = (3 cos t, 3 sin t)
    geom.x2 = @(t) 3*cos(t);  geom.y2 = @(t) 3*sin(t);

    % Перші похідні:
    geom.dx2 = @(t) -3*sin(t); geom.dy2 = @(t) 3*cos(t);

    % Другі похідні:
    geom.ddx2 = @(t) -3*cos(t); geom.ddy2 = @(t) -3*sin(t);


    % 2. Точний розв'язок


    % Джерело y* обираємо поза областю D
    y_star = [5, 5];

    % Функція точного розв'язку: фундаментальний розв'язок рівняння Лапласа
    % Формула: u(x) = 1/(2pi) * ln(1 / |x - y*|)
    u_exact = @(x) (1/(2*pi)) * log(1 / norm(x - y_star));

    % Точна нормальна похідна (для умови Неймана)
    % Формула: du/dn = (grad u, n) = -1/(2pi) * ((x - y*) * n) / |x - y*|^2
    du_dn_exact = @(x, nu) -1/(2*pi) * dot(x - y_star, nu) / sum((x - y_star).^2);


    % 3. Тестові точки

    % Точки лежать всередині кільця (між R=1 і R=3)
    test_points = [2.0, 0.0; 1.0, 1.0; 0.0, -2.0];


    % 4. Серце


    % N - параметр дискретизації з таблиці
    % Реальна кількість вузлів на одній кривій M = 2*N
    N_values = [4, 8, 16, 32, 64];


    fprintf('\n !ТАБЛИЦЯ ПОХИБОК!:\n');
    fprintf('------------------------------------------------------------\n');
    fprintf(' N  | 2N (Вузлів) |  Err(2,0)   |  Err(1,1)   |  Err(0,-2) \n');
    fprintf('------------------------------------------------------------\n');

    res = struct(); % Структура для зберігання результатів

    for N = N_values
        % 4.1: Підготовка сітки

        M = 2 * N;          % Загальна кількість вузлів на одній кривій
        h = 2*pi / M;       % Крок сітки по параметру t
        t = (0 : M-1)' * h; % Вектор значень параметра t_j = j*h

        % Вага квадратури трапецій:
        % Інтеграл 1/(2pi) * int(...) dt наближається сумою: (1/2pi) * (2pi/M) * sum(...) = (1/M) * sum(...)
        w = 1 / M;

        % 4.2: Обчислення геометрії у вузлах

        % Геометрія Г1 (координати P1, похідні dP1/ddP1, нормаль nu1, якобіан jac1)
        [P1, dP1, ddP1, nu1, jac1] = get_geom(t, geom.x1, geom.y1, geom.dx1, geom.dy1, geom.ddx1, geom.ddy1);

        % !!! Стандартна формула дає нормаль назовні кола!!!
        % Для внутрішньої межі (Г1) зовнішня нормаль до області D має дивитися НАЗОВНІ, тому змінюємо знак вектора нормалі!
        nu1 = nu1;

        % Геометрія Г2 (тут нормаль правильна - ззовні від D)
        [P2, dP2, ddP2, nu2, jac2] = get_geom(t, geom.x2, geom.y2, geom.dx2, geom.dy2, geom.ddx2, geom.ddy2);

        % 4.3: Формування СЛАР

        Dim = 2 * M;         % Розмірність матриці (2 криві * M точок)
        A = zeros(Dim, Dim); % Матриця системи
        RHS = zeros(Dim, 1); % Вектор правої частини

        % Подвійний цикл по точках спостереження (i) та інтегрування (j)
        for i = 1:M
            for j = 1:M

                % БЛОК A11: Г1 -> Г1 (Потенціал подвійного шару)
                if i == j
                    % Діагональ: Тут усувна особливість, тому формула: (x'' * nu) / (2 * |x'|)
                    val = dot(ddP1(i,:), nu1(i,:)) / (2 * jac1(i));
                else
                    % Поза діагоналі:
                    % Формула: (x - y) * nu(y) / |x - y|^2 * jac(y)
                    r = P1(i,:) - P1(j,:);
                    val = dot(r, nu1(j,:)) / sum(r.^2) * jac1(j);
                end
                A(i, j) = val * w; % Множимо на вагу квадратури

                % !Стрибок!: Для рівняння на діагоналі додається доданок +0.5 * psi (через перехід на межу)
                if i == j, A(i, j) += 0.5; end

                % БЛОК A12: Г2 -> Г1 (Потенціал простого шару)
                % Тут немає особливостей, тому ядро: ln(1/|x - y|) * jac(y)
                r = P1(i,:) - P2(j,:);
                A(i, j+M) = log(1 / norm(r)) * jac2(j) * w;

                % БЛОК A21: Г1 -> Г2 (Мішана похідна - найскладніша)
                % Формула: -( (nu_x * nu_y) - 2((r*nu_x)*(r*nu_y))/r^2 ) / r^2 * jac(y)
                r_vec = P2(i,:) - P1(j,:);
                d2 = sum(r_vec.^2);

                t1 = dot(nu2(i,:), nu1(j,:));
                t2 = 2 * dot(r_vec, nu2(i,:)) * dot(r_vec, nu1(j,:)) / d2;

                val = -(t1 - t2) / d2 * jac1(j);
                A(i+M, j) = val * w;

                % БЛОК A22: Г2 -> Г2 (Похідна простого шару)
                if i == j
                    % Діагональ: Та сама, що і в A11
                    val = dot(ddP2(i,:), nu2(i,:)) / (2 * jac2(i));
                else
                    % Поза діагоналю: (y - x) * nu(x) / |x - y|^2 * jac(y)
                    % !Тут вектор r = P2(j,:) - P2(i,:) (це y - x)!
                    r = P2(j,:) - P2(i,:);
                    val = dot(r, nu2(i,:)) / sum(r.^2) * jac2(j);
                end
                A(i+M, j+M) = val * w;

                % Стрибок: Для похідної простого шару стрибок +0.5
                if i == j, A(i+M, j+M) += 0.5; end
            end

            % Права частина
            % На Г1 задана умова Діріхле: u = u_exact
            RHS(i) = u_exact(P1(i,:));
            % На Г2 задана умова Неймана: du/dn = du_dn_exact
            RHS(i+M) = du_dn_exact(P2(i,:), nu2(i,:));
        end

        % Крок 4.4: Розв'язування системи
        if M < 8
             Psi = pinv(A)*RHS; % Псевдообернення для стабільності при дуже малих N (Підказав ШІ)
        else
             Psi = A\RHS;
        end

        % Розбиваємо вектор розв'язку на дві частини для Г1 і Г2
        psi1 = Psi(1:M);
        psi2 = Psi(M+1:end);

        % Збереження результатів останньої ітерації для графіків
        if N == 64
            res.psi1 = psi1; res.psi2 = psi2;
            res.P1 = P1; res.P2 = P2; res.nu1 = nu1; res.jac1 = jac1; res.jac2 = jac2;
        end

        % Крок 4.5: Перевірка похибки
        errs = zeros(1, 3);
        for k = 1:3
            % Обчислюємо наближений розв'язок у тестовій точці
            u_val = get_u(test_points(k,:), psi1, psi2, P1, nu1, jac1, P2, jac2);
            % Порівнюємо з точним
            errs(k) = abs(u_val - u_exact(test_points(k,:)));
        end


        fprintf(' %-2d | %-10d  | %.3e  | %.3e  | %.3e \n', N, M, errs);
    end
    fprintf('------------------------------------------------------------\n');


    % --- 5. ВІЗУАЛІЗАЦІЯ ---

    % Графік 1: Геометрія
    figure(1); clf; hold on; axis equal; grid on; box on; set(gca, 'FontSize', 12);
    tt = linspace(0, 2*pi, 200);

    % Обчислюємо координати меж для малювання
    p1_draw = [geom.x1(tt); geom.y1(tt)]';
    p2_draw = [geom.x2(tt); geom.y2(tt)]';

    plot(p1_draw(:,1), p1_draw(:,2), 'b-', 'LineWidth', 2);
    plot(p2_draw(:,1), p2_draw(:,2), 'r-', 'LineWidth', 2);
    plot(test_points(:,1), test_points(:,2), 'ko', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
    plot(y_star(1), y_star(2), 'm*', 'MarkerSize', 10);
    title('Геометрія'); legend('Г_1', 'Г_2', 'Points', 'Source');

    % Підготовка даних для полів
    n_grid = 80; % Достатньо густа сітка
    bound = 3.5;
    [X, Y] = meshgrid(linspace(-bound, bound, n_grid));
    Q = nan(size(X));

    fprintf('Побудова графіків...\n');
    for k=1:numel(X)
        pt = [X(k), Y(k)];
        r = norm(pt);

        % Перевірка приналежності: строго між колами R=1 і R=3
        % Додаємо невеликий відступ (epsilon = 0.05), щоб прибрати артефакти на краях
        if r > 1.05 && r < 2.95
            Q(k) = get_u(pt, res.psi1, res.psi2, res.P1, res.nu1, res.jac1, res.P2, res.jac2);
        end
    end

    % Графік 2: MESH (3D Сітка)
    figure(2); clf;
    mesh(X, Y, real(Q));
    view(-30, 60);
    axis([-bound bound -bound bound]);
    colormap jet;
    title('Графік наближеного розв''язку (Mesh)');
    xlabel('x'); ylabel('y'); zlabel('u(x)');

    % Графік 3: CONTOUR (Лінії рівня)
    figure(3); clf; hold on; axis equal; box on;
    [C, h] = contour(X, Y, real(Q), 25, 'LineWidth', 1.5); % 25 ліній рівня

    % Домальовуємо межі поверх
    plot(p1_draw(:,1), p1_draw(:,2), 'k-', 'LineWidth', 2);
    plot(p2_draw(:,1), p2_draw(:,2), 'k-', 'LineWidth', 2);

    colorbar; colormap jet;
    axis([-bound bound -bound bound]);
    title('Лінії рівня');
end


% Допоміжні функції


% Функція обчислення геометрії: повертає точки P, похідні dP, ddP, нормаль nu, якобіан jac
function [P, dP, ddP, nu, jac] = get_geom(t, fx, fy, fdx, fdy, fddx, fddy)
    N = length(t);
    P = [fx(t), fy(t)];         % Координати (x, y)
    dP = [fdx(t), fdy(t)];      % Перші похідні (x', y')
    ddP = [fddx(t), fddy(t)];   % Другі похідні (x'', y'')

    jac = sqrt(sum(dP.^2, 2));  % Якобіан = sqrt((x')^2 + (y')^2)
    nu = [dP(:,2), -dP(:,1)] ./ jac; % Нормаль (y', -x')/J
end

% Функція відновлення розв'язку
function u = get_u(pt, psi1, psi2, P1, nu1, jac1, P2, jac2)
    M = length(psi1); % Кількість точок (це 2N)
    w = 1 / M;        % Вага 1/(2N)

    % Доданок 1: Інтеграл по Г1 (Потенціал подвійного шару)
    % Ядро L1 = (x - y) * nu(y) / |x - y|^2
    R1 = pt - P1;
    dist1_sq = sum(R1.^2, 2);
    term1 = sum(R1 .* nu1, 2) ./ dist1_sq;
    % Інтегральна сума: sum(psi1 * Kernel * Jacobian)
    int1 = sum(psi1 .* term1 .* jac1);

    % Доданок 2: Інтеграл по Г2 (Потенціал простого шару)
    % Ядро L2 = ln(1 / |x - y|)
    R2 = pt - P2;
    dist2_sq = sum(R2.^2, 2);
    term2 = 0.5 * log(1 ./ dist2_sq); % Це те саме, що log(1/r)
    int2 = sum(psi2 .* term2 .* jac2);

    % Загальна формула: u = w * (int1 + int2)
    u = w * (int1 + int2);
end
