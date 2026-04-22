function solve_all_examples_full()
% =========================================================================
% ПРОГРАМА ДЛЯ ЧИСЕЛЬНОГО РОЗВ'ЯЗУВАННЯ ЗАДАЧІ ДІРІХЛЕ-НЕЙМАНА
% Метод: Нистрьом для системи інтегральних рівнянь 2-го роду
% =========================================================================
    clc;
    close all;

    % --- МЕНЮ ВИБОРУ ПРИКЛАДУ ---
    fprintf('============================================================\n');
    fprintf('Оберіть варіант задачі для розрахунку:\n');
    fprintf('------------------------------------------------------------\n');
    fprintf('1 - Приклад №1: Дві концентричні кола (R=1, R=3)\n');
    fprintf('2 - Приклад №2: Еліпс всередині кола\n');
    fprintf('3 - Приклад №3: "Сплюснутий змій" всередині еліпса\n');
    fprintf('------------------------------------------------------------\n');

    choice_str = input('Ваш вибір (1-3): ', 's');
    choice = str2double(choice_str);
    if isnan(choice) || choice < 1 || choice > 3
        choice = 1;
        fprintf('Вибрано за замовчуванням: Приклад №1\n');
    end

    % --- ПАРАМЕТРИ МЕТОДУ ---
    N_values = [4, 8, 16, 32, 64]; % Кількість точок дискретизації

    % --- ТОЧНИЙ РОЗВ'ЯЗОК (ДЛЯ ПЕРЕВІРКИ) ---
    % Використовуємо фундаментальний розв'язок з джерелом поза областю
    y_star = [5.0, 5.0];

    fprintf('\n>>> ЗАПУСК РОЗРАХУНКУ <<<\n');
    fprintf('Точний розв''язок: u(x) = 1/(2pi) * ln(1/|x - y*|)\n');
    fprintf('Джерело y*: (%g, %g)\n', y_star(1), y_star(2));

    % --- ТЕСТОВІ ТОЧКИ ---
    % Вибираємо точки, які гарантовано лежать всередині області D
    % для всіх трьох прикладів (приблизно посередині між межами)
    test_points = [
        2.0,  0.0;
        0.0,  1.8;
       -1.5, -0.5
    ];

    fprintf('\nТАБЛИЦЯ АБСОЛЮТНИХ ПОХИБОК:\n');
    fprintf('------------------------------------------------------------\n');
    fprintf(' N  |  Err(Pt 1)   |  Err(Pt 2)   |  Err(Pt 3) \n');
    fprintf('------------------------------------------------------------\n');

    % Структура для зберігання даних останньої ітерації (для графіків)
    ResultData = struct();

    % --- ОСНОВНИЙ ЦИКЛ ПО N ---
    for k = 1:length(N_values)
        N = N_values(k);

        % 1. Побудова рівномірної сітки
        h = 2 * pi / N;
        t = (0 : N-1)' * h;

        % 2. Обчислення геометрії (координати, похідні, нормалі)
        % Г1 (Внутрішня)
        [P1, dP1, ddP1, nu1, jac1] = compute_geometry(t, 1, choice);

        % Г2 (Зовнішня)
        [P2, dP2, ddP2, nu2, jac2] = compute_geometry(t, 2, choice);

        % 3. Формування матриці системи та вектора правих частин
        [MatrixA, VectorRHS] = assemble_linear_system(N, P1, nu1, dP1, ddP1, jac1, ...
                                                         P2, nu2, dP2, ddP2, jac2, ...
                                                         y_star);

        % 4. Розв'язування системи лінійних рівнянь
        % Використовуємо псевдообернення для малих N (стійкість)
        if N < 8
            Psi_vector = pinv(MatrixA) * VectorRHS;
        else
            Psi_vector = MatrixA \ VectorRHS;
        end

        % Розбиваємо вектор розв'язку на дві частини
        psi1 = Psi_vector(1:N);
        psi2 = Psi_vector(N+1:end);

        % Збереження результатів для N=64
        if N == 64
            ResultData.psi1 = psi1;
            ResultData.psi2 = psi2;
            ResultData.P1 = P1;
            ResultData.P2 = P2;
            ResultData.nu1 = nu1;
            ResultData.jac1 = jac1;
            ResultData.jac2 = jac2;
        end

        % 5. Перевірка точності в тестових точках
        errors = zeros(1, 3);
        for p_idx = 1:3
            pt = test_points(p_idx, :);

            % Обчислене значення
            u_calc = compute_field_at_point(pt, psi1, psi2, P1, nu1, jac1, P2, jac2);

            % Точне значення
            u_true = exact_solution(pt, y_star);

            errors(p_idx) = abs(u_calc - u_true);
        end

        % Виведення рядка таблиці
        fprintf(' %-2d | %.4e   | %.4e   | %.4e \n', N, errors);
    end
    fprintf('------------------------------------------------------------\n');

    % --- ВІЗУАЛІЗАЦІЯ РЕЗУЛЬТАТІВ ---
    plot_results(ResultData, test_points, y_star, choice);
end


% =========================================================================
% ФУНКЦІЯ: ЗАДАННЯ ГЕОМЕТРІЇ (ПАРАМЕТРИЗАЦІЯ)
% =========================================================================
function [P, dP, ddP, nu, jac] = compute_geometry(t, boundary_id, example_id)
    % boundary_id: 1 (Внутрішня Г1), 2 (Зовнішня Г2)
    % example_id: Номер прикладу

    x = zeros(size(t)); y = zeros(size(t));
    dx = zeros(size(t)); dy = zeros(size(t));
    ddx = zeros(size(t)); ddy = zeros(size(t));

    % --- ПРИКЛАД 1: ДВА КОЛА (R=1, R=3) ---
    if example_id == 1
        if boundary_id == 1 % Г1 (Внутрішня)
            R = 1.0;
        else                % Г2 (Зовнішня)
            R = 3.0;
        end
        x = R * cos(t);       dx = -R * sin(t);      ddx = -R * cos(t);
        y = R * sin(t);       dy =  R * cos(t);      ddy = -R * sin(t);

    % --- ПРИКЛАД 2: ЕЛІПС ВСЕРЕДИНІ КОЛА ---
    elseif example_id == 2
        if boundary_id == 1 % Г1: Еліпс
            a = 1.5; b = 0.7;
            x = a * cos(t);   dx = -a * sin(t);    ddx = -a * cos(t);
            y = b * sin(t);   dy =  b * cos(t);    ddy = -b * sin(t);
        else                % Г2: Коло R=3
            R = 3.0;
            x = R * cos(t);   dx = -R * sin(t);    ddx = -R * cos(t);
            y = R * sin(t);   dy =  R * cos(t);    ddy = -R * sin(t);
        end

    % --- ПРИКЛАД 3: СПЛЮСНУТИЙ ЗМІЙ (АСТРОЇДА) В ЕЛІПСІ ---
    elseif example_id == 3
        if boundary_id == 1 % Г1: Модифікована Астроїда (Змій)
            % Формула: x = a*cos^3(t), y = b*sin^3(t) (згладжена)
            % Для гладкості додамо лінійні члени, щоб не було гострих кутів
            % x = cos(t) + 0.3*cos(3t)
            % y = sin(t) - 0.3*sin(3t)

            x = cos(t) + 0.3 * cos(3*t);
            y = sin(t) - 0.3 * sin(3*t);

            dx = -sin(t) - 0.9 * sin(3*t);
            dy =  cos(t) - 0.9 * cos(3*t);

            ddx = -cos(t) - 2.7 * cos(3*t);
            ddy = -sin(t) + 2.7 * sin(3*t);

        else                % Г2: Еліпс (Витягнутий)
            a_out = 4.0; b_out = 3.0;
            x = a_out * cos(t);   dx = -a_out * sin(t);    ddx = -a_out * cos(t);
            y = b_out * sin(t);   dy =  b_out * cos(t);    ddy = -b_out * sin(t);
        end
    end

    % Формування вихідних масивів
    P = [x, y];
    dP = [dx, dy];
    ddP = [ddx, ddy];

    % Якобіан (довжина дотичного вектора)
    jac = sqrt(dx.^2 + dy.^2);

    % Нормаль
    % Стандартна формула (dy, -dx) дає вектор праворуч від руху.
    % При обході проти годинникової стрілки це вектор НАЗОВНІ кривої.
    nu = [dy, -dx] ./ jac;

    % ВАЖЛИВО: Для внутрішньої границі (Г1) зовнішня нормаль до області D
    % повинна дивитися ВНУТР кривої (в дірку). Тому змінюємо знак.
    if boundary_id == 1
        nu = -nu;
    end
end


% =========================================================================
% ФУНКЦІЯ: ТОЧНИЙ РОЗВ'ЯЗОК ТА ЙОГО ПОХІДНА
% =========================================================================
function val = exact_solution(x, y_star)
    % u(x) = 1/(2pi) * ln(1 / |x - y*|)
    dist = norm(x - y_star);
    val = (1 / (2 * pi)) * log(1 / dist);
end

function val = exact_normal_derivative(x, nu, y_star)
    % du/dn = (grad u, nu)
    % grad u = -1/(2pi) * (x - y*) / |x - y*|^2
    diff_vec = x - y_star;
    dist_sq = sum(diff_vec.^2);
    grad_u = (-1 / (2 * pi)) * diff_vec / dist_sq;
    val = dot(grad_u, nu);
end


% =========================================================================
% ФУНКЦІЯ: ЗБІРКА СИСТЕМИ ЛІНІЙНИХ РІВНЯНЬ
% =========================================================================
function [A, RHS] = assemble_linear_system(N, P1, n1, dP1, ddP1, J1, ...
                                              P2, n2, dP2, ddP2, J2, y_star)
    Dim = 2 * N;
    A = zeros(Dim, Dim);
    RHS = zeros(Dim, 1);

    % Вага квадратурної формули трапецій:
    % Інтеграл (1/2pi) * int(0..2pi) ~ (1/2pi) * (2pi/2N) * sum = (1/2N) * sum
    weight = 1 / (2 * N);

    for i = 1:N
        for j = 1:N
            % -------------------------------------------------------------
            % Блок A11: Г1 -> Г1 (Потенціал подвійного шару)
            % Ядро L11. На діагоналі - кривина.
            % -------------------------------------------------------------
            if i == j
                numerator = dot(ddP1(i,:), n1(i,:));
                val_L11 = numerator / (2 * J1(i));
            else
                r_vec = P1(i,:) - P1(j,:);
                dist2 = sum(r_vec.^2);
                val_L11 = dot(r_vec, n1(j,:)) / dist2 * J1(j);
            end

            A(i, j) = val_L11 * weight;

            % Додаємо стрибок потенціалу на діагоналі (-0.5 * I)
            if i == j
                A(i, j) = A(i, j) - 0.5;
            end

            % -------------------------------------------------------------
            % Блок A12: Г2 -> Г1 (Потенціал простого шару)
            % Ядро L12 = ln(1/r)
            % -------------------------------------------------------------
            r_vec = P1(i,:) - P2(j,:);
            dist = norm(r_vec);
            val_L12 = log(1 / dist) * J2(j);

            A(i, j + N) = val_L12 * weight;

            % -------------------------------------------------------------
            % Блок A21: Г1 -> Г2 (Мішана похідна / Гіперсингулярне)
            % Ядро L21
            % -------------------------------------------------------------
            r_vec = P2(i,:) - P1(j,:);
            dist2 = sum(r_vec.^2);

            term1 = dot(n2(i,:), n1(j,:));
            term2 = 2 * dot(r_vec, n2(i,:)) * dot(r_vec, n1(j,:)) / dist2;

            % Формула: - (n(x)n(y) - 2...)/r^2
            val_L21 = -(term1 - term2) / dist2 * J1(j);

            A(i + N, j) = val_L21 * weight;

            % -------------------------------------------------------------
            % Блок A22: Г2 -> Г2 (Похідна простого шару)
            % Ядро L22. На діагоналі - кривина.
            % -------------------------------------------------------------
            if i == j
                numerator = dot(ddP2(i,:), n2(i,:));
                val_L22 = numerator / (2 * J2(i));
            else
                r_vec = P2(j,:) - P2(i,:); % Вектор y - x
                dist2 = sum(r_vec.^2);
                val_L22 = dot(r_vec, n2(i,:)) / dist2 * J2(j);
            end

            A(i + N, j + N) = val_L22 * weight;

            % Додаємо стрибок похідної на діагоналі (+0.5 * I)
            if i == j
                A(i + N, j + N) = A(i + N, j + N) + 0.5;
            end
        end

        % --- ПРАВА ЧАСТИНА ---
        % Г1 (Діріхле): значення u_exact
        RHS(i) = exact_solution(P1(i,:), y_star);

        % Г2 (Нейман): значення нормальної похідної
        RHS(i + N) = exact_normal_derivative(P2(i,:), n2(i,:), y_star);
    end
end


% =========================================================================
% ФУНКЦІЯ: ВІДНОВЛЕННЯ ПОЛЯ В ТОЧЦІ
% =========================================================================
function u_val = compute_field_at_point(pt, psi1, psi2, P1, nu1, jac1, P2, jac2)
    N = length(psi1);
    weight = 1 / (2 * N);

    % Внесок від Г1 (Потенціал подвійного шару)
    R1 = pt - P1;
    dist1_sq = sum(R1.^2, 2);
    kernel1 = sum(R1 .* nu1, 2) ./ dist1_sq;
    integral1 = sum(psi1 .* kernel1 .* jac1);

    % Внесок від Г2 (Потенціал простого шару)
    R2 = pt - P2;
    dist2_sq = sum(R2.^2, 2);
    kernel2 = -0.5 * log(dist2_sq); % це ln(1/r)
    integral2 = sum(psi2 .* kernel2 .* jac2);

    u_val = weight * (integral1 + integral2);
end


% =========================================================================
% ФУНКЦІЯ: ВІЗУАЛІЗАЦІЯ РЕЗУЛЬТАТІВ
% =========================================================================
function plot_results(Data, test_points, y_star, example_id)

    % --- РИС. 1: ГЕОМЕТРІЯ ЗАДАЧІ ---
    figure(1); clf;
    hold on; axis equal; grid on; box on;
    set(gca, 'FontSize', 12);

    % Малюємо гладкі лінії для меж
    tt = linspace(0, 2*pi, 300);
    [p1_plt, ~, ~, ~, ~] = compute_geometry(tt', 1, example_id);
    [p2_plt, ~, ~, ~, ~] = compute_geometry(tt', 2, example_id);

    plot(p1_plt(:,1), p1_plt(:,2), 'b-', 'LineWidth', 2);
    plot(p2_plt(:,1), p2_plt(:,2), 'r-', 'LineWidth', 2);

    % Точки
    plot(test_points(:,1), test_points(:,2), 'ko', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
    plot(y_star(1), y_star(2), 'm*', 'MarkerSize', 10, 'LineWidth', 2);

    legend('\Gamma_1 (Внутрішня)', '\Gamma_2 (Зовнішня)', 'Тестові точки', 'Джерело y^*');
    title('Геометрія розрахункової області');
    xlabel('x'); ylabel('y');


    % --- ПІДГОТОВКА ДАНИХ ДЛЯ 3D ГРАФІКІВ ---
    fprintf('Побудова графіків... ');

    % Визначаємо межі для графіків
    x_min = min(p2_plt(:,1)) - 0.5; x_max = max(p2_plt(:,1)) + 0.5;
    y_min = min(p2_plt(:,2)) - 0.5; y_max = max(p2_plt(:,2)) + 0.5;

    n_grid = 70;
    [X, Y] = meshgrid(linspace(x_min, x_max, n_grid), linspace(y_min, y_max, n_grid));
    U_grid = nan(size(X));

    for k = 1:numel(X)
        pt = [X(k), Y(k)];

        % Перевірка приналежності до області D (ray casting algorithm)
        in_outer = inpolygon(pt(1), pt(2), p2_plt(:,1), p2_plt(:,2));
        in_inner = inpolygon(pt(1), pt(2), p1_plt(:,1), p1_plt(:,2));

        % Додатковий відступ від межі для краси (щоб не було артефактів)
        dist_to_outer = min(sqrt(sum((p2_plt - pt).^2, 2)));
        dist_to_inner = min(sqrt(sum((p1_plt - pt).^2, 2)));

        if in_outer && ~in_inner && dist_to_outer > 0.05 && dist_to_inner > 0.05
            U_grid(k) = compute_field_at_point(pt, Data.psi1, Data.psi2, ...
                                               Data.P1, Data.nu1, Data.jac1, ...
                                               Data.P2, Data.jac2);
        end
    end
    fprintf('Готово.\n');

    % --- РИС. 2: 3D MESH PLOT ---
    figure(2); clf;
    mesh(X, Y, real(U_grid));
    view(-35, 45);
    axis tight;
    colormap jet;
    title('Графік наближеного розв''язку (Mesh)');
    xlabel('x'); ylabel('y'); zlabel('u(x)');
    grid on;

    % --- РИС. 3: CONTOUR PLOT ---
    figure(3); clf; hold on; axis equal; box on;
    set(gca, 'FontSize', 12);

    % Лінії рівня
    [C, h] = contour(X, Y, real(U_grid), 20, 'LineWidth', 1.5);
    colorbar;
    colormap jet;

    % Межі області (поверх контурів)
    plot(p1_plt(:,1), p1_plt(:,2), 'k-', 'LineWidth', 2);
    plot(p2_plt(:,1), p2_plt(:,2), 'k-', 'LineWidth', 2);

    title('Лінії рівня розв''язку');
    xlabel('x'); ylabel('y');
end
