function solve_dirichlet_neumann()
% =========================================================================
% ПРОГРАМА: Розв'язування задачі Діріхле-Неймана (Метод Нистрьома)
% Покращена версія з автоматичними похідними, обробкою сингулярностей
% =========================================================================
    clc; close all;
    fprintf('============================================================\n');
    fprintf('   ЧИСЕЛЬНЕ РОЗВ''ЯЗУВАННЯ МІШАНОЇ ЗАДАЧІ (МЕТОД НИСТРЬОМА)\n');
    fprintf('============================================================\n');
    fprintf('ВАЖЛИВО:\n');
    fprintf('1. Вводьте формули для t від 0 до 2*pi.\n');
    fprintf('2. Похідні обчислюються автоматично.\n');
    fprintf('3. Програма сама виправить напрямок нормалі для внутрішньої границі.\n\n');

    % --- 1. ВВЕДЕННЯ КРИВИХ З АВТОМАТИЧНИМИ ПОХІДНИМИ ---
    
    % Г1 - Внутрішня крива (Діріхле)
    fprintf('>>> Введіть параметри ВНУТРІШНЬОЇ кривої Г1 (Діріхле) <<<\n');
    def_x1 = 'cos(t)';
    def_y1 = 'sin(t)';
    s_x1 = input(sprintf('x1(t) [Enter = "%s"]: ', def_x1), 's');
    if isempty(s_x1), s_x1 = def_x1; end
    s_y1 = input(sprintf('y1(t) [Enter = "%s"]: ', def_y1), 's');
    if isempty(s_y1), s_y1 = def_y1; end

    % Створюємо анонімні функції та автоматично обчислюємо похідні
    fprintf('   Обчислення похідних для Г1...\n');
    [geom.x1, geom.y1, geom.dx1, geom.dy1, geom.ddx1, geom.ddy1] = ...
        create_curve_with_derivatives(s_x1, s_y1);

    % Г2 - Зовнішня крива (Нейман)
    fprintf('\n>>> Введіть параметри ЗОВНІШНЬОЇ кривої Г2 (Нейман) <<<\n');
    def_x2 = '3*cos(t)';
    def_y2 = '3*sin(t)';
    s_x2 = input(sprintf('x2(t) [Enter = "%s"]: ', def_x2), 's');
    if isempty(s_x2), s_x2 = def_x2; end
    s_y2 = input(sprintf('y2(t) [Enter = "%s"]: ', def_y2), 's');
    if isempty(s_y2), s_y2 = def_y2; end

    fprintf('   Обчислення похідних для Г2...\n');
    [geom.x2, geom.y2, geom.dx2, geom.dy2, geom.ddx2, geom.ddy2] = ...
        create_curve_with_derivatives(s_x2, s_y2);

    % --- 2. ТЕСТОВІ ТОЧКИ ---
    fprintf('\n>>> Введіть 3 тестові точки всередині області <<<\n');
    test_points = zeros(3, 2);
    def_pts = [2, 0; 1, 1; 0, -2];
    use_def = input('Використати стандартні точки (2,0), (1,1), (0,-2)? y/n [y]: ', 's');
    if isempty(use_def) || lower(use_def) == 'y'
        test_points = def_pts;
    else
        for k=1:3
            s_pt = input(sprintf('Точка %d [x, y]: ', k), 's');
            test_points(k,:) = eval(['[', s_pt, ']']);
        end
    end

    % --- 3. ТОЧНИЙ РОЗВ'ЯЗОК (для генерації RHS) ---
    % Фундаментальний розв'язок: Phi(x,y) = -1/(2*pi) * ln(|x-y|)
    % Для задачі Лапласа в 2D
    y_star = [5, 5]; % Джерело поза областю
    u_exact = @(x) -1/(2*pi) * log(sqrt(sum((x - y_star).^2)));
    
    % Нормальна похідна: du/dn = grad(u) * n
    % grad(u) = -1/(2*pi) * (x-y*) / |x-y*|^2
    du_dn_exact = @(x, nu) -1/(2*pi) * dot(x - y_star, nu) / sum((x - y_star).^2);

    % --- 4. ОБЧИСЛЕННЯ ---
    N_values = [8, 16, 32, 64];
    
    % Структура для збереження всіх результатів
    results = struct();
    results.N_values = N_values;
    results.test_points = test_points;
    results.errors = cell(length(N_values), 1);
    results.solutions = cell(length(N_values), 1);

    fprintf('\n========================================================================\n');
    fprintf(' ТАБЛИЦЯ ПОХИБОК (Абсолютна похибка |u_approx - u_exact|)\n');
    fprintf('------------------------------------------------------------------------\n');
    fprintf(' %-4s | x=(%g,%g) | x=(%g,%g) | x=(%g,%g) | Відносна похибка\n', ...
        'N', test_points(1,:), test_points(2,:), test_points(3,:));
    fprintf('------------------------------------------------------------------------\n');

    for idx = 1:length(N_values)
        N = N_values(idx);
        h = 2*pi / N;
        t = (0 : N-1)' * h;

        % Отримуємо геометрію. is_hole=1 для Г1 (переверне нормаль)
        [geom1, singularity_flag1] = get_geometry(t, geom.x1, geom.y1, ...
            geom.dx1, geom.dy1, geom.ddx1, geom.ddy1, 1);
        [geom2, singularity_flag2] = get_geometry(t, geom.x2, geom.y2, ...
            geom.dx2, geom.dy2, geom.ddx2, geom.ddy2, 0);

        % Обробка сингулярностей
        if singularity_flag1 || singularity_flag2
            fprintf('\n⚠️  ВИЯВЛЕНО СИНГУЛЯРНІСТЬ!\n');
            if singularity_flag1
                fprintf('   Проблема на Г1: нульова похідна або дуже малий jacobian\n');
            end
            if singularity_flag2
                fprintf('   Проблема на Г2: нульова похідна або дуже малий jacobian\n');
            end
            cont = input('   Продовжити обчислення? (y/n) [y]: ', 's');
            if ~isempty(cont) && lower(cont) == 'n'
                fprintf('Обчислення перервано користувачем.\n');
                return;
            end
            fprintf('   Продовжуємо з обережністю...\n');
        end

        % Створюємо анонімні функції для ядер
        kernel_funcs = create_kernel_functions();

        Dim = 2*N;
        A = zeros(Dim, Dim);
        RHS = zeros(Dim, 1);
        w = 2*pi / N; % Вага квадратури (метод Нистрьома: рівномірна сітка)

        % Побудова системи
        for i = 1:N
            for j = 1:N
                % L11: Г1 -> Г1 (Подвійний шар)
                A(i, j) = kernel_funcs.L11(i, j, geom1, geom1) * w;
                if i==j
                    A(i,j) = A(i,j) - 0.5; % Стрибок подвійного шару
                end

                % L12: Г2 -> Г1 (Простий шар)
                A(i, j+N) = kernel_funcs.L12(i, j, geom1, geom2) * w;

                % L21: Г1 -> Г2 (Мішана похідна - гіперсингулярне ядро)
                A(i+N, j) = kernel_funcs.L21(i, j, geom2, geom1) * w;

                % L22: Г2 -> Г2 (Похідна простого шару)
                A(i+N, j+N) = kernel_funcs.L22(i, j, geom2, geom2) * w;
                if i==j
                    A(i+N, j+N) = A(i+N, j+N) + 0.5; % Стрибок похідної простого шару
                end
            end

            % Праві частини
            RHS(i) = u_exact(geom1.P(i,:));
            RHS(i+N) = du_dn_exact(geom2.P(i,:), geom2.nu(i,:));
        end

        % Розв'язок системи
        try
            Psi = A \ RHS;
        catch ME
            fprintf('\n❌ ПОМИЛКА ПРИ РОЗВ''ЯЗУВАННІ СИСТЕМИ:\n');
            fprintf('   %s\n', ME.message);
            fprintf('   Умовність матриці: cond(A) = %.2e\n', cond(A));
            cont = input('   Спробувати з іншим N? (y/n) [n]: ', 's');
            if isempty(cont) || lower(cont) ~= 'y'
                return;
            end
            continue;
        end
        
        psi1 = Psi(1:N);
        psi2 = Psi(N+1:end);

        % Збереження повних результатів
        sol = struct();
        sol.N = N;
        sol.psi1 = psi1;
        sol.psi2 = psi2;
        sol.geom1 = geom1;
        sol.geom2 = geom2;
        sol.A = A;
        sol.RHS = RHS;
        results.solutions{idx} = sol;

        % Обчислення похибок
        errs = zeros(1, 3);
        u_exact_vals = zeros(1, 3);
        for k=1:3
            u_val = recover_solution(test_points(k,:), psi1, psi2, geom1, geom2);
            u_exact_vals(k) = u_exact(test_points(k,:));
            errs(k) = abs(u_val - u_exact_vals(k));
        end
        
        % Відносна похибка
        rel_err = mean(errs ./ (abs(u_exact_vals) + 1e-10));
        
        results.errors{idx} = struct('absolute', errs, 'relative', rel_err, ...
            'u_exact', u_exact_vals);
        
        fprintf(' %-4d | %.4e   | %.4e   | %.4e   | %.2e\n', ...
            N, errs(1), errs(2), errs(3), rel_err);
    end
    fprintf('------------------------------------------------------------------------\n');

    % Збереження останнього розв'язку для візуалізації
    final_sol = results.solutions{end};
    results.final_solution = final_sol;

    % --- 5. ВІЗУАЛІЗАЦІЯ ---
    fprintf('\nПобудова графіків...\n');
    
    % Графік 1: Геометрія
    figure(1); clf; hold on; axis equal; grid on; box on;
    title('Геометрія області та тестові точки', 'FontSize', 12);
    tt = linspace(0, 2*pi, 200)';
    plot(geom.x1(tt), geom.y1(tt), 'b-', 'LineWidth', 2);
    plot(geom.x2(tt), geom.y2(tt), 'r-', 'LineWidth', 2);
    plot(test_points(:,1), test_points(:,2), 'ko', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
    legend('\Gamma_1 (Діріхле)', '\Gamma_2 (Нейман)', 'Test Points', 'Location', 'best');
    xlabel('x'); ylabel('y');

    % Графік 2: Наближений розв'язок та лінії рівня
    figure(2); clf; hold on; axis equal;
    title('Наближений розв''язок u(x) та лінії рівня', 'FontSize', 12);
    
    % Створюємо сітку
    x_min = min(final_sol.geom2.P(:,1)) - 0.5;
    x_max = max(final_sol.geom2.P(:,1)) + 0.5;
    y_min = min(final_sol.geom2.P(:,2)) - 0.5;
    y_max = max(final_sol.geom2.P(:,2)) + 0.5;
    
    nx = 80; ny = 80;
    [X, Y] = meshgrid(linspace(x_min, x_max, nx), linspace(y_min, y_max, ny));
    U_surf = nan(size(X));

    fprintf('   Обчислення розв''язку на сітці...\n');
    for k = 1:numel(X)
        pt = [X(k), Y(k)];
        % Перевірка на "всередині Г2 і ззовні Г1"
        in_outer = inpolygon(pt(1), pt(2), final_sol.geom2.P(:,1), final_sol.geom2.P(:,2));
        in_inner = inpolygon(pt(1), pt(2), final_sol.geom1.P(:,1), final_sol.geom1.P(:,2));

        if in_outer && ~in_inner
            try
                U_surf(k) = recover_solution(pt, final_sol.psi1, final_sol.psi2, ...
                    final_sol.geom1, final_sol.geom2);
            catch
                U_surf(k) = NaN;
            end
        end
    end
    
    % Контурний графік з лініями рівня
    contourf(X, Y, U_surf, 25);
    colorbar;
    c = colorbar;
    c.Label.String = 'u(x,y)';
    
    % Додаємо контури границь
    plot(geom.x1(tt), geom.y1(tt), 'w-', 'LineWidth', 2);
    plot(geom.x2(tt), geom.y2(tt), 'k-', 'LineWidth', 2);
    
    % Додаємо тестові точки
    plot(test_points(:,1), test_points(:,2), 'wo', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
    
    xlabel('x'); ylabel('y');
    
    % Графік 3: Збіжність
    figure(3); clf; hold on; grid on; box on;
    title('Збіжність методу', 'FontSize', 12);
    abs_errors = cell2mat(cellfun(@(e) e.absolute, results.errors, 'UniformOutput', false));
    semilogy(N_values, abs_errors(:,1), 'o-', 'LineWidth', 2, 'MarkerSize', 8);
    semilogy(N_values, abs_errors(:,2), 's-', 'LineWidth', 2, 'MarkerSize', 8);
    semilogy(N_values, abs_errors(:,3), '^-', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('N (кількість вузлів)');
    ylabel('Абсолютна похибка');
    legend(sprintf('Точка 1: (%.1f, %.1f)', test_points(1,:)), ...
           sprintf('Точка 2: (%.1f, %.1f)', test_points(2,:)), ...
           sprintf('Точка 3: (%.1f, %.1f)', test_points(3,:)), ...
           'Location', 'best');
    
    fprintf('\n✅ Обчислення завершено!\n');
    fprintf('   Результати збережено в структурі results.\n');
end

% =========================================================================
% ДОПОМІЖНІ ФУНКЦІЇ
% =========================================================================

function [fx, fy, fdx, fdy, fddx, fddy] = create_curve_with_derivatives(s_x, s_y)
    % Створює анонімні функції для кривої та її похідних
    % Використовує символьне обчислення, якщо доступне, інакше чисельне
    
    % Основні функції
    fx = str2func(['@(t) ' s_x]);
    fy = str2func(['@(t) ' s_y]);
    
    % Спробуємо символьне обчислення
    try
        if exist('syms', 'file')
            syms t_sym;
            x_sym = eval(s_x);
            y_sym = eval(s_y);
            dx_sym = diff(x_sym, t_sym);
            dy_sym = diff(y_sym, t_sym);
            ddx_sym = diff(dx_sym, t_sym);
            ddy_sym = diff(dy_sym, t_sym);
            
            fdx = matlabFunction(dx_sym, 'Vars', t_sym);
            fdy = matlabFunction(dy_sym, 'Vars', t_sym);
            fddx = matlabFunction(ddx_sym, 'Vars', t_sym);
            fddy = matlabFunction(ddy_sym, 'Vars', t_sym);
            return;
        end
    catch
        % Якщо символьне обчислення не вдалося, використовуємо чисельне
    end
    
    % Чисельне обчислення похідних
    h_num = 1e-6;
    fdx = @(t) (fx(t + h_num) - fx(t - h_num)) / (2 * h_num);
    fdy = @(t) (fy(t + h_num) - fy(t - h_num)) / (2 * h_num);
    
    % Другі похідні
    fddx = @(t) (fx(t + h_num) - 2*fx(t) + fx(t - h_num)) / (h_num^2);
    fddy = @(t) (fy(t + h_num) - 2*fy(t) + fy(t - h_num)) / (h_num^2);
end

function [geom, singularity_flag] = get_geometry(t, fx, fy, fdx, fdy, fddx, fddy, is_hole)
    % Обчислює геометрію з перевіркою сингулярностей
    % is_hole = 1 -> Нормаль розвертається всередину отвору
    
    x = fx(t);
    y = fy(t);
    dx = fdx(t);
    dy = fdy(t);
    ddx = fddx(t);
    ddy = fddy(t);

    P = [x(:), y(:)];
    dP = [dx(:), dy(:)];
    ddP = [ddx(:), ddy(:)];
    jac = sqrt(dx.^2 + dy.^2);

    % Перевірка сингулярностей
    singularity_flag = false;
    if any(jac < 1e-10)
        singularity_flag = true;
        warning('Виявлено нульовий або дуже малий jacobian!');
    end

    % Нормаль: (dy, -dx) / |dP|
    nu = [dy(:), -dx(:)] ./ jac(:);

    if is_hole
        nu = -nu; % Для дірки зовнішня нормаль до D дивиться в центр дірки
    end

    geom = struct();
    geom.P = P;
    geom.dP = dP;
    geom.ddP = ddP;
    geom.nu = nu;
    geom.jac = jac(:);
end

function kernels = create_kernel_functions()
    % Створює структуру анонімних функцій для ядер
    
    % Фундаментальний розв'язок: Phi(x,y) = -1/(2*pi) * ln(|x-y|)
    % Для 2D задачі Лапласа
    
    kernels = struct();
    
    % L11: Подвійний шар на Г1
    % K(x,y) = dPhi/dn_y = (x-y)*n_y / |x-y|^2 * (1/(2*pi))
    kernels.L11 = @(i, j, geom_obs, geom_src) ...
        compute_L11(i, j, geom_obs, geom_src);
    
    % L12: Простий шар Г2->Г1
    % K(x,y) = Phi(x,y) = -1/(2*pi) * ln(|x-y|)
    kernels.L12 = @(i, j, geom_obs, geom_src) ...
        compute_L12(i, j, geom_obs, geom_src);
    
    % L21: Мішана похідна Г1->Г2 (гіперсингулярне)
    % K(x,y) = d^2Phi/(dn_x dn_y) = -1/(2*pi) * (n_x*n_y - 2*(r*n_x)*(r*n_y)/r^2) / r^2
    kernels.L21 = @(i, j, geom_obs, geom_src) ...
        compute_L21(i, j, geom_obs, geom_src);
    
    % L22: Похідна простого шару Г2->Г2
    % K(x,y) = dPhi/dn_x = (y-x)*n_x / |x-y|^2 * (1/(2*pi))
    kernels.L22 = @(i, j, geom_obs, geom_src) ...
        compute_L22(i, j, geom_obs, geom_src);
end

function val = compute_L11(i, j, geom_obs, geom_src)
    % L11: Подвійний шар
    % K(x,y) = (1/(2*pi)) * (x-y)*n_y / |x-y|^2
    if i == j
        % Границя через кривину: kappa / (2*|dP|)
        % де kappa = (dP'' * n) / |dP|
        curv = dot(geom_src.ddP(i,:), geom_src.nu(i,:));
        val = (1/(2*pi)) * curv / (2 * geom_src.jac(i));
    else
        r = geom_obs.P(i,:) - geom_src.P(j,:);
        d2 = sum(r.^2);
        if d2 < 1e-14
            val = 0;
        else
            val = (1/(2*pi)) * dot(r, geom_src.nu(j,:)) / d2 * geom_src.jac(j);
        end
    end
end

function val = compute_L12(i, j, geom_obs, geom_src)
    % L12: Простий шар
    % K(x,y) = -1/(2*pi) * ln(|x-y|)
    r = geom_obs.P(i,:) - geom_src.P(j,:);
    dist = sqrt(sum(r.^2));
    if dist < 1e-14
        val = 0; % Обробка сингулярності
    else
        val = -1/(2*pi) * log(dist) * geom_src.jac(j);
    end
end

function val = compute_L21(i, j, geom_obs, geom_src)
    % L21: Мішана похідна (гіперсингулярне ядро)
    % K(x,y) = -1/(2*pi) * (n_x*n_y - 2*(r*n_x)*(r*n_y)/r^2) / r^2
    r = geom_obs.P(i,:) - geom_src.P(j,:);
    d2 = sum(r.^2);
    
    if d2 < 1e-14
        % Для сингулярного випадку потрібна спеціальна обробка
        % Тут повертаємо 0, але в реальності потрібна regularization
        val = 0;
    else
        term1 = dot(geom_obs.nu(i,:), geom_src.nu(j,:));
        term2 = 2 * dot(r, geom_obs.nu(i,:)) * dot(r, geom_src.nu(j,:)) / d2;
        val = -1/(2*pi) * (term1 - term2) / d2 * geom_src.jac(j);
    end
end

function val = compute_L22(i, j, geom_obs, geom_src)
    % L22: Похідна простого шару
    % K(x,y) = (1/(2*pi)) * (y-x)*n_x / |x-y|^2
    if i == j
        % Границя через кривину (аналогічно L11)
        curv = dot(geom_obs.ddP(i,:), geom_obs.nu(i,:));
        val = (1/(2*pi)) * curv / (2 * geom_obs.jac(i));
    else
        r = geom_src.P(j,:) - geom_obs.P(i,:); % y - x
        d2 = sum(r.^2);
        if d2 < 1e-14
            val = 0;
        else
            val = (1/(2*pi)) * dot(r, geom_obs.nu(i,:)) / d2 * geom_src.jac(j);
        end
    end
end

function u = recover_solution(pt, psi1, psi2, geom1, geom2)
    % Відновлення розв'язку u(x) за формулою представлення
    % u(x) = int_Г1 K1(x,y)*psi1(y) ds(y) + int_Г2 K2(x,y)*psi2(y) ds(y)
    
    N = length(psi1);
    w = 2*pi / N;
    
    % Інтеграл по Г1 (Подвійний шар): K = (1/(2*pi)) * (x-y)*n_y / |x-y|^2
    R1 = pt - geom1.P;
    d1_sq = sum(R1.^2, 2);
    % Уникаємо ділення на нуль
    d1_sq(d1_sq < 1e-14) = 1e-14;
    term1 = sum(R1 .* geom1.nu, 2) ./ d1_sq;
    int1 = (1/(2*pi)) * sum(psi1 .* term1 .* geom1.jac) * w;
    
    % Інтеграл по Г2 (Простий шар): K = -1/(2*pi) * ln(|x-y|)
    R2 = pt - geom2.P;
    dist2 = sqrt(sum(R2.^2, 2));
    dist2(dist2 < 1e-14) = 1e-14;
    term2 = -1/(2*pi) * log(dist2);
    int2 = sum(psi2 .* term2 .* geom2.jac) * w;
    
    u = int1 + int2;
end
