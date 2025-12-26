function generate_complete_fitting_results(density_params_init, density_params_opt, exp_data_multi, grid_model, input_params)
% GENERATE_COMPLETE_FITTING_RESULTS 生成完整的拟合结果报告
% 统一格式：MF_R评估、台阶显示、χ/λ²、完整兼容输出、基底约束信息

    %% 创建结果文件夹
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    result_folder = fullfile('fittingresult', timestamp);
    if ~exist('fittingresult', 'dir')
        mkdir('fittingresult');
    end
    mkdir(result_folder);
    fprintf('结果保存到: %s\n', result_folder);
    
    %% 构建初始和优化后的模型
    grid_model_init = update_density_distribution(grid_model, input_params, density_params_init);
    grid_model_opt = update_density_distribution(grid_model, input_params, density_params_opt);
    
    %% 计算MF_R (纯反射率Merit Function，不加权重)
    fprintf('计算MF_R评估指标...\n');
    MF_R_individual = zeros(length(exp_data_multi), 1);
    sum_squared_residuals = 0;
    total_points = 0;
    
    for i = 1:length(exp_data_multi)
        theta_exp = exp_data_multi(i).theta;
        R_exp = exp_data_multi(i).R;
        energy_index = exp_data_multi(i).energy_index;
        
        R_calc = calculate_reflectivity(grid_model_opt, input_params, energy_index, theta_exp);
        residuals = (R_exp - R_calc) ./ R_exp;
        MF_R_individual(i) = sum(residuals.^2);
        
        sum_squared_residuals = sum_squared_residuals + sum(residuals.^2);
        total_points = total_points + length(theta_exp);
    end
    MF_R_total = sum_squared_residuals / total_points;
    
    fprintf('总体MF_R: %.6f\n', MF_R_total);
    
    %% 检查基底约束效果 ← 新增
    if input_params.substrate_constraint_weight > 0
        fprintf('检查基底约束效果...\n');
        [substrate_analysis] = analyze_substrate_constraint_effect(density_params_opt, grid_model, input_params);
    end
    
    %% 第一部分：反射率拟合对比图（包含初始vs拟合）
    fprintf('生成反射率拟合图...\n');
    fig1 = figure('Name', '多波长反射率拟合结果', 'Position', [100, 100, 1200, 800]);
    
    % 存储反射率数据用于导出
    reflectivity_data = struct();
    
    for i = 1:length(exp_data_multi)
        subplot(2, 3, i);
        
        theta_exp = exp_data_multi(i).theta;
        R_exp = exp_data_multi(i).R;
        energy_index = exp_data_multi(i).energy_index;
        wavelength = input_params.wavelengths(energy_index);
        energy = input_params.energies(energy_index);
        
        % 计算连续角度的反射率用于绘图
        theta_max = max(theta_exp) * 1.1;
        theta_calc = linspace(0, theta_max, 300);
        R_init_calc = calculate_reflectivity(grid_model_init, input_params, energy_index, theta_calc);
        R_opt_calc = calculate_reflectivity(grid_model_opt, input_params, energy_index, theta_calc);
        R_calc_exp = calculate_reflectivity(grid_model_opt, input_params, energy_index, theta_exp);
        
        % 绘图（添加初始计算结果）
        semilogy(theta_calc, R_init_calc, 'r--', 'LineWidth', 1.5, 'DisplayName', '初始计算');
        hold on;
        semilogy(theta_calc, R_opt_calc, 'b-', 'LineWidth', 2, 'DisplayName', '拟合结果');
        semilogy(theta_exp, R_exp, 'ko', 'MarkerSize', 4, 'DisplayName', '实验数据');
        
        xlabel('角度 (度)');
        ylabel('反射率');
        title(sprintf('%.1f eV (%.3f nm)\nMF_R = %.4f', energy, wavelength, MF_R_individual(i)));
        legend('Location', 'best', 'FontSize', 8);
        grid on;
        xlim([0, theta_max]);
        hold off;
        
        % 存储数据
        field_name = sprintf('wl_%d_nm', round(wavelength*1000));
        reflectivity_data.(field_name).theta_exp = theta_exp;
        reflectivity_data.(field_name).R_exp = R_exp;
        reflectivity_data.(field_name).R_calc = R_calc_exp;
        reflectivity_data.(field_name).theta_smooth = theta_calc;
        reflectivity_data.(field_name).R_smooth = R_opt_calc;
        reflectivity_data.(field_name).MF_R = MF_R_individual(i);
        reflectivity_data.(field_name).energy = energy;
        reflectivity_data.(field_name).wavelength = wavelength;
    end
    
    sgtitle(sprintf('多波长反射率拟合结果 (总体MF_R = %.6f)', MF_R_total), ...
            'FontSize', 14, 'FontWeight', 'bold');
    
    % 保存图片和数据
    savefig(fig1, fullfile(result_folder, 'reflectivity_fitting.fig'));
    
    % 导出反射率数据
    fid = fopen(fullfile(result_folder, 'reflectivity_data.txt'), 'w');
    fprintf(fid, '%% 多波长反射率拟合数据\n');
    fprintf(fid, '%% 生成时间: %s\n', datestr(now));
    fprintf(fid, '%% 总体MF_R: %.6f\n', MF_R_total);
    fprintf(fid, '%% MF_R计算: sum((R_exp-R_calc)/R_exp)^2，不加权重\n\n');
    
    for i = 1:length(exp_data_multi)
        wavelength = input_params.wavelengths(exp_data_multi(i).energy_index);
        energy = input_params.energies(exp_data_multi(i).energy_index);
        fprintf(fid, '%% 波长: %.3f nm, 能量: %.1f eV, MF_R: %.6f\n', ...
                wavelength, energy, MF_R_individual(i));
        fprintf(fid, 'theta_exp\tR_exp\tR_calc\n');
        
        field_name = sprintf('wl_%d_nm', round(wavelength*1000));
        data = reflectivity_data.(field_name);
        for j = 1:length(data.theta_exp)
            fprintf(fid, '%.6f\t%.6e\t%.6e\n', ...
                    data.theta_exp(j), data.R_exp(j), data.R_calc(j));
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
    close(fig1);
    
    %% 第二部分：密度分布对比图（初始vs拟合，台阶显示）
    fprintf('生成密度分布对比图...\n');
    
    z_grid = grid_model.z_grid;
    density_grid_init = grid_model_init.density_grid;
    density_grid_opt = grid_model_opt.density_grid;
    
    % 检测实际存在的拟合元素
    existing_elements = {};
    for i = 1:length(input_params.fit_elements)
        elem = input_params.fit_elements{i};
        if max(density_grid_opt.(elem)) > 1e-6
            existing_elements{end+1} = elem;
        end
    end
    
    % 检测非拟合但有密度的元素
    all_elements = {'Si', 'Co', 'Mg', 'O', 'C'};
    other_elements = setdiff(all_elements, input_params.fit_elements);
    for i = 1:length(other_elements)
        elem = other_elements{i};
        if isfield(density_grid_opt, elem) && max(density_grid_opt.(elem)) > 1e-6
            existing_elements{end+1} = elem;
        end
    end
    
    if isempty(existing_elements)
        fprintf('警告：未检测到任何有效的元素密度分布\n');
        return;
    end
    
    % 创建密度分布对比图
    fig2 = figure('Name', '密度分布对比 (初始vs拟合)', 'Position', [200, 200, 800, 600]);
    
    % 颜色映射
    element_colors = containers.Map({'Si', 'Co', 'Mg', 'O', 'C'}, ...
                                   {[1,0,0], [0,0,1], [0,1,0], [1,0.5,0], [0.5,0,1]});
    
    for i = 1:length(existing_elements)
        elem = existing_elements{i};
        density_init = density_grid_init.(elem);
        density_opt = density_grid_opt.(elem);
        
        if isKey(element_colors, elem)
            color = element_colors(elem);
        else
            color = rand(1,3);
        end
        
        % 使用stairs显示台阶（体现dz分层特性）
        if ismember(elem, input_params.fit_elements)
            % 拟合元素：显示初始和拟合
            stairs(z_grid, density_init, '--', 'Color', color, 'LineWidth', 1.5, ...
                   'DisplayName', sprintf('%s (初始)', elem));
            hold on;
            stairs(z_grid, density_opt, '-', 'Color', color, 'LineWidth', 2.5, ...
                   'DisplayName', sprintf('%s (拟合)', elem));
        else
            % 非拟合元素：只显示固定值
            stairs(z_grid, density_opt, ':', 'Color', color, 'LineWidth', 1, ...
                   'DisplayName', sprintf('%s (固定)', elem));
            hold on;
        end
    end
    
    % 标记基底约束区域 ← 新增
    if input_params.substrate_constraint_weight > 0 && ismember('Si', existing_elements)
        num_constrained = input_params.num_substrate_constraint_layers;
        constraint_start_idx = length(z_grid) - num_constrained + 1;
        if constraint_start_idx >= 1
            constraint_z_start = z_grid(constraint_start_idx);
            constraint_z_end = max(z_grid);
            y_lim = ylim;
            fill([constraint_z_start, constraint_z_end, constraint_z_end, constraint_z_start], ...
                 [y_lim(1), y_lim(1), y_lim(2), y_lim(2)], ...
                 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none', ...
                 'DisplayName', sprintf('基底约束区(最后%d层)', num_constrained));
        end
    end
    
    xlabel('深度 (nm)');
    ylabel('密度 (g/cm³)');
    title('密度分布对比 (初始vs拟合)');
    legend('Location', 'best');
    grid on;
    xlim([0, max(z_grid)]);
    ylim([0, max(ylim)*1.05]);
    hold off;
    
    % 保存图片
    savefig(fig2, fullfile(result_folder, 'density_distribution.fig'));
    
    % 导出密度数据
    fid = fopen(fullfile(result_folder, 'density_distribution.txt'), 'w');
    fprintf(fid, '%% 密度分布对比数据 (初始vs拟合)\n');
    fprintf(fid, '%% 生成时间: %s\n', datestr(now));
    fprintf(fid, '%% 网格步长: %.3f nm\n', grid_model.dz);
    
    % 添加基底约束信息 ← 新增
    if input_params.substrate_constraint_weight > 0
        fprintf(fid, '%% 基底约束: 最后%d层Si密度 → %.2f g/cm³ (权重=%.1f)\n', ...
                input_params.num_substrate_constraint_layers, ...
                input_params.substrate_density, ...
                input_params.substrate_constraint_weight);
    end
    
    % 写入表头
    fprintf(fid, 'depth(nm)');
    density_fields = {};
    for i = 1:length(existing_elements)
        elem = existing_elements{i};
        if ismember(elem, input_params.fit_elements)
            density_fields{end+1} = sprintf('%s_init', elem);
            density_fields{end+1} = sprintf('%s_fit', elem);
            fprintf(fid, '\t%s_init(g/cm3)\t%s_fit(g/cm3)', elem, elem);
        else
            density_fields{end+1} = sprintf('%s_fixed', elem);
            fprintf(fid, '\t%s_fixed(g/cm3)', elem);
        end
    end
    fprintf(fid, '\n');
    
    % 写入数据
    for i = 1:length(z_grid)
        fprintf(fid, '%.6f', z_grid(i));
        for j = 1:length(existing_elements)
            elem = existing_elements{j};
            if ismember(elem, input_params.fit_elements)
                fprintf(fid, '\t%.6f\t%.6f', density_grid_init.(elem)(i), density_grid_opt.(elem)(i));
            else
                fprintf(fid, '\t%.6f', density_grid_opt.(elem)(i));
            end
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
    close(fig2);
    
    %% 保存标准density_distribution.mat文件
    fprintf('保存density_distribution.mat文件...\n');
    
    % 构建与层参数项目兼容的数据结构
    density_data = struct();
    density_data.z_grid = z_grid;
    density_data.timestamp = timestamp;
    density_data.grid_step = grid_model.dz;
    density_data.calculation_depth = max(z_grid);
    density_data.source = 'partial_density_fitting';
    
    % 基底约束信息 ← 新增
    if input_params.substrate_constraint_weight > 0
        density_data.substrate_constraint.enabled = true;
        density_data.substrate_constraint.layers = input_params.num_substrate_constraint_layers;
        density_data.substrate_constraint.weight = input_params.substrate_constraint_weight;
        density_data.substrate_constraint.target_density = input_params.substrate_density;
    else
        density_data.substrate_constraint.enabled = false;
    end
    
    % 按照标准字段名保存密度数据
    standard_elements = {'Si', 'Co', 'Mg', 'O', 'C'};
    fprintf('\n元素密度统计:\n');
    
    for i = 1:length(standard_elements)
        elem = standard_elements{i};
        field_name = sprintf('density_%s', elem);
        if isfield(density_grid_opt, elem)
            density_data.(field_name) = density_grid_opt.(elem);
        else
            density_data.(field_name) = zeros(size(z_grid));
        end
        
        % 统计信息
        max_density = max(density_data.(field_name));
        nonzero_points = sum(density_data.(field_name) > 1e-6);
        fprintf('  %s: 最大密度=%.4f g/cm³, 非零点=%d\n', elem, max_density, nonzero_points);
    end
    
    % 保存MAT文件
    save(fullfile(result_folder, 'density_distribution.mat'), 'density_data', '-v7.3');
    fprintf('✅ density_distribution.mat 已保存\n');
    
    %% 第三部分：极化率/λ²分布图（台阶显示）
    fprintf('生成极化率/λ²分布图...\n');
    fig3 = figure('Name', '极化率/λ²分布 (台阶)', 'Position', [300, 300, 1000, 600]);
    
    % 为每个波长计算并显示极化率/λ²
    colors = lines(length(input_params.wavelengths));
    chi_data = struct();
    
    for i = 1:length(input_params.wavelengths)
        wavelength = input_params.wavelengths(i);
        energy = input_params.energies(i);
        
        % 计算极化率
        [chi_real, chi_imag] = density_to_chi(grid_model_opt, input_params, i);
        
        % 计算χ/λ²（消除波长依赖性）
        chi_real_over_lambda2 = chi_real / (wavelength^2);
        chi_imag_over_lambda2 = chi_imag / (wavelength^2);
        
        % 存储数据
        field_name = sprintf('wl_%d_nm', round(wavelength*1000));
        chi_data.(field_name).chi_real_over_lambda2 = chi_real_over_lambda2;
        chi_data.(field_name).chi_imag_over_lambda2 = chi_imag_over_lambda2;
        chi_data.(field_name).energy = energy;
        chi_data.(field_name).wavelength = wavelength;
        
        % 绘制实部/λ²（使用stairs显示台阶）
        subplot(1, 2, 1);
        stairs(z_grid, chi_real_over_lambda2, 'Color', colors(i,:), 'LineWidth', 2, ...
               'DisplayName', sprintf('%.1f eV', energy));
        hold on;
        
        % 绘制虚部/λ²（取负号，使用stairs显示台阶）
        subplot(1, 2, 2);
        stairs(z_grid, -chi_imag_over_lambda2, 'Color', colors(i,:), 'LineWidth', 2, ...
               'DisplayName', sprintf('%.1f eV', energy));
        hold on;
    end
    
    subplot(1, 2, 1);
    xlabel('深度 (nm)');
    ylabel('Re(χ/λ²) (nm⁻²)');
    title('极化率实部/λ²');
    legend('Location', 'best');
    grid on;
    
    subplot(1, 2, 2);
    xlabel('深度 (nm)');
    ylabel('-Im(χ/λ²) (nm⁻²)');
    title('极化率虚部/λ²');
    legend('Location', 'best');
    grid on;
    
    sgtitle('极化率/λ² 分布 (台阶显示)', 'FontSize', 14, 'FontWeight', 'bold');
    
    % 保存图片
    savefig(fig3, fullfile(result_folder, 'chi_over_lambda2.fig'));
    
    % 导出极化率/λ²数据
    fid = fopen(fullfile(result_folder, 'chi_over_lambda2.txt'), 'w');
    fprintf(fid, '%% 极化率/λ²分布数据 (台阶显示)\n');
    fprintf(fid, '%% 生成时间: %s\n', datestr(now));
    fprintf(fid, '%% 说明: 极化率除以波长平方，消除波长依赖性，体现dz分层特性\n\n');
    
    % 写入表头
    fprintf(fid, 'depth(nm)');
    for i = 1:length(input_params.wavelengths)
        energy = input_params.energies(i);
        fprintf(fid, '\tRe_chi/λ²_%.1feV\t-Im_chi/λ²_%.1feV', energy, energy);
    end
    fprintf(fid, '\n');
    
    % 写入数据
    for j = 1:length(z_grid)
        fprintf(fid, '%.6f', z_grid(j));
        for i = 1:length(input_params.wavelengths)
            wavelength = input_params.wavelengths(i);
            field_name = sprintf('wl_%d_nm', round(wavelength*1000));
            fprintf(fid, '\t%.6e\t%.6e', ...
                    chi_data.(field_name).chi_real_over_lambda2(j), ...
                    -chi_data.(field_name).chi_imag_over_lambda2(j));
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
    close(fig3);
    
    %% 生成总结报告
    fid = fopen(fullfile(result_folder, 'fitting_summary.txt'), 'w');
    fprintf(fid, '多波长X射线反射率部分密度拟合结果总结\n');
    fprintf(fid, '=====================================\n');
    fprintf(fid, '生成时间: %s\n', datestr(now));
    fprintf(fid, '结果文件夹: %s\n\n', result_folder);
    
    fprintf(fid, '拟合质量评估 (MF_R - 纯反射率Merit Function):\n');
    fprintf(fid, '计算公式: MF_R = sum((R_exp - R_calc)/R_exp)^2 (不加权重)\n');
    fprintf(fid, '总体MF_R: %.6f\n\n', MF_R_total);
    
    fprintf(fid, '各波长MF_R详情:\n');
    for i = 1:length(exp_data_multi)
        wavelength = input_params.wavelengths(exp_data_multi(i).energy_index);
        energy = input_params.energies(exp_data_multi(i).energy_index);
        fprintf(fid, '%.3f nm (%.1f eV): MF_R = %.6f\n', ...
                wavelength, energy, MF_R_individual(i));
    end
    
    fprintf(fid, '\n拟合元素: %s\n', strjoin(input_params.fit_elements, ', '));
    fprintf(fid, '检测到的元素: %s\n', strjoin(existing_elements, ', '));
    
    % 基底约束信息 ← 新增
    fprintf(fid, '\n基底约束设置:\n');
    if input_params.substrate_constraint_weight > 0
        fprintf(fid, '约束状态: 已启用\n');
        fprintf(fid, '约束层数: 最后%d层\n', input_params.num_substrate_constraint_layers);
        fprintf(fid, '目标密度: %.2f g/cm³ (基底Si)\n', input_params.substrate_density);
        fprintf(fid, '约束权重: %.1f\n', input_params.substrate_constraint_weight);
        
        % 添加基底约束效果分析
        if exist('substrate_analysis', 'var')
            fprintf(fid, '约束效果: 最后%d层平均偏差=%.4f g/cm³\n', ...
                    substrate_analysis.num_layers, substrate_analysis.avg_deviation);
        end
    else
        fprintf(fid, '约束状态: 未启用\n');
    end
    
    fprintf(fid, '\n显示特性:\n');
    fprintf(fid, '- 密度分布: 台阶显示，体现dz=%.3f nm分层特性\n', grid_model.dz);
    fprintf(fid, '- 极化率分布: χ/λ²显示，消除波长依赖性\n');
    fprintf(fid, '- 反射率拟合: 包含初始(虚线) vs 拟合(实线)对比\n');
    if input_params.substrate_constraint_weight > 0
        fprintf(fid, '- 基底约束: 红色阴影区域标记约束层\n');
    end
    
    fprintf(fid, '\n生成文件列表:\n');
    fprintf(fid, '- reflectivity_fitting.fig: 反射率拟合对比图 (含初始vs拟合)\n');
    fprintf(fid, '- reflectivity_data.txt: 反射率拟合数据\n');
    fprintf(fid, '- density_distribution.fig: 密度分布对比图 (台阶显示)\n');
    fprintf(fid, '- density_distribution.txt: 密度分布对比数据\n');
    fprintf(fid, '- density_distribution.mat: 密度分布MAT文件（标准兼容格式）\n');
    fprintf(fid, '- chi_over_lambda2.fig: 极化率/λ²分布图 (台阶显示)\n');
    fprintf(fid, '- chi_over_lambda2.txt: 极化率/λ²分布数据\n');
    fprintf(fid, '- fitting_summary.txt: 本总结报告\n');
    
    fclose(fid);
    
    fprintf('\n=== 拟合结果生成完成 ===\n');
    fprintf('总体MF_R: %.6f\n', MF_R_total);
    fprintf('结果已保存到: %s\n', result_folder);
    fprintf('拟合元素: %s\n', strjoin(input_params.fit_elements, ', '));
    
    % 基底约束总结 ← 新增
    if input_params.substrate_constraint_weight > 0
        fprintf('基底约束: 最后%d层Si → %.2f g/cm³ (权重=%.1f)\n', ...
                input_params.num_substrate_constraint_layers, ...
                input_params.substrate_density, ...
                input_params.substrate_constraint_weight);
    else
        fprintf('基底约束: 未启用\n');
    end
    
    fprintf('特点: 台阶显示(dz=%.3f nm)、χ/λ²、初始vs拟合对比、完整兼容输出\n', grid_model.dz);
    
    % 显示各波长MF_R
    fprintf('\n各波长MF_R详情:\n');
    for i = 1:length(exp_data_multi)
        wavelength = input_params.wavelengths(exp_data_multi(i).energy_index);
        energy = input_params.energies(exp_data_multi(i).energy_index);
        fprintf('%.3f nm (%.1f eV): MF_R = %.6f\n', ...
                wavelength, energy, MF_R_individual(i));
    end
end

function substrate_analysis = analyze_substrate_constraint_effect(density_params_opt, grid_model, input_params)
% 分析基底约束的效果
    substrate_analysis = struct();
    
    % 找到Si元素索引
    si_elem_idx = find(strcmp(input_params.fit_elements, 'Si'));
    
    if ~isempty(si_elem_idx)
        num_points = grid_model.num_points;
        start_idx = (si_elem_idx-1) * num_points + 1;
        end_idx = si_elem_idx * num_points;
        si_densities = density_params_opt(start_idx:end_idx);
        
        num_constrained = input_params.num_substrate_constraint_layers;
        target_density = input_params.substrate_density;
        
        % 计算最后几层的平均偏差
        deviations = [];
        for i = 1:num_constrained
            layer_idx = num_points - i + 1;
            if layer_idx >= 1
                deviation = abs(si_densities(layer_idx) - target_density);
                deviations = [deviations; deviation];
            end
        end
        
        substrate_analysis.num_layers = length(deviations);
        substrate_analysis.avg_deviation = mean(deviations);
        substrate_analysis.max_deviation = max(deviations);
        substrate_analysis.target_density = target_density;
        
        fprintf('基底约束效果分析:\n');
        fprintf('  目标密度: %.2f g/cm³\n', target_density);
        fprintf('  约束层数: %d\n', substrate_analysis.num_layers);
        fprintf('  平均偏差: %.4f g/cm³\n', substrate_analysis.avg_deviation);
        fprintf('  最大偏差: %.4f g/cm³\n', substrate_analysis.max_deviation);
        
        % 显示具体的约束层密度
        fprintf('  约束层密度: ');
        for i = 1:min(num_constrained, 3)
            layer_idx = num_points - i + 1;
            if layer_idx >= 1
                fprintf('%.3f ', si_densities(layer_idx));
            end
        end
        fprintf('g/cm³\n');
    else
        substrate_analysis.num_layers = 0;
        substrate_analysis.avg_deviation = 0;
        substrate_analysis.max_deviation = 0;
        substrate_analysis.target_density = input_params.substrate_density;
    end
end