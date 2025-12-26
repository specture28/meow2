function plot_fitting_results(density_params_init, density_params_opt, exp_data_multi, grid_model, input_params)
% PLOT_FITTING_RESULTS 绘制拟合结果的综合图表
%   包含：反射率对比、密度分布、χ分布
%   新增功能：固定元素颜色、χ最大值归一化、导出TXT文件

    fprintf('=== 生成拟合结果图表 ===\n');
    
    % 生成时间戳（用于文件保存）
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    results_dir = 'fitting_results';
    if ~exist(results_dir, 'dir')
        mkdir(results_dir);
    end
    
    % 1. 反射率对比图 (2x3子图布局)
    figure('Name', '反射率拟合对比', 'Position', [100, 100, 1200, 800]);
    reflectivity_data = plot_reflectivity_comparison(density_params_init, density_params_opt, exp_data_multi, grid_model, input_params);
    
    % 2. 密度分布图
    figure('Name', '密度分布', 'Position', [200, 200, 800, 600]);
    density_data = plot_density_profiles(density_params_init, density_params_opt, grid_model, input_params);
    
    % 3. χ分布图 (实部和虚部，最大值归一化)
    figure('Name', 'χ分布', 'Position', [300, 300, 800, 600]);
    chi_data = plot_chi_profiles_normalized(density_params_opt, grid_model, input_params);
    
    % 4. 导出TXT文件
    export_data_to_txt(reflectivity_data, density_data, chi_data, results_dir, timestamp);
    
    fprintf('所有图表生成完成！\n');
end

function reflectivity_data = plot_reflectivity_comparison(density_params_init, density_params_opt, exp_data_multi, grid_model, input_params)
% 绘制反射率对比图 - 6个子图，每个波长一个
    
    % 初始化数据结构
    reflectivity_data = struct();
    
    for i = 1:length(exp_data_multi)
        subplot(2, 3, i);
        
        % 获取实验数据
        theta_exp = exp_data_multi(i).theta;
        R_exp = exp_data_multi(i).R;
        energy_index = exp_data_multi(i).energy_index;
        energy = input_params.energies(energy_index);
        
        % 计算初始反射率
        grid_init = update_density_distribution(grid_model, input_params, density_params_init);
        R_init = calculate_reflectivity(grid_init, input_params, energy_index, theta_exp);
        
        % 计算拟合后反射率
        grid_opt = update_density_distribution(grid_model, input_params, density_params_opt);
        R_fit = calculate_reflectivity(grid_opt, input_params, energy_index, theta_exp);
        
        % 绘制对数坐标图
        semilogy(theta_exp, R_exp, 'ko', 'MarkerSize', 4, 'DisplayName', '实验数据');
        hold on;
        semilogy(theta_exp, R_init, 'r--', 'LineWidth', 1.5, 'DisplayName', '初始计算');
        semilogy(theta_exp, R_fit, 'b-', 'LineWidth', 2, 'DisplayName', '拟合结果');
        
        xlabel('角度 (度)');
        ylabel('反射率');
        title(sprintf('%.1f eV (%.3f nm)', energy, input_params.wavelengths(energy_index)));
        legend('Location', 'best', 'FontSize', 8);
        grid on;
        xlim([min(theta_exp), max(theta_exp)]);
        hold off;
        
        % 保存数据用于导出
        reflectivity_data(i).energy = energy;
        reflectivity_data(i).wavelength = input_params.wavelengths(energy_index);
        reflectivity_data(i).theta = theta_exp;
        reflectivity_data(i).R_exp = R_exp;
        reflectivity_data(i).R_init = R_init;
        reflectivity_data(i).R_fit = R_fit;
    end
    
    sgtitle('反射率拟合对比 - 多波长', 'FontSize', 14, 'FontWeight', 'bold');
end

function density_data = plot_density_profiles(density_params_init, density_params_opt, grid_model, input_params)
% 绘制密度分布图 - 显示拟合前后各元素密度变化（阶梯状）
% 使用固定颜色映射
    
    % 固定元素颜色映射
    element_colors = containers.Map({'Si', 'Co', 'Mg', 'O', 'C'}, ...
                                   {[1,0,0], [0,0,1], [0,1,0], [1,0.5,0], [0.5,0,1]});
    % 红色Si, 蓝色Co, 绿色Mg, 橙色O, 紫色C
    
    % 获取深度网格
    z_grid = grid_model.z_grid;
    
    % 获取初始和优化后的密度分布
    grid_init = update_density_distribution(grid_model, input_params, density_params_init);
    grid_opt = update_density_distribution(grid_model, input_params, density_params_opt);
    
    % 构建阶梯绘图用的深度数组
    z_plot = create_stairs_depth_array(z_grid, grid_model.dz);
    
    % 初始化密度数据结构
    density_data = struct();
    density_data.z_grid = z_grid;
    density_data.z_plot = z_plot;
    
    % 绘制每个拟合元素的密度分布
    for i = 1:length(input_params.fit_elements)
        elem = input_params.fit_elements{i};
        
        density_init = grid_init.density_grid.(elem);
        density_opt = grid_opt.density_grid.(elem);
        
        % 构建阶梯数据
        density_init_plot = create_stairs_data(density_init);
        density_opt_plot = create_stairs_data(density_opt);
        
        % 使用固定颜色
        color = element_colors(elem);
        
        % 初始密度（虚线阶梯）
        stairs(z_plot, density_init_plot, '--', 'Color', color, 'LineWidth', 1.5, ...
               'DisplayName', sprintf('%s (初始)', elem));
        hold on;
        
        % 拟合后密度（实线阶梯）
        stairs(z_plot, density_opt_plot, '-', 'Color', color, 'LineWidth', 2.5, ...
               'DisplayName', sprintf('%s (拟合)', elem));
        
        % 保存数据
        density_data.(sprintf('%s_init', elem)) = density_init;
        density_data.(sprintf('%s_opt', elem)) = density_opt;
    end
    
    % 添加非拟合元素的密度分布（如果有）
    non_fit_elements = setdiff(grid_model.all_elements, input_params.fit_elements);
    if ~isempty(non_fit_elements)
        for i = 1:length(non_fit_elements)
            elem = non_fit_elements{i};
            density_profile = grid_opt.density_grid.(elem);
            if max(density_profile) > 1e-6  % 只显示有意义的密度
                density_plot = create_stairs_data(density_profile);
                if isKey(element_colors, elem)
                    color = element_colors(elem);
                else
                    color = [0.5, 0.5, 0.5];
                end
                stairs(z_plot, density_plot, ':', 'Color', color, 'LineWidth', 1, ...
                       'DisplayName', sprintf('%s (固定)', elem));
                hold on;
                
                % 保存数据
                density_data.(sprintf('%s_fixed', elem)) = density_profile;
            end
        end
    end
    
    xlabel('深度 (nm)');
    ylabel('密度 (g/cm³)');
    title('元素密度分布');
    legend('Location', 'best');
    grid on;
    xlim([0, max(z_plot)]);
    ylim([0, max(ylim)*1.05]);
    hold off;
    
    % 添加基底信息文本
    text(0.02, 0.98, sprintf('基底: %s (%.2f g/cm³)', ...
         grid_model.substrate.material, grid_model.substrate.density), ...
         'Units', 'normalized', 'VerticalAlignment', 'top', ...
         'BackgroundColor', 'white', 'EdgeColor', 'black');
end

function chi_data = plot_chi_profiles_normalized(density_params_opt, grid_model, input_params)
% 绘制最大值归一化的χ分布图 - 实部和虚部随深度变化
% 每个能量分别用自己的最大值归一化，突出起伏模式
    
    % 获取深度网格
    z_grid = grid_model.z_grid;
    
    % 更新密度分布
    grid_opt = update_density_distribution(grid_model, input_params, density_params_opt);
    
    % 构建阶梯绘图用的深度数组
    z_plot = create_stairs_depth_array(z_grid, grid_model.dz);
    
    % 颜色映射
    colors = jet(length(input_params.energies));
    
    % 初始化数据结构
    chi_data = struct();
    chi_data.z_grid = z_grid;
    chi_data.z_plot = z_plot;
    chi_data.energies = input_params.energies;
    
    % 上子图：χ实部（最大值归一化）
    subplot(2, 1, 1);
    for i = 1:length(input_params.energies)
        [chi_real, ~] = density_to_chi(grid_opt, input_params, i);
        
        % 最大值归一化（避免除零）
        chi_real_max = max(abs(chi_real));
        if chi_real_max > 1e-15
            chi_real_normalized = chi_real / chi_real_max;
        else
            chi_real_normalized = chi_real;  % 如果最大值接近零，保持原值
        end
        
        % 构建阶梯数据
        chi_real_plot = create_stairs_data(chi_real_normalized);
        
        stairs(z_plot, chi_real_plot, 'Color', colors(i,:), 'LineWidth', 1.5, ...
               'DisplayName', sprintf('%.1f eV (max=%.2e)', ...
               input_params.energies(i), chi_real_max));
        hold on;
        
        % 保存数据
        chi_data.(sprintf('chi_real_normalized_%d', i)) = chi_real_normalized;
        chi_data.(sprintf('chi_real_max_%d', i)) = chi_real_max;
    end
    xlabel('深度 (nm)');
    ylabel('χ实部 (归一化)');
    title('χ实部分布 (各能量最大值归一化)');
    legend('Location', 'best', 'FontSize', 8);
    grid on;
    xlim([0, max(z_plot)]);
    ylim([-1.1, 1.1]);  % 固定在-1到1范围
    hold off;
    
    % 下子图：χ虚部（最大值归一化）
    subplot(2, 1, 2);
    for i = 1:length(input_params.energies)
        [~, chi_imag] = density_to_chi(grid_opt, input_params, i);
        
        % 最大值归一化（避免除零）
        chi_imag_max = max(abs(chi_imag));
        if chi_imag_max > 1e-15
            chi_imag_normalized = chi_imag / chi_imag_max;
        else
            chi_imag_normalized = chi_imag;  % 如果最大值接近零，保持原值
        end
        
        % 构建阶梯数据
        chi_imag_plot = create_stairs_data(chi_imag_normalized);
        
        stairs(z_plot, chi_imag_plot, 'Color', colors(i,:), 'LineWidth', 1.5, ...
               'DisplayName', sprintf('%.1f eV (max=%.2e)', ...
               input_params.energies(i), chi_imag_max));
        hold on;
        
        % 保存数据
        chi_data.(sprintf('chi_imag_normalized_%d', i)) = chi_imag_normalized;
        chi_data.(sprintf('chi_imag_max_%d', i)) = chi_imag_max;
    end
    xlabel('深度 (nm)');
    ylabel('χ虚部 (归一化)');
    title('χ虚部分布 (各能量最大值归一化)');
    legend('Location', 'best', 'FontSize', 8);
    grid on;
    xlim([0, max(z_plot)]);
    ylim([-1.1, 1.1]);  % 固定在-1到1范围
    hold off;
    
    sgtitle('最大值归一化χ分布 - 多波长起伏对比', 'FontSize', 14, 'FontWeight', 'bold');
    
    fprintf('χ最大值归一化完成\n');
end

function export_data_to_txt(reflectivity_data, density_data, chi_data, results_dir, timestamp)
% 导出三种数据到TXT文件
    
    fprintf('=== 导出数据到TXT文件 ===\n');
    
    %% 1. 导出反射率数据
    ref_filename = fullfile(results_dir, sprintf('reflectivity_data_%s.txt', timestamp));
    fid = fopen(ref_filename, 'w');
    
    % 写入标题
    fprintf(fid, '%% 反射率拟合数据\n');
    fprintf(fid, '%% 时间戳: %s\n', timestamp);
    fprintf(fid, '%%\n');
    
    for i = 1:length(reflectivity_data)
        fprintf(fid, '%% 能量: %.1f eV, 波长: %.3f nm\n', ...
                reflectivity_data(i).energy, reflectivity_data(i).wavelength);
        fprintf(fid, '%% 角度(度)\t实验反射率\t初始计算\t拟合结果\n');
        
        for j = 1:length(reflectivity_data(i).theta)
            fprintf(fid, '%.6f\t%.6e\t%.6e\t%.6e\n', ...
                    reflectivity_data(i).theta(j), ...
                    reflectivity_data(i).R_exp(j), ...
                    reflectivity_data(i).R_init(j), ...
                    reflectivity_data(i).R_fit(j));
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
    fprintf('  反射率数据: %s\n', ref_filename);
    
    %% 2. 导出密度分布数据
    density_filename = fullfile(results_dir, sprintf('density_profiles_%s.txt', timestamp));
    fid = fopen(density_filename, 'w');
    
    % 写入标题
    fprintf(fid, '%% 密度分布数据\n');
    fprintf(fid, '%% 时间戳: %s\n', timestamp);
    fprintf(fid, '%% 深度(nm)');
    
    % 动态写入列标题
    field_names = fieldnames(density_data);
    density_fields = {};
    for i = 1:length(field_names)
        if ~ismember(field_names{i}, {'z_grid', 'z_plot'})
            density_fields{end+1} = field_names{i};
            fprintf(fid, '\t%s(g/cm³)', field_names{i});
        end
    end
    fprintf(fid, '\n');
    
    % 写入数据
    for i = 1:length(density_data.z_grid)
        fprintf(fid, '%.6f', density_data.z_grid(i));
        for j = 1:length(density_fields)
            field = density_fields{j};
            fprintf(fid, '\t%.6f', density_data.(field)(i));
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
    fprintf('  密度分布数据: %s\n', density_filename);
    
    %% 3. 导出χ分布数据
    chi_filename = fullfile(results_dir, sprintf('chi_profiles_%s.txt', timestamp));
    fid = fopen(chi_filename, 'w');
    
    % 写入标题
    fprintf(fid, '%% χ分布数据 (最大值归一化)\n');
    fprintf(fid, '%% 时间戳: %s\n', timestamp);
    fprintf(fid, '%% 深度(nm)');
    
    % 写入能量标题
    for i = 1:length(chi_data.energies)
        fprintf(fid, '\tχ_real_norm_%.1feV\tχ_imag_norm_%.1feV', ...
                chi_data.energies(i), chi_data.energies(i));
    end
    fprintf(fid, '\n');
    
    % 写入归一化因子
    fprintf(fid, '%% 归一化因子(最大值): ');
    for i = 1:length(chi_data.energies)
        real_max_field = sprintf('chi_real_max_%d', i);
        imag_max_field = sprintf('chi_imag_max_%d', i);
        fprintf(fid, '\t%.6e\t%.6e', ...
                chi_data.(real_max_field), chi_data.(imag_max_field));
    end
    fprintf(fid, '\n');
    
    % 写入数据
    for i = 1:length(chi_data.z_grid)
        fprintf(fid, '%.6f', chi_data.z_grid(i));
        for j = 1:length(chi_data.energies)
            chi_real_field = sprintf('chi_real_normalized_%d', j);
            chi_imag_field = sprintf('chi_imag_normalized_%d', j);
            fprintf(fid, '\t%.6f\t%.6f', ...
                    chi_data.(chi_real_field)(i), chi_data.(chi_imag_field)(i));
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
    fprintf('  χ分布数据: %s\n', chi_filename);
    
    fprintf('TXT文件导出完成！\n');
end

function z_plot = create_stairs_depth_array(z_grid, dz)
% 为阶梯图创建深度数组
% 每个网格点扩展为网格边界
    z_plot = [];
    for i = 1:length(z_grid)
        z_start = z_grid(i) - dz/2;
        z_end = z_grid(i) + dz/2;
        if i == 1
            z_start = 0;  % 第一个点从0开始
        end
        z_plot = [z_plot; z_start; z_end];
    end
end

function data_plot = create_stairs_data(data)
% 为阶梯图创建数据数组
% 每个数据点重复两次以形成水平线段
    data_plot = [];
    for i = 1:length(data)
        data_plot = [data_plot; data(i); data(i)];
    end
end