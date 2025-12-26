function residuals = mf_wrapper_multi(density_params, exp_data_multi, grid_model, input_params)
% MF_WRAPPER_MULTI 多波长拟合的目标函数
%   计算所有波长的拟合残差 + 平滑性惩罚 + 基底约束

    % 硬编码的权重数组（可手动调节）
    weights = [7, 2.6, 2.5, 1, 1, 2];  % 6个波长的权重
    
    % 确保权重数量匹配
    if length(weights) ~= length(exp_data_multi)
        error('权重数组长度与数据组数不匹配');
    end
    
    % 更新密度分布
    grid_model_updated = update_density_distribution(grid_model, input_params, density_params);
    
    % 1. 计算所有波长的反射率残差
    all_residuals = [];
    
    for i = 1:length(exp_data_multi)
        % 获取当前波长的数据
        theta_exp = exp_data_multi(i).theta;
        R_exp = exp_data_multi(i).R;
        energy_index = exp_data_multi(i).energy_index;
        weight = weights(i);
        
        % 计算理论反射率
        R_calc = calculate_reflectivity(grid_model_updated, input_params, energy_index, theta_exp);
        
        % 计算加权残差（线性形式）
        R_residuals = weight * (R_exp - R_calc) ./ R_exp;
        
        % 确保是列向量并累加
        all_residuals = [all_residuals; R_residuals(:)];
    end
    
    % 2. 计算平滑性残差（线性形式）
    if input_params.Q > 0
        smoothness_residuals = calculate_smoothness_residuals(density_params, grid_model, input_params);
        all_residuals = [all_residuals; smoothness_residuals];
    end
    
    % 3. 计算基底约束残差 ← 新增
    if input_params.substrate_constraint_weight > 0
        substrate_residuals = calculate_substrate_constraint_residuals(density_params, grid_model, input_params);
        all_residuals = [all_residuals; substrate_residuals];
    end
    
    residuals = all_residuals;
end

function smoothness_residuals = calculate_smoothness_residuals(density_params, grid_model, input_params)
% 计算平滑性残差（基于相邻梯度差，即二阶导数形式）
    smoothness_residuals = [];
    num_points = grid_model.num_points;
    dz = grid_model.dz;  % 网格步长
    
    for elem_idx = 1:length(input_params.fit_elements)
        % 提取该元素的密度参数向量
        start_idx = (elem_idx-1) * num_points + 1;
        end_idx = elem_idx * num_points;
        densities = density_params(start_idx:end_idx);
        
        % 计算相邻梯度差：dρ/dz[i+1] - dρ/dz[i]
        % 其中 dρ/dz[i] = (ρ[i] - ρ[i-1]) / h （后向差分）
        % 相邻梯度差 = (ρ[i+1] - 2ρ[i] + ρ[i-1]) / h
        for i = 1:(length(densities)-1)
            if i == 1
                rho_prev = 0;  % 空气密度 ρ[0] = 0
            else
                rho_prev = densities(i-1);
            end
            
            % 计算二阶导数形式的梯度差
            gradient_diff = (densities(i+1) - 2*densities(i) + rho_prev) / dz;
            gradient_diff_residual = sqrt(input_params.Q) * gradient_diff;
            smoothness_residuals = [smoothness_residuals; gradient_diff_residual];
        end
    end
end

function substrate_residuals = calculate_substrate_constraint_residuals(density_params, grid_model, input_params)
% 计算基底约束残差：强制最后几层的Si密度接近基底密度
    substrate_residuals = [];
    
    % 约束参数
    num_constrained_layers = input_params.num_substrate_constraint_layers; % 约束最后几层
    substrate_weight = sqrt(input_params.substrate_constraint_weight);     % 约束权重
    target_density = input_params.substrate_density;                       % 目标基底密度
    
    % 找到Si元素在拟合元素中的位置
    si_elem_idx = find(strcmp(input_params.fit_elements, 'Si'));
    
    if isempty(si_elem_idx)
        % 如果Si不在拟合元素中，不施加约束
        fprintf('警告: Si不在拟合元素中，跳过基底约束\n');
        return;
    end
    
    % 提取Si元素的密度参数
    num_points = grid_model.num_points;
    start_idx = (si_elem_idx-1) * num_points + 1;
    end_idx = si_elem_idx * num_points;
    si_densities = density_params(start_idx:end_idx);
    
    % 对最后几层施加约束
    for i = 1:num_constrained_layers
        layer_idx = num_points - i + 1; % 从最后一层开始倒数
        if layer_idx >= 1
            % 约束：该层Si密度应接近基底密度
            density_deviation = si_densities(layer_idx) - target_density;
            constraint_residual = substrate_weight * density_deviation;
            substrate_residuals = [substrate_residuals; constraint_residual];
        end
    end
    
    % 调试信息（可选）
    % if ~isempty(substrate_residuals)
    %     fprintf('基底约束: 约束最后%d层Si密度 → %.2f g/cm³ (权重=√%.1f=%.2f)\n', ...
    %             num_constrained_layers, target_density, input_params.substrate_constraint_weight, substrate_weight);
    %     fprintf('  当前最后几层Si密度: ');
    %     for i = 1:min(num_constrained_layers, 3)  % 显示前3层
    %         layer_idx = num_points - i + 1;
    %         if layer_idx >= 1
    %             fprintf('%.3f ', si_densities(layer_idx));
    %         end
    %     end
    %     fprintf('g/cm³\n');
    % end
end