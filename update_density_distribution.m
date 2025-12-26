function grid_model_updated = update_density_distribution(grid_model, input_params, density_params)
% UPDATE_DENSITY_DISTRIBUTION 根据优化参数更新密度分布

    grid_model_updated = grid_model;
    
    fit_elements = input_params.fit_elements;
    num_points = grid_model.num_points;
    num_fit_elements = length(fit_elements);
    
    % 检查参数数量
    expected_params = num_points * num_fit_elements;
    if length(density_params) ~= expected_params
        error('密度参数数量不匹配：期望 %d，实际 %d', expected_params, length(density_params));
    end
    
    % 按元素分组更新密度：[elem1_all_points, elem2_all_points, ...]
    for i = 1:num_fit_elements
        elem = fit_elements{i};
        start_idx = (i-1) * num_points + 1;
        end_idx = i * num_points;
        
        % 更新该元素的密度分布，并确保非负
        grid_model_updated.density_grid.(elem) = max(0, density_params(start_idx:end_idx));
    end
end