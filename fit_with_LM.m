function density_params_opt = fit_with_LM(density_params_init, exp_data_multi, grid_model, input_params)
% FIT_WITH_LM 使用Levenberg-Marquardt算法拟合多波长数据
%   density_params_opt = FIT_WITH_LM(density_params_init, exp_data_multi, grid_model, input_params)
%
%   输入参数:
%   density_params_init - 初始密度参数向量
%   exp_data_multi - 多波长实验数据结构体数组
%   grid_model - 网格模型结构体
%   input_params - 输入参数结构体

    fprintf('=== 开始Levenberg-Marquardt多波长拟合 ===\n');
    fprintf('优化参数数量: %d\n', length(density_params_init));
    fprintf('实验数据组数: %d\n', length(exp_data_multi));
    
    % 设置优化选项
    options = optimoptions('lsqnonlin', ...
                         'Display', 'iter', ...
                         'MaxIterations', 15, ...
                         'MaxFunctionEvaluations', 1000000, ...
                         'FunctionTolerance', 1e-12, ...
                         'StepTolerance', 1e-12, ...
                         'Algorithm', 'levenberg-marquardt');
    
    % 设置参数边界（密度应为非负）
    lb = zeros(size(density_params_init));
    ub = [];
    
    % 定义目标函数
    objective = @(params) mf_wrapper_multi(params, exp_data_multi, grid_model, input_params);
    
    % 执行优化
    [density_params_opt, resnorm, residual, exitflag, output] = ...
        lsqnonlin(objective, density_params_init, lb, ub, options);
    
    % 输出结果
    fprintf('\n=== 优化完成 ===\n');
    fprintf('退出标志: %d\n', exitflag);
    fprintf('函数评估次数: %d\n', output.funcCount);
    fprintf('最终残差平方和: %.6e\n', resnorm);
end