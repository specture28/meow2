function main_multi_wavelength_fitting()
% MAIN_MULTI_WAVELENGTH_FITTING 多波长X射线反射率拟合主程序
%   完整的工作流程：参数输入 → 初始化 → 数据加载 → 拟合 → 结果展示

    clc; clear; close all;
    
    fprintf('================================================\n');
    fprintf('    多波长X射线反射率部分密度拟合程序\n');
    fprintf('================================================\n');
    
    try
        %% 步骤1：参数输入
        fprintf('\n【步骤1/6】参数输入\n');
        input_params = input_parameters();
        
        %% 步骤2：初始化网格模型
        fprintf('\n【步骤2/6】初始化网格模型\n');
        grid_model = initialize_grid_model(input_params);
        
        %% 步骤3：加载实验数据
        fprintf('\n【步骤3/6】加载多波长实验数据\n');
        exp_data_multi = load_experimental_data_multi(input_params);
        
        %% 步骤4：构造初始密度参数向量
        fprintf('\n【步骤4/6】构造初始密度参数\n');
        density_params_init = [];
        for i = 1:length(input_params.fit_elements)
            elem = input_params.fit_elements{i};
            initial_density = grid_model.density_grid.(elem);
            density_params_init = [density_params_init; initial_density(:)];
        end
        fprintf('初始参数向量长度: %d\n', length(density_params_init));
        
        %% 步骤5：执行拟合
        fprintf('\n【步骤5/6】开始多波长拟合\n');
        density_params_opt = fit_with_LM(density_params_init, exp_data_multi, grid_model, input_params);
        
        %% 步骤6：生成完整结果报告 (统一输出格式)
        fprintf('\n【步骤6/6】生成完整结果报告 (包含MF_R评估)\n');
        generate_complete_fitting_results(density_params_init, density_params_opt, exp_data_multi, grid_model, input_params);

        fprintf('\n================================================\n');
        fprintf('           多波长拟合程序运行完成！\n');
        fprintf('================================================\n');
        
    catch ME
        fprintf('\n❌ 程序执行出错：%s\n', ME.message);
        fprintf('错误位置：%s (第%d行)\n', ME.stack(1).name, ME.stack(1).line);
        fprintf('请检查输入参数和数据文件\n');
    end
end