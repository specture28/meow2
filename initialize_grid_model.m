function grid_model = initialize_grid_model(input_params)
    fprintf('初始化网格模型和密度分布...\n');
    
    %% 所有可能的元素
    all_elements = {'Si', 'Co', 'Mg', 'O', 'C'};
    
    if input_params.use_saved_density
        %% 加载模式
        fprintf('从文件加载密度分布: %s\n', input_params.density_file);
        
        try
            loaded_data = load(input_params.density_file);
            density_data = loaded_data.density_data;
            
            % 提取网格和密度信息
            z_grid = density_data.z_grid;
            old_grid_step = z_grid(2) - z_grid(1);
            calculation_depth = max(z_grid);
            new_grid_step = input_params.grid_step;
            
            fprintf('原网格步长: %.3f nm, 新网格步长: %.3f nm\n', old_grid_step, new_grid_step);
            
            % 检查是否需要重采样
            if abs(old_grid_step - new_grid_step) > 1e-6
                fprintf('网格步长不同，开始重采样...\n');
                
                % 生成新网格
                z_grid_new = 0:new_grid_step:calculation_depth;
                
                % 对每个元素进行重采样
                density_grid = struct();
                for i = 1:length(all_elements)
                    elem = all_elements{i};
                    field_name = sprintf('density_%s', elem);
                    if isfield(density_data, field_name)
                        old_density = density_data.(field_name);
                        % 使用右端点策略重采样
                        new_density = resample_density_profile(z_grid, old_density, z_grid_new);
                        density_grid.(elem) = new_density;
                    else
                        density_grid.(elem) = zeros(size(z_grid_new));
                        fprintf('警告: 文件中未找到%s的密度数据，初始化为0\n', elem);
                    end
                end
                
                % 更新网格信息
                z_grid = z_grid_new;
                grid_step = new_grid_step;
                num_points = length(z_grid);
                
                fprintf('重采样完成: %d → %d 网格点\n', length(density_data.z_grid), num_points);
                
            else
                % 网格步长相同，直接使用原有数据
                fprintf('网格步长相同，直接加载数据\n');
                
                % 构建密度分布结构
                density_grid = struct();
                for i = 1:length(all_elements)
                    elem = all_elements{i};
                    field_name = sprintf('density_%s', elem);
                    if isfield(density_data, field_name)
                        density_grid.(elem) = density_data.(field_name);
                    else
                        % 如果文件中没有某个元素，初始化为0
                        density_grid.(elem) = zeros(size(z_grid));
                        fprintf('警告: 文件中未找到%s的密度数据，初始化为0\n', elem);
                    end
                end
                
                grid_step = old_grid_step;
                num_points = length(z_grid);
            end
            
            fprintf('成功加载密度分布，网格点数: %d\n', num_points);
            
        catch ME
            error('加载密度文件失败: %s', ME.message);
        end
        
    else
        %% 新建模式 - 根据层结构初始化
        fprintf('根据层结构生成初始密度分布\n');
        
        % 提取参数
        num_layers = input_params.num_layers;
        layer_materials = input_params.layer_materials;
        layer_thickness = input_params.layer_thickness;
        grid_step = input_params.grid_step;
        surface_buffer = input_params.surface_buffer;
        bulk_density = input_params.bulk_density;
        
        % 计算网格（只包含薄膜部分）
        total_layer_thickness = sum(layer_thickness);
        calculation_depth = surface_buffer + total_layer_thickness;
        z_grid = 0:grid_step:calculation_depth;
        num_points = length(z_grid);
        
        fprintf('薄膜计算深度: %.1f nm, 网格点数: %d\n', calculation_depth, num_points);
        
        % 初始化所有元素的密度分布为0
        density_grid = struct();
        for i = 1:length(all_elements)
            elem = all_elements{i};
            density_grid.(elem) = zeros(size(z_grid));
        end
        
        % 根据层结构分配密度
        fprintf('分配各层密度:\n');
        for i = 1:num_points
            z_pos = z_grid(i);
            
            if z_pos < surface_buffer
                % 空气层 - 所有密度保持0
                continue;
            else
                % 确定属于哪一层
                current_depth = surface_buffer;
                for layer = 1:num_layers
                    layer_end = current_depth + layer_thickness(layer);
                    if z_pos >= current_depth && z_pos < layer_end
                        % 在第layer层内，给对应材料赋值块体密度
                        material = layer_materials{layer};
                        density_grid.(material)(i) = bulk_density(material);
                        break;
                    end
                    current_depth = layer_end;
                end
            end
        end
        
        % 显示每层的设置
        current_depth = surface_buffer;
        for layer = 1:num_layers
            material = layer_materials{layer};
            density_val = bulk_density(material);
            fprintf('  第%d层 (%s): %.1f-%.1f nm, 密度=%.2f g/cm³\n', ...
                    layer, material, current_depth, current_depth+layer_thickness(layer), density_val);
            current_depth = current_depth + layer_thickness(layer);
        end
    end
    
    %% 基底信息（固定为单晶Si，不参与网格计算）
    substrate_info = struct();
    substrate_info.material = 'Si';
    substrate_info.density = input_params.substrate_density;  % 2.42 g/cm³
    substrate_info.start_depth = calculation_depth;
    
    %% 构建输出结构
    grid_model = struct();
    grid_model.z_grid = z_grid;  % 只包含薄膜网格
    grid_model.density_grid = density_grid;  % 包含所有5个元素的密度分布
    grid_model.num_points = num_points;
    grid_model.calculation_depth = calculation_depth;  % 薄膜结束深度
    grid_model.substrate = substrate_info;  % 基底作为半无限介质
    grid_model.dz = grid_step;
    grid_model.all_elements = all_elements;
    
    %% 验证拟合元素的密度分布
    fprintf('\n拟合元素密度统计:\n');
    for i = 1:length(input_params.fit_elements)
        elem = input_params.fit_elements{i};
        density_profile = density_grid.(elem);
        max_density = max(density_profile);
        nonzero_points = sum(density_profile > 1e-6);
        fprintf('  %s: 最大密度=%.2f g/cm³, 非零点数=%d\n', elem, max_density, nonzero_points);
    end
    
    fprintf('基底: %s (%.2f g/cm³, 固定, 不参与网格计算)\n', ...
            substrate_info.material, substrate_info.density);
    fprintf('薄膜密度分布初始化完成！\n');
end

function new_density = resample_density_profile(old_z, old_density, new_z)
% 使用右端点策略重采样密度分布
    new_density = zeros(size(new_z));
    old_dz = old_z(2) - old_z(1);
    
    for i = 1:length(new_z)
        z_point = new_z(i);
        
        % 找到该点属于原网格的哪个区间
        % 原网格区间: [old_z(j), old_z(j+1))
        grid_index = floor(z_point / old_dz) + 1;
        
        % 边界处理：如果恰好在边界上，使用右端点策略
        if abs(mod(z_point, old_dz)) < 1e-10 && grid_index > 1
            grid_index = grid_index; % 取右侧区间的密度
        end
        
        % 确保索引在有效范围内
        if grid_index >= 1 && grid_index <= length(old_density)
            new_density(i) = old_density(grid_index);
        end
    end
end