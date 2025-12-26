function input_params = input_parameters()
    fprintf('=== 多波长部分密度拟合参数输入 ===\n');
    
    %% 物理常数
    re = 2.82e-13;  % 经典电子半径 (cm)
    NA = 6.02214076e23;  % 阿伏伽德罗常数 (atoms/mol)
    
    %% 能量和波长对应关系
    energies = [90, 220, 410, 640, 1125, 1625];  % eV
    wavelengths = [13.78, 5.636, 3.024, 1.937, 1.102, 0.763];  % nm
    fprintf('使用%d个波长进行多波长拟合\n', length(wavelengths));
    
    %% ASF数据表（f1, f2对应6个能量点）
    asf_data = struct();
    % Si
    asf_data.Si.f1 = [0.581, 9.63, 13.12, 13.36, 12.78, 11.66];
    asf_data.Si.f2 = [0.436, 8.56, 4.43, 2.26, 0.859, 0.461];
    % Co
    asf_data.Co.f1 = [8.83, 15.74, 17.31, 15.03, 23.71, 26.62];
    asf_data.Co.f2 = [8.89, 7.13, 4.21, 2.61, 11.56, 6.81];
    % Mg
    asf_data.Mg.f1 = [2.24, 10.04, 11.07, 10.83, 9.54, 11.19];
    asf_data.Mg.f2 = [8.11, 5.25, 2.29, 1.11, 0.333, 3.23];
    % O
    asf_data.O.f1 = [5.93, 6.57, 5.81, 6.62, 8.32, 8.36];
    asf_data.O.f2 = [2.81, 0.919, 0.352, 3.49, 1.45, 0.772];
    % C
    asf_data.C.f1 = [4.26, 3.59, 5.59, 6.28, 6.31, 6.22];
    asf_data.C.f2 = [0.784, 0.245, 2.62, 1.34, 0.512, 0.257];
    
    %% 原子质量和块体密度数据库（薄膜材料）
    atomic_mass = containers.Map({'Si', 'Co', 'Mg', 'O', 'C'}, ...
                                {28.085, 58.933, 24.305, 15.999, 12.011});  % g/mol
    bulk_density = containers.Map({'Si', 'Co', 'Mg', 'O', 'C'}, ...
                                 {2.33, 8.86, 1.74, 1.43, 2.27});  % g/cm³
    
    %% 基底参数（固定）
    substrate_density = 2.33;  % 单晶硅密度 g/cm³
    
    %% 选择初始化方式
    use_saved = input('是否从文件加载初始密度分布？(y/n): ', 's');
    use_saved_density = strcmpi(use_saved, 'y');
    
    if use_saved_density
        %% 加载模式
        [density_file, density_path] = uigetfile('*.mat', '选择密度分布文件');
        if isequal(density_file, 0)
            error('未选择文件，程序退出');
        end
        density_file = fullfile(density_path, density_file);
        fprintf('已选择文件: %s\n', density_file);
        
        %% 网格参数输入（新增部分）
        fprintf('\n=== 网格参数 ===\n');
        grid_step = input('新的网格步长(nm, 默认0.1): ');
        if isempty(grid_step), grid_step = 0.1; end
        
    else
        %% 新建模式 - 输入层结构
        fprintf('\n=== 层结构输入 ===\n');
        fprintf('可用薄膜材料: Si, Co, Mg, O, C\n');
        fprintf('基底固定为单晶硅(密度=%.2f g/cm³)\n', substrate_density);
        
        num_layers = input('薄膜层数: ');
        
        layer_materials = cell(num_layers, 1);
        layer_thickness = zeros(num_layers, 1);
        
        for i = 1:num_layers
            fprintf('\n--- 第 %d 层 ---\n', i);
            while true
                material = input(sprintf('第%d层材料: ', i), 's');
                if isKey(bulk_density, material)
                    layer_materials{i} = material;
                    fprintf('选择%s (块体密度: %.2f g/cm³)\n', material, bulk_density(material));
                    break;
                else
                    fprintf('材料 "%s" 不可用，请从[Si, Co, Mg, O, C]中选择\n', material);
                end
            end
            layer_thickness(i) = input(sprintf('第%d层厚度(nm): ', i));
        end
        
        %% 网格参数
        fprintf('\n=== 网格参数 ===\n');
        grid_step = input('网格步长(nm, 默认0.1): ');
        if isempty(grid_step), grid_step = 0.1; end
        surface_buffer = input('表面缓冲(nm, 默认2.0): ');
        if isempty(surface_buffer), surface_buffer = 2.0; end
    end
    
    %% 选择拟合元素
    available_elements = {'Si', 'Co', 'Mg', 'O', 'C'};
    fprintf('\n=== 拟合元素选择 ===\n');
    fprintf('可用元素: %s\n', strjoin(available_elements, ', '));
    fprintf('注意: 基底Si密度固定，不参与拟合\n');
    fit_elements_str = input('选择用于拟合的元素(空格分隔): ', 's');
    fit_elements = strsplit(fit_elements_str);
    
    % 验证选择的元素
    for i = 1:length(fit_elements)
        if ~ismember(fit_elements{i}, available_elements)
            error('元素 "%s" 不可用', fit_elements{i});
        end
    end
    
    %% 密度优化边界
    fprintf('\n=== 密度边界设置 ===\n');
    density_bounds = containers.Map();
    for i = 1:length(fit_elements)
        elem = fit_elements{i};
        bulk_dens = bulk_density(elem);
        fprintf('\n%s 块体密度: %.2f g/cm³\n', elem, bulk_dens);
        lower = input(sprintf('下界(默认0): '));
        upper = input(sprintf('上界(默认%.2f): ', bulk_dens*1.1));
        if isempty(lower), lower = 0; end
        if isempty(upper), upper = bulk_dens * 1.1; end
        density_bounds(elem) = [lower, upper];
    end
    
    %% 平滑性约束参数
    fprintf('\n=== 平滑性约束设置 ===\n');
    Q = input('输入Q值(默认1): ');
    if isempty(Q), Q = 1; end
    
    %% 基底约束设置 ← 新增
    fprintf('\n=== 基底约束设置 ===\n');
    fprintf('基底Si密度: %.2f g/cm³\n', substrate_density);
    
    % 检查Si是否在拟合元素中
    if ismember('Si', fit_elements)
        fprintf('检测到Si在拟合元素中，可以启用基底约束\n');
        enable_substrate_constraint = input('是否启用基底约束？(y/n, 默认y): ', 's');
        if isempty(enable_substrate_constraint), enable_substrate_constraint = 'y'; end
        
        if strcmpi(enable_substrate_constraint, 'y')
            num_substrate_constraint_layers = input('约束最后几层？(默认2): ');
            if isempty(num_substrate_constraint_layers), num_substrate_constraint_layers = 2; end
            
            substrate_constraint_weight = input('基底约束权重(默认100): ');
            if isempty(substrate_constraint_weight), substrate_constraint_weight = 100; end
            
            fprintf('✅ 基底约束已启用: 最后%d层Si密度 → %.2f g/cm³ (权重=%.1f)\n', ...
                    num_substrate_constraint_layers, substrate_density, substrate_constraint_weight);
        else
            num_substrate_constraint_layers = 0;
            substrate_constraint_weight = 0;
            fprintf('❌ 未启用基底约束\n');
        end
    else
        fprintf('Si不在拟合元素中，自动禁用基底约束\n');
        num_substrate_constraint_layers = 0;
        substrate_constraint_weight = 0;
    end
    
    %% 输出结构
    input_params = struct();
    input_params.re = re;
    input_params.NA = NA;
    input_params.energies = energies;
    input_params.wavelengths = wavelengths;
    input_params.asf_data = asf_data;
    input_params.atomic_mass = atomic_mass;
    input_params.bulk_density = bulk_density;
    input_params.substrate_density = substrate_density;  % 基底密度
    input_params.fit_elements = fit_elements;
    input_params.density_bounds = density_bounds;
    input_params.use_saved_density = use_saved_density;
    input_params.Q = Q;  % 平滑性权重参数
    
    % 基底约束参数 ← 新增
    input_params.num_substrate_constraint_layers = num_substrate_constraint_layers;
    input_params.substrate_constraint_weight = substrate_constraint_weight;
    
    if use_saved_density
        input_params.density_file = density_file;
        input_params.grid_step = grid_step;  % 新增
    else
        input_params.num_layers = num_layers;
        input_params.layer_materials = layer_materials;
        input_params.layer_thickness = layer_thickness;
        input_params.grid_step = grid_step;
        input_params.surface_buffer = surface_buffer;
    end
    
    %% 确认显示
    fprintf('\n=== 参数确认 ===\n');
    fprintf('拟合元素: %s\n', strjoin(fit_elements, ', '));
    fprintf('基底: 单晶硅(%.2f g/cm³, 固定)\n', substrate_density);
    fprintf('平滑性权重 Q: %.1f\n', Q);
    
    % 显示基底约束信息
    if substrate_constraint_weight > 0
        fprintf('基底约束: 最后%d层Si → %.2f g/cm³ (权重=%.1f)\n', ...
                num_substrate_constraint_layers, substrate_density, substrate_constraint_weight);
    else
        fprintf('基底约束: 未启用\n');
    end
    
    if use_saved_density
        fprintf('加载文件: %s\n', density_file);
        fprintf('网格步长: %.3f nm\n', grid_step);
    else
        fprintf('薄膜层结构:\n');
        for i = 1:num_layers
            fprintf('  第%d层: %s, %.1f nm\n', i, layer_materials{i}, layer_thickness(i));
        end
        fprintf('薄膜总厚度: %.1f nm\n', sum(layer_thickness));
        fprintf('网格步长: %.3f nm\n', grid_step);
    end
    fprintf('使用%d个波长进行拟合\n', length(wavelengths));
    
    fprintf('参数输入完成！\n');
end