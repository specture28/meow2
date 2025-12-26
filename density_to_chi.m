function [chi_real, chi_imag] = density_to_chi(grid_model, input_params, energy_index)
% DENSITY_TO_CHI 将密度分布转换为极化率分布
%   [chi_real, chi_imag] = DENSITY_TO_CHI(grid_model, input_params, energy_index)
%   
%   使用严谨公式直接计算极化率：
%   χ(λ) = (re/π) × λ² × Σ nⱼ f⁽ʲ⁾(λ)
%   其中 f⁽ʲ⁾(λ) = f₁⁽ʲ⁾(λ) - if₂⁽ʲ⁾(λ)

    %% 提取参数
    density_grid = grid_model.density_grid;
    all_elements = grid_model.all_elements;
    num_points = grid_model.num_points;
    
    % 物理常数
    re = input_params.re;           % 经典电子半径 (cm)
    NA = input_params.NA;           % 阿伏伽德罗常数
    
    % 当前能量的参数
    wavelength = input_params.wavelengths(energy_index);  % nm
    lambda_cm = wavelength * 1e-7;  % 转换为 cm
    
    asf_data = input_params.asf_data;
    atomic_mass = input_params.atomic_mass;
    
    %% 公式中的公共因子
    % χ(λ) = (re/π) × λ² × Σ nⱼ f⁽ʲ⁾(λ)
    prefactor = (re / pi) * lambda_cm^2;
    
    %% 初始化输出数组
    chi_real = zeros(num_points, 1);
    chi_imag = zeros(num_points, 1);
    
    %% 遍历每个网格点
    for i = 1:num_points
        
        chi_point_real = 0;
        chi_point_imag = 0;
        
        % 遍历所有元素
        for j = 1:length(all_elements)
            elem = all_elements{j};
            
            % 获取该元素在该点的密度 (g/cm³)
            rho_j = density_grid.(elem)(i);
            
            if rho_j > 1e-10  % 只计算有意义的密度
                
                % 获取原子质量 μⱼ (g/mol)
                mu_j = atomic_mass(elem);
                
                % 计算数密度 nⱼ (atoms/cm³)
                n_j = rho_j * NA / mu_j;
                
                % 获取该元素在当前能量的ASF
                f1_j = asf_data.(elem).f1(energy_index);
                f2_j = asf_data.(elem).f2(energy_index);
                
                % f⁽ʲ⁾(λ) = f₁⁽ʲ⁾(λ) - if₂⁽ʲ⁾(λ)
                % 所以：
                % Re[f⁽ʲ⁾] = f₁⁽ʲ⁾
                % Im[f⁽ʲ⁾] = -f₂⁽ʲ⁾
                
                % 累加该元素对χ的贡献
                chi_point_real = chi_point_real + n_j * f1_j;
                chi_point_imag = chi_point_imag + n_j * (-f2_j);
            end
        end
        
        % 应用公共因子
        chi_real(i) = prefactor * chi_point_real;
        chi_imag(i) = prefactor * chi_point_imag;
    end
    
    %% 调试信息（可选）
    if any(isnan(chi_real)) || any(isnan(chi_imag))
        warning('density_to_chi: 检测到NaN值，请检查输入数据');
    end
    
    % 显示统计信息
    if nargout == 0
        fprintf('density_to_chi 转换完成 (直接公式):\n');
        fprintf('  能量索引: %d (%.1f eV, %.3f nm)\n', ...
                energy_index, input_params.energies(energy_index), wavelength);
        fprintf('  λ = %.3f nm, λ² = %.2e cm²\n', wavelength, lambda_cm^2);
        fprintf('  公共因子 (re/π)λ² = %.2e\n', prefactor);
        fprintf('  χ_real 范围: %.2e ~ %.2e\n', min(chi_real), max(chi_real));
        fprintf('  χ_imag 范围: %.2e ~ %.2e\n', min(chi_imag), max(chi_imag));
    end
end