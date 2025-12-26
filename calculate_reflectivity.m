function R = calculate_reflectivity(grid_model, input_params, energy_index, theta)
% CALCULATE_REFLECTIVITY 基于密度分布计算X射线反射率
%   R = CALCULATE_REFLECTIVITY(grid_model, input_params, energy_index, theta)
%
%   输入参数:
%   grid_model - 网格模型结构体（包含密度分布）
%   input_params - 输入参数结构体（包含ASF数据等）
%   energy_index - 能量点索引
%   theta - 入射角度数组(度)
%
%   输出:
%   R - 反射率数组
%
%   计算流程:
%   密度分布 → χ值分布 → 层结构 → Parratt算法 → 反射率

    %% 1. 将密度分布转换为χ值分布
    [chi_real, chi_imag] = density_to_chi(grid_model, input_params, energy_index);
    
    % 合成复数χ
    chi_values = chi_real + 1i * chi_imag;
    
    %% 2. 获取波长
    wavelength = input_params.wavelengths(energy_index);  % nm
    
    %% 3. 构建层结构
    dz = grid_model.dz;
    n = length(chi_values);
    
    layers = struct('chi', cell(n, 1), 'thickness', cell(n, 1));
    
    for i = 1:n
        layers(i).chi = chi_values(i);
        layers(i).thickness = dz;
    end
    
    %% 4. 调用Parratt算法计算反射率
    R = parratt(theta, wavelength, layers);
    
    %% 调试信息（可选）
    if any(isnan(R)) || any(R < 0) || any(R > 1)
        warning('calculate_reflectivity: 反射率数值异常');
        fprintf('  χ统计: real=[%.2e, %.2e], imag=[%.2e, %.2e]\n', ...
                min(chi_real), max(chi_real), min(chi_imag), max(chi_imag));
        fprintf('  R统计: [%.2e, %.2e], NaN个数=%d\n', ...
                min(R), max(R), sum(isnan(R)));
    end
end