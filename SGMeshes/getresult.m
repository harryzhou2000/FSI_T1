%% 从文本文件中导入数据
% 用于从以下文本文件中导入数据的脚本:
%
%    filename: D:\Documents\FEM\Mesh\poissionresult.txt
%
% 由 MATLAB 于 2020-05-16 18:00:08 自动生成

%% 设置导入选项
opts = delimitedTextImportOptions("NumVariables", 3);

% 指定范围和分隔符
opts.DataLines = [1, Inf];
opts.Delimiter = "\t";

% 指定列名称和类型
opts.VariableNames = ["VarName1", "VarName2", "VarName3"];
opts.VariableTypes = ["double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% 导入数据
result = readtable("D:\Documents\FEM\Mesh\poissionresult.txt", opts);

%% 转换为输出类型
result = table2array(result);

%% 清除临时变量
clear opts