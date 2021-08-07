function PE2Report = getReport(filename, dataLines)
%IMPORTFILE 从文本文件中导入数据
%  PE2REPORT = IMPORTFILE(FILENAME)读取文本文件 FILENAME 中默认选定范围的数据。  以表形式返回数据。
%
%  PE2REPORT = IMPORTFILE(FILE, DATALINES)按指定行间隔读取文本文件 FILENAME
%  中的数据。对于不连续的行间隔，请将 DATALINES 指定为正整数标量或 N×2 正整数标量数组。
%
%  示例:
%  PE2Report = importfile("E:\Documents\FEM\Mesh\PE2_Report.txt", [1, Inf]);
%
%  另请参阅 READTABLE。
%
% 由 MATLAB 于 2021-03-11 11:56:00 自动生成

%% 输入处理

% 如果不指定 dataLines，请定义默认范围
if nargin < 2
    dataLines = [1, Inf];
end

%% 设置导入选项
opts = delimitedTextImportOptions("NumVariables", 13);

% 指定范围和分隔符
opts.DataLines = dataLines;
opts.Delimiter = "\t";

% 指定列名称和类型
opts.VariableNames = ["x1", "x2", "x3", "y1", "y2", "y3", "u1", "u2", "u3", "v1", "v2", "v3", "vms"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% 导入数据
PE2Report = readtable(filename, opts);

end