function SG0012T1O2 = ImportExtracted(filename, dataLines)
%IMPORTFILE 从文本文件中导入数据
%  SG0012T1O2 = IMPORTFILE(FILENAME)读取文本文件 FILENAME 中默认选定范围的数据。
%  以表形式返回数据。
%
%  SG0012T1O2 = IMPORTFILE(FILE, DATALINES)按指定行间隔读取文本文件 FILENAME
%  中的数据。对于不连续的行间隔，请将 DATALINES 指定为正整数标量或 N×2 正整数标量数组。
%
%  示例:
%  SG0012T1O2 = importfile("E:\Documents\FSI_T1\Unstructured\UGSolver_Submit\dout\extracted\SG_0012_T1_O2.dat", [17, Inf]);
%
%  另请参阅 READTABLE。
%
% 由 MATLAB 于 2021-06-24 23:03:25 自动生成

%% 输入处理

% 如果不指定 dataLines，请定义默认范围
if nargin < 2
    dataLines = [17, Inf];
end

%% 设置导入选项
opts = delimitedTextImportOptions("NumVariables", 11);

% 指定范围和分隔符
opts.DataLines = dataLines;
opts.Delimiter = " ";

% 指定列名称和类型
opts.VariableNames = ["x", "y", "rho", "u", "v", "e", "um", "p", "a", "ma", "Var11"];
opts.SelectedVariableNames = ["x", "y", "rho", "u", "v", "e", "um", "p", "a", "ma"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string"];
opts = setvaropts(opts, 11, "WhitespaceRule", "preserve");
opts = setvaropts(opts, 11, "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% 导入数据
SG0012T1O2 = readtable(filename, opts);

end