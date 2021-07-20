function SG0012T1O2 = ImportExtracted(filename, dataLines)
%IMPORTFILE ���ı��ļ��е�������
%  SG0012T1O2 = IMPORTFILE(FILENAME)��ȡ�ı��ļ� FILENAME ��Ĭ��ѡ����Χ�����ݡ�
%  �Ա���ʽ�������ݡ�
%
%  SG0012T1O2 = IMPORTFILE(FILE, DATALINES)��ָ���м����ȡ�ı��ļ� FILENAME
%  �е����ݡ����ڲ��������м�����뽫 DATALINES ָ��Ϊ������������ N��2 �������������顣
%
%  ʾ��:
%  SG0012T1O2 = importfile("E:\Documents\FSI_T1\Unstructured\UGSolver_Submit\dout\extracted\SG_0012_T1_O2.dat", [17, Inf]);
%
%  ������� READTABLE��
%
% �� MATLAB �� 2021-06-24 23:03:25 �Զ�����

%% ���봦��

% �����ָ�� dataLines���붨��Ĭ�Ϸ�Χ
if nargin < 2
    dataLines = [17, Inf];
end

%% ���õ���ѡ��
opts = delimitedTextImportOptions("NumVariables", 11);

% ָ����Χ�ͷָ���
opts.DataLines = dataLines;
opts.Delimiter = " ";

% ָ�������ƺ�����
opts.VariableNames = ["x", "y", "rho", "u", "v", "e", "um", "p", "a", "ma", "Var11"];
opts.SelectedVariableNames = ["x", "y", "rho", "u", "v", "e", "um", "p", "a", "ma"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string"];
opts = setvaropts(opts, 11, "WhitespaceRule", "preserve");
opts = setvaropts(opts, 11, "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% ��������
SG0012T1O2 = readtable(filename, opts);

end