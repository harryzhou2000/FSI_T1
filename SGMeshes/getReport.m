function PE2Report = getReport(filename, dataLines)
%IMPORTFILE ���ı��ļ��е�������
%  PE2REPORT = IMPORTFILE(FILENAME)��ȡ�ı��ļ� FILENAME ��Ĭ��ѡ����Χ�����ݡ�  �Ա���ʽ�������ݡ�
%
%  PE2REPORT = IMPORTFILE(FILE, DATALINES)��ָ���м����ȡ�ı��ļ� FILENAME
%  �е����ݡ����ڲ��������м�����뽫 DATALINES ָ��Ϊ������������ N��2 �������������顣
%
%  ʾ��:
%  PE2Report = importfile("E:\Documents\FEM\Mesh\PE2_Report.txt", [1, Inf]);
%
%  ������� READTABLE��
%
% �� MATLAB �� 2021-03-11 11:56:00 �Զ�����

%% ���봦��

% �����ָ�� dataLines���붨��Ĭ�Ϸ�Χ
if nargin < 2
    dataLines = [1, Inf];
end

%% ���õ���ѡ��
opts = delimitedTextImportOptions("NumVariables", 13);

% ָ����Χ�ͷָ���
opts.DataLines = dataLines;
opts.Delimiter = "\t";

% ָ�������ƺ�����
opts.VariableNames = ["x1", "x2", "x3", "y1", "y2", "y3", "u1", "u2", "u3", "v1", "v2", "v3", "vms"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% ��������
PE2Report = readtable(filename, opts);

end