function PE_A = getA(filename, dataLines)
%IMPORTFILE ���ı��ļ��е�������
%  PE_A = IMPORTFILE(FILENAME)��ȡ�ı��ļ� FILENAME ��Ĭ��ѡ����Χ�����ݡ�  ������ֵ���ݡ�
%
%  PE_A = IMPORTFILE(FILE, DATALINES)��ָ���м����ȡ�ı��ļ� FILENAME
%  �е����ݡ����ڲ��������м�����뽫 DATALINES ָ��Ϊ������������ N��2 �������������顣
%
%  ʾ��:
%  PE_A = importfile("E:\Documents\FEM\Mesh\PE2_A.txt", [1, Inf]);
%
%  ������� READTABLE��
%
% �� MATLAB �� 2021-03-11 01:07:28 �Զ�����

%% ���봦��

% �����ָ�� dataLines���붨��Ĭ�Ϸ�Χ
if nargin < 2
    dataLines = [1, Inf];
end

%% ���õ���ѡ��
opts = delimitedTextImportOptions("NumVariables", 3);

% ָ����Χ�ͷָ���
opts.DataLines = dataLines;
opts.Delimiter = "\t";

% ָ�������ƺ�����
opts.VariableNames = ["VarName1", "VarName2", "e09"];
opts.VariableTypes = ["double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% ��������
PE_A = readtable(filename, opts);

%% ת��Ϊ�������
PE_A = table2array(PE_A);
end