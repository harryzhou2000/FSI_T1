%% ���ı��ļ��е�������
% ���ڴ������ı��ļ��е������ݵĽű�:
%
%    filename: D:\Documents\FEM\Mesh\poissionresult.txt
%
% �� MATLAB �� 2020-05-16 18:00:08 �Զ�����

%% ���õ���ѡ��
opts = delimitedTextImportOptions("NumVariables", 3);

% ָ����Χ�ͷָ���
opts.DataLines = [1, Inf];
opts.Delimiter = "\t";

% ָ�������ƺ�����
opts.VariableNames = ["VarName1", "VarName2", "VarName3"];
opts.VariableTypes = ["double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% ��������
result = readtable("D:\Documents\FEM\Mesh\poissionresult.txt", opts);

%% ת��Ϊ�������
result = table2array(result);

%% �����ʱ����
clear opts