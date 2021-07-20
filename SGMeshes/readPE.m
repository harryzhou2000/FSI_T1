function resultPE2 = readPE(filename, dataLines)
%IMPORTFILE ���ı��ļ��е�������
%  RESULTPE2 = IMPORTFILE(FILENAME)��ȡ�ı��ļ� FILENAME ��Ĭ��ѡ����Χ�����ݡ�  �Ա���ʽ�������ݡ�
%
%  RESULTPE2 = IMPORTFILE(FILE, DATALINES)��ָ���м����ȡ�ı��ļ� FILENAME
%  �е����ݡ����ڲ��������м�����뽫 DATALINES ָ��Ϊ������������ N��2 �������������顣
%
%  ʾ��:
%  resultPE2 = importfile("E:\Documents\FEM\Mesh\resultPE2.txt", [1, Inf]);
%
%  ������� READTABLE��
%
% �� MATLAB �� 2021-03-11 00:41:00 �Զ�����

%% ���봦��

% �����ָ�� dataLines���붨��Ĭ�Ϸ�Χ
if nargin < 2
    dataLines = [1, Inf];
end

%% ���õ���ѡ��
opts = delimitedTextImportOptions("NumVariables", 4);

% ָ����Χ�ͷָ���
opts.DataLines = dataLines;
opts.Delimiter = "\t";

% ָ�������ƺ�����
opts.VariableNames = ["ux", "uy", "x", "y"];
opts.VariableTypes = ["double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% ��������
resultPE2 = readtable(filename, opts);

end