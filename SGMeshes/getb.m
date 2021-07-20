function b = getb(filename, dataLines)
%IMPORTFILE ���ı��ļ��е�������
%  B = IMPORTFILE(FILENAME)��ȡ�ı��ļ� FILENAME ��Ĭ��ѡ����Χ�����ݡ�  �Ա���ʽ�������ݡ�
%
%  B = IMPORTFILE(FILE, DATALINES)��ָ���м����ȡ�ı��ļ� FILENAME
%  �е����ݡ����ڲ��������м�����뽫 DATALINES ָ��Ϊ������������ N��2 �������������顣
%
%  ʾ��:
%  b = importfile("E:\Documents\FEM\Mesh\PE2_b.txt", [1, Inf]);
%
%  ������� READTABLE��
%
% �� MATLAB �� 2021-03-11 01:13:19 �Զ�����

%% ���봦��

% �����ָ�� dataLines���붨��Ĭ�Ϸ�Χ
if nargin < 2
    dataLines = [1, Inf];
end

%% ���õ���ѡ��
opts = delimitedTextImportOptions("NumVariables", 1);

% ָ����Χ�ͷָ���
opts.DataLines = dataLines;
opts.Delimiter = ",";

% ָ�������ƺ�����
opts.VariableNames = "VarName1";
opts.VariableTypes = "double";
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% ��������
b = readtable(filename, opts);

end