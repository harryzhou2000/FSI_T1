function resultPE2Eigen = getReportEigen(filename, dataLines)
%IMPORTFILE ���ı��ļ��е�������
%  RESULTPE2EIGEN = IMPORTFILE(FILENAME)��ȡ�ı��ļ� FILENAME ��Ĭ��ѡ����Χ�����ݡ�
%  �Ա���ʽ�������ݡ�
%
%  RESULTPE2EIGEN = IMPORTFILE(FILE, DATALINES)��ָ���м����ȡ�ı��ļ� FILENAME
%  �е����ݡ����ڲ��������м�����뽫 DATALINES ָ��Ϊ������������ N��2 �������������顣
%
%  ʾ��:
%  resultPE2Eigen = importfile("E:\VS 03\StructrualGrid\SGSolverA\Mesh\resultPE2_Eigen.txt", [2, Inf]);
%
%  ������� READTABLE��
%
% �� MATLAB �� 2021-04-03 16:58:53 �Զ�����

%% ���봦��

% �����ָ�� dataLines���붨��Ĭ�Ϸ�Χ
if nargin < 2
    dataLines = [2, Inf];
end

%% ���õ���ѡ��
opts = delimitedTextImportOptions("NumVariables", 12);

% ָ����Χ�ͷָ���
opts.DataLines = dataLines;
opts.Delimiter = "\t";

% ָ�������ƺ�����
opts.VariableNames = ["x1", "x2", "x3", "y1", "y2", "y3", "u1", "u2", "u3", "v1", "v2", "v3"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% ��������
resultPE2Eigen = readtable(filename, opts);

end