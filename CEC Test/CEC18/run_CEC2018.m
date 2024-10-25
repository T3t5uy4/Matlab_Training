% cec2018���Ժ���-��ע΢�Ź��ںš��Ż��㷨����
clc;clear;close all
%%
nPop=50; % ��Ⱥ��
Max_iter=300; % ����������
dim = 10; % ά�ȣ���ѡ 2, 10, 20��30��50��100

%%  ѡ����
for Function_name=1:30 % �������� 1 - 30

    % lb->���ޣ�ub->���ޣ�fobj->Ŀ�꺯��
    lb=-100*ones(1,dim);
    ub=100*ones(1,dim);
    fobj = @(x) cec18_func(x',Function_name); % ����cec2018 ���Ժ���
    
    %% �����㷨-����ʽһ����
    Optimal_results={}; % ������浽 Optimal results
    index = 1;
    % WOA
    tic
    [Best_score,Best_x,cg_curve]=WOA(nPop,Max_iter,lb,ub,dim,fobj);
    Optimal_results{1,index}="WOA";         % �㷨����
    Optimal_results{2,index}=cg_curve;      % ��������
    Optimal_results{3,index}=Best_score;   % ���ź���ֵ
    Optimal_results{4,index}=Best_x;          % ���ű���
    Optimal_results{5,index}=toc;               % ����ʱ��
    index = index +1;
    % HHO
    tic
    [Best_score,Best_x,cg_curve]=HHO(nPop,Max_iter,lb,ub,dim,fobj);
    Optimal_results{1,index}="HHO";
    Optimal_results{2,index}=cg_curve;
    Optimal_results{3,index}=Best_score;
    Optimal_results{4,index}=Best_x;
    Optimal_results{5,index}=toc;
    index = index +1;
    % GWO
    tic
    [Best_score,Best_x,cg_curve]=GWO(nPop,Max_iter,lb,ub,dim,fobj);
    Optimal_results{1,index}="GWO";
    Optimal_results{2,index}=cg_curve;
    Optimal_results{3,index}=Best_score;
    Optimal_results{4,index}=Best_x;
    Optimal_results{5,index}=toc;
    index = index +1;
    % ���и����㷨����������ʽ�ڴ˴����
    
    % 
    %% plot ��������
    figure(Function_name)
    for i = 1:size(Optimal_results, 2)
        %     plot(Optimal_results{2, i},'Linewidth',2)
        semilogy(Optimal_results{2, i},'Linewidth',2)
        hold on
    end
    title(['CEC2018 Convergence curve, Dim=' num2str(dim)])
    xlabel('Iteration');
    ylabel(['Best score on F' num2str(Function_name) ]);
    axis tight
    grid on;box on
    set(gcf,'Position',[400 200 400 250])
    legend(Optimal_results{1, :})
    saveas(gcf,['F' num2str(Function_name) '.jpg']) % ����ͼ��Ϊ jpg��ʽ
    %% ���浽excel
    filename = ['cec2018_Results_D' num2str(dim) '.xlsx']; % ������ļ�����
    sheet = 1; % ��һ��sheet
    func =['F' num2str(Function_name)]; % ��������
    xlswrite(filename, [Optimal_results{1, :}], sheet, 'B1:D1') % �㷨����
    xlswrite(filename, {func}, sheet, ['A' num2str(Function_name+1)]) % ��������
    xlswrite(filename, [Optimal_results{3, :}], sheet, ['B' num2str(Function_name+1)]) % ���Ž��

end