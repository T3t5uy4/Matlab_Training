%% Simulated Annealing��SA
function [Best_score, Best_pos, curve] = SA(Mmax, l, u, dim, fobj)
    % ����:
    %        fobj = ��Ӧ�Ⱥ���
    %        x0 = ������Ⱥ
    %        l = ��Ⱥ�±߽�
    %        u = ��Ⱥ�ϱ߽�
    %        Mmax = ����¶�
    %        TolFun = �Ż��仯���̶�
    %
    %
    % ���:
    %        x0 = ����Ż������Ⱥ
    %        f0 = ����Ż������Ⱥ����Ӧ��ֵ
    TolFun = 10E-10; %ģ���˻����̶�
    x0 = (u - l) .* rand(1, dim) + l; %�����ʼ��ģ���˻�;
    f = fobj; %��Ӧ�Ⱥ���
    x = x0;
    fx = feval(f, x); %������Ӧ��ֵ
    f0 = fx;
    count = 1; %���ڼ�¼�������߱��
    %ģ���˻���Ҫ����
    for m = 1:Mmax
        T = m / Mmax; %�¶�
        mu = 10 ^ (T * 1000);
        %For each temperature we take 100 test points to simulate reach termal
        for k = 0:100
            dx = mu_inv(2 * rand(1, dim) - 1, mu) .* (u - l);
            x1 = x + dx;
            %�߽紦����ֹԽ��
            x1 = (x1 < l) .* l + (l <= x1) .* (x1 <= u) .* x1 + (u < x1) .* u;
            %���㵱ǰλ����Ӧ��ֵ����Ӧ��ֵƫ��
            fx1 = feval(f, x1); df = fx1 - fx;
            % ���df<0����ܸý⣬�������0 ������Metropolis׼������ж��Ƿ����
            if (df < 0 || rand < exp(-T * df / (abs(fx) + eps) / TolFun)) == 1
                x = x1; fx = fx1;
            end

            %�жϵ�ǰ���Ƿ���ţ����������.
            if fx1 < f0 == 1
                x0 = x1; f0 = fx1;
            end

        end

        curve(count) = f0;
        count = count + 1;
    end

    Best_pos = x0;
    Best_score = f0;
end

function x = mu_inv(y, mu)
    %ģ���˻������λ��ƫ��
    x = (((1 + mu) .^ abs(y) - 1) / mu) .* sign(y);
end
