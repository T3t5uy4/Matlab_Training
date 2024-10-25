%% ��������
function [ BestPosition,BestFitness ,fes ] = Con_OEDCanSolution_For_OMGSCA( Xi,gBest,dim,Q,F,fobj ,fes )
%UNTITLED3 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[ L,M ] = Get_QDOA( Q,F );
BestPosition=zeros(1,dim);
BestFitness=inf;
UB=max(Xi,gBest);
LB=min(Xi,gBest);
if dim>F
    a=2:dim-1;
    b=randperm(length(a));
    K=a(b(1:F-1));
    K=sort(K);
    for i=1:M
        for j=1:F
            if j==1
                Position(i,1:K(j))=LB(1:K(j))+(L(i,j)-1)*(UB(1:K(j))-LB(1:K(j)))./(Q-1);
            elseif j==F
                Position(i,K(j-1)+1:dim)=LB(K(j-1)+1:dim)+(L(i,j)-1)*(UB(K(j-1)+1:dim)-LB(K(j-1)+1:dim))./(Q-1);
            else
                Position(i,K(j-1)+1:K(j))=LB(K(j-1)+1:K(j))+(L(i,j)-1)*(UB(K(j-1)+1:K(j))-LB(K(j-1)+1:K(j)))./(Q-1);
            end
        end
        Fitness(1,i)=fobj( Position(i,:));
        fes = fes + 1;
    end
    for q=1:Q
        for i=1:M
            for j=1:F
                if L(i,j)==q
                    Z(i,j)=1;
                else
                    Z(i,j)=0;
                end
            end
        end
        S(q,:)=(Fitness*Z)./sum(Z);
    end
    [~,Index]=min(S);
    for j=1:F
        if j==1
            Position(M+1,1:K(j))=LB(1:K(j))+(Index(j)-1)*(UB(1:K(j))-LB(1:K(j)))./(Q-1);
        elseif j==F
            Position(M+1,K(j-1)+1:dim)=LB(K(j-1)+1:dim)+(Index(j)-1)*(UB(K(j-1)+1:dim)-LB(K(j-1)+1:dim))./(Q-1);
        else
            Position(M+1,K(j-1)+1:K(j))=LB(K(j-1)+1:K(j))+(Index(j)-1)*(UB(K(j-1)+1:K(j))-LB(K(j-1)+1:K(j)))./(Q-1);
        end
    end
    for i=1:M+1
        Fitness(1,i)=fobj( Position(i,:));
        fes = fes + 1;
        if Fitness(1,i)<BestFitness
            BestFitness=Fitness(1,i);
            BestPosition=Position(i,:);
        end
    end
else
	for i=1:M
        for j=1:dim
            Position(i,j)=LB(j)+(L(i,j)-1)*(UB(j)-LB(j))./(Q-1);
        end
         Fitness(1,i)=fobj( Position(i,:));
         fes = fes + 1;
    end
    for q=1:Q
        for i=1:M
            for j=1:dim
                if L(i,j)==q
                    Z(i,j)=1;
                else
                    Z(i,j)=0;
                end
            end
        end
        S(q,:)=(Fitness*Z)./sum(Z);
    end
    [~,Index]=min(S);
    for j=1:dim
        Position(M+1,j)=LB(j)+(Index(j)-1)*(UB(j)-LB(j))./(Q-1);
    end
    for i=1:M+1
        Fitness(1,i)=fobj( Position(i,:));
        fes = fes + 1;
        if Fitness(1,i)<BestFitness
            BestFitness=Fitness(1,i);
            BestPosition=Position(i,:);
        end
    end
end
end
