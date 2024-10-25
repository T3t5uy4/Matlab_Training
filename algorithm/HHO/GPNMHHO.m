% Developed in MATLAB R2013b
%
%  Author and programmer: Ali Asghar Heidari, University of Tehran, 09-28-2017
%         e-Mail: as_heidari@ut.ac.ir, aliasghar68@gmail.com

% Harris's hawk optimizer
% The HHO Algorithm: In this algorithm, Harris' hawks try to
% catch the rabbit.

% T: maximum iterations, N: populatoin size, CNVG: Convergence curve
function [N,Rabbit_Location, CNVG]=GPNMHHO(N,MaxFES,lb,ub,dim,fobj)

disp('HHO is now tackling your problem')
tic
% initialize the location and Energy of the rabbit
Rabbit_Location=zeros(1,dim);
Rabbit_Energy=inf;

%Initialize the locations of Harris' hawks
X=initialization(N,dim,ub,lb);
N = size(X,1);
r = zeros(size(X));
Boundary_size = size(ub,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%佳点集进行初始化%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:N
    for j = 1:size(X,2)
        r(i,j) = i * exp(j) - floor(i * exp(j));
    end
end
for i = 1:N
    for j = 1:size(X,2)
        if Boundary_size > 1
             X(i,j) = lb(j) + r(i,j) * (ub(j) - lb(j));
        else
            X(i,j) = lb + r(i,j) * (ub - lb);
        end
    end   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CNVG=[];
FES=0; % FES counter

fitcount=0; % Loop counter

while FES<MaxFES
    for i=1:size(X,1)
        % Check boundries
        FU=X(i,:)>ub;FL=X(i,:)<lb;X(i,:)=(X(i,:).*(~(FU+FL)))+ub.*FU+lb.*FL;
        % fitness of locations
        fitness=fobj(X(i,:));
        FES=FES+1;
        % Update the location of Rabbit
        if fitness<Rabbit_Energy
            Rabbit_Energy=fitness;
            Rabbit_Location=X(i,:);
        end
        fitcount = fitcount + 1;
        CNVG(1,fitcount) = Rabbit_Energy;
    end
    
    c=2*(1-(FES/MaxFES)); % factor to show the decreaing energy of rabbit
    % Update the location of Harris' hawks
    for i=1:size(X,1)
        
        Escaping_Energy=c*(2*rand()-1);  % escaping energy of rabbit
        
        if abs(Escaping_Energy)>=1
            %% Exploration:
            % Harris' hawks search for the rabbit, observing the area to find a rabbit
            % Harris' hawks perch randomly based on 2 strategy:
            
            r=rand();
            rand_Hawk_index = floor(N*rand()+1);
            X_rand = X(rand_Hawk_index, :);
            if r<0.5
                % perch based on other family members
                X(i,:)=X_rand-rand()*abs(X_rand-2*rand()*X(i,:));
            elseif r>0.5
                % perch on a random tall tree (random site inside group's home range)
                X(i,:)=(Rabbit_Location(1,:)-mean(X))-rand()*((ub-lb)*rand+lb);
            end
            
        elseif abs(Escaping_Energy)<1
            %% Exploitation:
            % Attacking the rabbit using 4 strategies regarding the behavior of the rabbit
            
            %% phase 1: surprise pounce (seven kills): hawks go to kill when rabbit is surprised
            % surprise pounce (seven kills): multiple, short rapid dives by different hawks
            
            r=rand(); % probablity of each event % r>0.5 : rabbit cannot escape -- r<0.5 : rabbit can escape
            
            if r>0.5 && abs(Escaping_Energy)<0.5 % Hard besiege % When rabbit cannot or do not try to escape
                X(i,:)=(Rabbit_Location)-Escaping_Energy*abs(Rabbit_Location-X(i,:));
            end
            
            if r>0.5 && abs(Escaping_Energy)>0.5  % Soft besiege % When rabbit cannot or do not try to escape
                Jump_strength=2*(1-rand()); % random jump strength of the rabbit
                X(i,:)=(Rabbit_Location-X(i,:))-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:));
            end
            
            %% phase 2: performing team rapid dives (leapfrog movements) when prey escapes
            if r<0.5 && abs(Escaping_Energy)>0.5, % Soft besiege % rabbit try to escape by many zigzag deceptive motions
                
                Jump_strength=2*(1-rand());
                X1=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:));
                fx1 = fobj(X1);
                fxi = fobj(X(i,:));
                if fx1 < Rabbit_Energy
                    Rabbit_Energy=fx1;
                    Rabbit_Location=X1;
                end
                FES = FES + 1;
                fitcount = fitcount + 1;
                CNVG(1,fitcount) = Rabbit_Energy;
                if fxi < Rabbit_Energy
                    Rabbit_Energy=fxi;
                    Rabbit_Location=X(i,:);
                end
                FES = FES + 1;
                fitcount = fitcount + 1;
                CNVG(1,fitcount) = Rabbit_Energy;
                if fx1<fxi % improved move?
                    X(i,:)=X1;
                else % hawks perform levy-based short rapid dives around the rabbit
                    X2=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:))+rand(1,dim).*Levy(dim);
                    fx2 = fobj(X2);
                    if fx2 < Rabbit_Energy
                        Rabbit_Energy=fx2;
                        Rabbit_Location=X2;
                    end
                    FES = FES + 1;
                    fitcount = fitcount + 1;
                    CNVG(1,fitcount) = Rabbit_Energy;
                    if fx2<fxi, % improved move?
                        X(i,:)=X2;
                    end
                end
            end
            
            if r<0.5 && abs(Escaping_Energy)<0.5, % Hard besiege % rabbit try to escape by many zigzag deceptive motions
                % hawks try to decrease their average location with the rabbit
                Jump_strength=2*(1-rand());
                X1=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-mean(X));
                fx1 = fobj(X(i,:));
                if fx1 < Rabbit_Energy
                    Rabbit_Energy=fx1;
                    Rabbit_Location=X1;
                end
                FES = FES + 1;
                fitcount = fitcount + 1;
                CNVG(1,fitcount) = Rabbit_Energy;
                fxi = fobj(X(i,:));
                if fxi < Rabbit_Energy
                    Rabbit_Energy=fxi;
                    Rabbit_Location=X(i,:);
                end
                FES = FES + 1;
                fitcount = fitcount + 1;
                CNVG(1,fitcount) = Rabbit_Energy;
%                 FES=FES+2;
                if fx1<fxi % improved move?
                    X(i,:)=X1;
                else % Perform levy-based short rapid dives around the rabbit
                    X2=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-mean(X))+rand(1,dim).*Levy(dim);
                    fx2 = fobj(X2);
                    if fx2 < Rabbit_Energy
                        Rabbit_Energy=fx2;
                        Rabbit_Location=X2;
                    end
                    FES = FES + 1;
                    fitcount = fitcount + 1;
                    CNVG(1,fitcount) = Rabbit_Energy;
                    if fx2<fxi, % improved move?
                        X(i,:)=X2;
                    end
                end
            end
            %%
        end
    end
    options = optimset('MaxFunEvals', floor(MaxFES * 0.1));
    [x, fval, ~, output]  = fminsearchbnd(fobj,Rabbit_Location,lb,ub,options);
    if fval < Rabbit_Energy
        Rabbit_Energy = fval;
        Rabbit_Location = x;
    end
   FES = FES + output.funcCount;
   CNVG(:,fitcount + 1:fitcount + output.funcCount) = Rabbit_Energy;
   fitcount = fitcount + output.funcCount;
end
CNVG = CNVG(:,1:MaxFES);
toc
end

function o=Levy(d)
beta=1.5;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;v=randn(1,d);step=u./abs(v).^(1/beta);
o=step;
end