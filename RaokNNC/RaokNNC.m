% ------------------------------------------------------------------------%
% Rao with k Nearest Neighbor Comparison Method
% <RaokNNC>
% Hoang-Anh Pham, 2024
% Department of Structural Mechanics, 
% Hanoi University of Civil Engineering
% Email: anhph2@huce.edu.com
% ------------------------------------------------------------------------%
% Reference:
% Pham, H. A., Dang, V. H., Vu, T. C., & Nguyen, B. D. (2024). 
% An Efficient k-NN-based Rao Optimization Method for 
% Optimal Discrete Sizing of Truss Structures. 
% Applied Soft Computing, 111373.
% ------------------------------------------------------------------------%
function [xopt,fopt,exitflag,output,X,Fitness,history,FE,DI,S] = RaokNNC(fname,fcons,nvars,LB,UB,DX,option,type,nnc,kv)

    % Setting algorithm parameters
    T_max = option(1); 
    Tol = option(2);
    NP = option(3);

    % Initialize the population/solutions
    [row,col] = size(DX);
    rnd = rand(NP,nvars);
    %rnd = lhsdesign(NP,nvars);
    lb = repmat(LB,NP,1);
    ub = repmat(UB,NP,1);
    X = lb + rnd.*(ub-lb);
    for i=1:NP
        for j=1:row, X(i,DX(j,1)) = dismap(X(i,DX(j,1)),DX(j,2:col)); end     
        % Evaluate fitness and constraint violation
        Fitness(i) = feval(fname,X(i,:)); 
        Const(i) = constraint(fcons,X(i,:),1e-4); 
    end    
    
    % Sort the population in ascending order of fitness       
    [~,ID] = consort(Fitness,Const,20,0); 
    xbest = X(ID(1),:); 
    fbest = Fitness(ID(1)); 
    cbest = Const(ID(1));    
    fmean = mean(Fitness);
    en=0;
    
    % Calculate diversity index of initial population
    DI_0 = d_index(X,LB,UB); 
    DI_t = DI_0;
    DI = DI_t; 
    
    % Initializing 
    num_FE = NP; num_succ = 0; num_fail = 0; 
    num_skip = 0; 
    S = num_skip; 
    
    N = [0,0]; history = []; 
    FE = num_FE;
    
%    fprintf('\n						 Best\t\t    Mean\n')
%    fprintf('Generation\tNFE\t\t f(x)\n')
    
    % Start the iterations
    iter=0; dev = Tol+1;
    while (iter < T_max) && (dev > Tol),
       
%       fprintf('\t%i\t\t%i\t\t\t %f\t %f\t\n',...
%            iter,num_FE,fbest,fmean)

       iter=iter+1;
       
       inum_FE=0; inum_skip = 0; 

       % Sort the population in ascending order of fitness       
       [~,ID] = consort(Fitness,Const,20,0); 
       
       xbest = X(ID(1),:);
       fbest = Fitness(ID(1)); 
       cbest = Const(ID(1)); 
       xold = xbest; 
       fold = fbest; 
       xworst = X(ID(end),:);

       X = X(ID(1:end),:); 
       Fitness = Fitness(ID(1:end));
       Const = Const(ID(1:end));
              
       % Loop over all solutions
       for k=1:NP
           % Take random vector from current solutions
           JK = randperm(round(NP)); JK(JK==k)=[]; 
           r1 = JK(1); 

           % Set direction for mutation
           ecom = ebetter(Fitness(k),Const(k),Fitness(r1),Const(r1),en); 
           if ecom==1, d=1; else d=-1; end;   
           
           %% Select method
           switch type
                % Rao-1
                case '1', 
                x_new = X(k,:) + rand(1,nvars).*(xold-xworst);
                % Rao-2
                case '2', 
                x_new = X(k,:) + rand(1,nvars).*(xold-xworst) + ...
                               d*rand(1,nvars).*(X(k,:)-X(r1,:));
                case 'c', 
                p =(k-1)/(NP-1);
                x_new = X(k,:) + rand(1,nvars).*(xold-xworst) + ...
                               p*d*rand(1,nvars).*(X(k,:)-X(r1,:));                           
           end          

           %% Check bound constraints
           x_new = boundConst(x_new,X(k,:),LB,UB);
           for j=1:row, x_new(DX(j,1)) = dismap(x_new(DX(j,1)),DX(j,2:col)); end     
           
           % Apply the k-NNC
           eval=1; 
           switch nnc  
                case 'knn' 
                eval=kNNC(k,x_new,X,Fitness,Const,en,kv);      
           end
           if eval==0,
              inum_skip=inum_skip+1;
           end
           
           % Check to skip usefuless evaluation
           if eval>0, eval=0;
              fnew = feval(fname,x_new); 
              
              inum_FE = inum_FE + 1; 

              if Const(k)>en, eval=1;
              else if fnew<Fitness(k), eval=1; end
              end
              history = [history,fbest];
           end
         
           if eval>0,           
                cnew = constraint(fcons,x_new,1e-4);
                
                % Select new solution as member if better
                ecom  = ebetter(fnew,cnew,Fitness(k),Const(k),en); 
                if (ecom == 1), 
                    X(k,:) = x_new; 
                    Fitness(k) = fnew; 
                    Const(k) = cnew;
                    
                    num_succ = num_succ + 1; 
                else
                    num_fail = num_fail + 1;
                end;                
           end      
           
           % Update the current global best
           ecom = ebetter(Fitness(k),Const(k),fbest,cbest,0);
           if ecom==1, xbest = X(k,:); fbest = Fitness(k); cbest = Const(k); end                          
           if fbest<fold, N=[iter,num_FE]; fold=fbest; end
           
       end % End for
       
       DI_t = d_index(X,LB,UB); DI=[DI,DI_t];

       num_FE = num_FE + inum_FE; 
       num_skip = num_skip + inum_skip;
       
       S = [S,num_skip/(iter*NP)];
       FE = [FE,num_FE];
       
       fmean = mean(Fitness); 
       dev = abs(fmean/fbest-1);
      
    end % End while
    
    % Return the optimal solution
    xopt = xbest; fopt = fbest;
    
    if iter>=T_max, 
        exitflag = 0;   massage = 'generations exceed max. generation';
    else exitflag = 1;  massage = 'diversity is lower than Tol';
    end
    
    output = struct('generations',iter,...
                    'skiprate',num_skip/(num_FE+num_skip)*100,...
                    'succskip',0,...                    
                    'funccount',num_FE,...
                    'constcount',num_succ+num_fail,...
                    'success',num_succ,...
                    'fail',num_fail,...
                    'minfunc',N,...
                    'maxconstraint',cbest,...
                    'diversity',d_index(X,LB,UB),...
                    'message',massage);    
end % Main function



