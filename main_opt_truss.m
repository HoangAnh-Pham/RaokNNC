% ------------------------------------------------------------------------%
% MAIN PROGRAM
% 10-bar truss sizing
% Hoang-Anh Pham, 2024
% Department of Structural Mechanics, 
% Hanoi University of Civil Engineering
% Email: anhph2@huce.edu.com
% ------------------------------------------------------------------------%
addpath('fem');
addpath('RaokNNC');

clc; close all; clear all; 
global nvars XB XV

%% Setting optimization parameters
Ng = 1000;      % No. iterations
Tol = 1e-6;     % Min. relative error
NP = 25;        % Population size
NoR = 1;        % No. runs
para = [Ng, Tol, NP];

%% Optimization problem
problem = {{'3-bar',@truss3obj,@truss3cons,@truss3data};        % 1
           {'10-bar',@truss_obj,@truss10cons,@truss10data}};    % 2

pb = 2;  % problem ID
tp = 2;  % 1- continuous; 2- discrete

%% Get structural data
truss_name = problem{pb}{1};
fname = problem{pb}{2};
fcons = problem{pb}{3};
data = problem{pb}{4};

feval(data);

switch tp
    case 1
    B = XB; DX = [];
    case 2
    B = XV;
    for i=1:nvars, DX(i,:) = [i,XV]; end
end

LB = min(B)*ones(1,nvars); % Lower bound
UB = max(B)*ones(1,nvars); % Upper bound

%% Optimization method
ID = [6,9]; % List of method ID
method={{'DE',          @dDEmRao,'rnd','',''};          % 1        
        {'Rao1',        @dDEmRao,'rao','',''};          % 2
        {'DE-Rao1',     @dDEmRao,'hb1','',''};          % 3
        {'DE-Rao1-DiC', @dDEmRao,'hb1','','dic'};       % 4
        {'dDEmRao-DiC', @dDEmRao,'hb2','d','dic'};      % 5
        {'Rao-1',       @RaokNNC,'1','',3};             % 6     
        {'Rao-2',       @RaokNNC,'2','',3};             % 7     
        {'cRao',        @RaokNNC,'c','',3};             % 8     
        {'Rao1-kNNC',   @RaokNNC,'1','knn',3};          % 9     
        {'Rao2-kNNC',   @RaokNNC,'2','knn',3};          % 10     
        {'cRao-kNNC',   @RaokNNC,'c','knn',3};          % 11     
        };
method = method(ID); 
mt = 1:length(method);

%% Run optimization
disp(['Problem: ',truss_name]);
for i=1:length(mt)
    tic;
    disp(['Method: ',method{mt(i)}{1}]);
    algorithm = method{mt(i)}{2};
    option = {method{mt(i)}{3:end}};
    varin = {algorithm,fname,fcons,nvars,LB,UB,DX,para,option{:}};
    
    for t=1:NoR
        disp(['Run: ',int2str(t)]);
        [xopt,fopt,exitflag,out,X,scores,V,FE,DI,S] = feval(varin{:});
        
        % show results
        disp(['Optimized weight:',num2str(fopt)]);
        disp(['CV:',num2str(max(feval(fcons,xopt)))]);
        disp(out);
        
        figure; hold all; box on;
        title([truss_name,', ',method{mt(i)}{1},' ','run ',num2str(t)]);
        plot(FE,'-b','LineWidth',1.5);  
        
        ylabel('No. function evaluations');
        xlabel('Iterations');
        hold off;
        
        figure; box on;
        plot(V,'LineWidth',1.5); 
        title([truss_name,', ',method{mt(i)}{1},' ','run ',num2str(t)]);
        xlabel('FEs'); ylabel('Weight');
        hold off;
        
        figure; box on;
        semilogy(DI,'LineWidth',1.5); 
        title([truss_name,', ',method{mt(i)}{1},' ','run ',num2str(t)]);
        xlabel('Iterations'); ylabel('Diversity index'); box on;
        hold off;
        
        % save results
        save([truss_name,'-',method{mt(i)}{1},'_T',num2str(t),'.mat']);
    end
    toc;
end


