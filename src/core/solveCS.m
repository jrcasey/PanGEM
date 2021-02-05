function [essentialRxns, results] = solveCS(model)


%% Setup problem for Mosek
S = model.S;
c = model.c;
beq = model.b;
lb = model.lb;
ub = model.ub;

% Size of the problem
[m,n] = size(S);

% Use the SIMPLEX algorithm rather than IP

% Feed options to Mosek
params = [];
params.relGap=0.4
linopts = getMILPParams_Mosek8(params)
linopts.MSK_IPAR_OPTIMIZER='MSK_OPTIMIZER_FREE_SIMPLEX';

[v,fval,exitflag,output] = linprog(-c,[],[],S,beq,lb,ub,[],linopts);

% Coding of lower/upper bounds as linear inequality with A_ulb v <= b_bb
A_ulb = [-diag(ones(n,1));diag(ones(n,1))];

b_bb = [-lb;ub];

% v and v2 should be equivalent!!! (and they are in terms of function value but not in terms of solution v/v2)
[v2,fval2,exitflag2,output2] = linprog(-c,A_ulb,b_bb,S,beq,[],[],[],linopts); % edited to new format (JRC)

% Check that fval and fval2 are close
if abs(fval2 - fval) < 1e-2
    result.LIQ = 1;
else result.LIQ = 0;
end

% Identify the column of S that relates to the objective function
indBiomass = find(c);

% Reformulate the problem by increasing the lower bound on indBiomass
% element
gVec = [0,0.001,0.002,0.003,0.004,0.005:0.005:-fval2];
nGs = length(gVec);

vMat = zeros(n,nGs);
fVec = zeros(1,nGs);

for i=1:nGs
    
    % Set the lower bound of
    lb(indBiomass) = gVec(i);
    b_bb = [-lb;ub];
    S_ex = [S,-S];
    A_pos = diag(-ones(2*n,1));
    b_pos = zeros(2*n,1);
    A_lb = [-diag(ones(n,1)),diag(ones(n,1))];
    A_ub = [diag(ones(n,1)),-diag(ones(n,1))];
    A = [A_pos;A_lb;A_ub];
    b = [b_pos;b_bb];
    c_lam = ones(2*n,1);
    [v_ex,fval_ex,exitflag_ex,output_ex] = linprog(c_lam,A,b,full(S_ex),beq,[],[],[],linopts);
    vL1 = v_ex(1:n)-v_ex(n+1:end);
    vMat(:,i) = vL1;
    fVec(i) = fval_ex;
end

essentialRxns = find(sum(abs(vMat(:,2:nGs))')>1e-12);
results.gVec = gVec;
results.vMat = vMat;
results.fVec = fVec;



