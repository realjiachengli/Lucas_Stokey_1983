addpath('/Users/jasonli/Course Space/2022 Spring/Numerical Methods/code directory/libm14022022')


tic;
% define parameters
param.nS = 4;       % num of states
param.beta = 0.98;  % discount factor
param.gamma = 0.7;  % weight on logL

% param.Trans = [
%     0.5,0.5;
%     0.5,0.5;
% ];

% param.Trans = [
%     0,0.5,0.5;
%     0.6,0.4,0;
%     0,0.2,0.8;
% ];

param.Trans = [
    0,0.5,0.2,0.3;
    0.6,0.2,0.1,0.1;
    0.2,0.3,0.2,0.3;
    0.2,0.7,0,0.1;
];

% param.Trans = [
%     0,0.5,0.2,0,0.2,0.1;
%     0.6,0.2,0.1,0.05,0,0.05;
%     0.2,0.3,0,0,0.2,0.3;
%     0.2,0.3,0.3,0.1,0,0.1;
%     0.2,0.3,0,0,0.2,0.3;
%     0.4,0.4,0.1,0.05,0,0.05;
% ];


% param.gfunc = [0.1,0.2];
param.gfunc = [0.1,0.2,0.4,0.2];


param.x_grid = linspace(-5.,5.,200);
param.REDUCE = 1;

param.tol = 1e-5;
param.MaxIter = 500;
param.opts = optimset('Display','off');

% extrapolation penalty factor
param.lambda_upp = 0.8;
param.lambda_low = 1.2;     
x0.s0 = 1;
x0.b0 = 0.5;

outcome = time0_problem(x0,param);
toc;