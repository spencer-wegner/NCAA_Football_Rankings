clear all; close all; clc;
%%

%% build (import/load) adjacency matrix
teams = {'Washington','Washington State','Stanford','California',...
    'Oregon State','Oregon','Colorado','USC','Utah','Arizona State',...
    'UCLA','Arizona'};
amat = xlsread('Amat.xlsx');
awaywinsmat = xlsread('AwayWinsMatrix.xlsx');
homelossmat = xlsread('HomeLossesMatrix.xlsx');
domlossmat = xlsread('DominatedLossesMatrix');
domwinsmat = xlsread('DominatingWinsMatrix');

A = amat(1:end-1,:);
B = awaywinsmat(1:end-1,:);
C = domwinsmat(1:end-1,:);
D = homelossmat(1:end-1,:);
E = domlossmat(1:end-1,:);

At = A';
b = ones(1,length(A));
TotalGames = (A + At)*b';
k = sum(TotalGames)/length(TotalGames);

% G = digraph(A,teams);
% plot(G)
% 
% [V,D] = eig(A);
% alpha_max = 1/max(max(real(D)));

alpha = 2*k/(k^2 - k);

%% Calc number of direct wins for team i
dw = zeros(1,length(A));
dl = zeros(1,length(A));
for i = 1:length(A)
    for j = 1:length(A)
        dw(i) = dw(i) + A(j,i);
        dl(i) = dl(i) + At(j,i);
    end
end
dw_mat = [1:length(A);dw];
dw_sorted = sortrows(dw_mat',2,'descend')';
ranked1 = {};
for i = 1:length(A)
   ranked1{i} = teams{dw_sorted(1,i)}; 
end

%% Calc number of indirect wins for team i
iw = zeros(1,length(A));
for i = 1:length(A)
    for j = 1:length(A)
        for k = 1:length(A)
            iw(i) = iw(i) + A(k,j)*A(j,i);
        end
    end
end

iw_mat = [1:length(A);iw];
iw_sorted = sortrows(iw_mat',2,'descend')';
ranked2 = {};
for i = 1:length(A)
   ranked2{i} = teams{iw_sorted(1,i)}; 
end

% ranking based on s = direct + indirect - losses
s2 = dw + iw - dl;
s2_mat = [1:length(A);s2];
s2_sorted = sortrows(s2_mat',2,'descend')';
ranked_std2 = {};
for i = 1:length(A)
   ranked_std2{i} = teams{s2_sorted(1,i)}; 
end

%% Ranking modifications

%% Team 1 wins

gamma = 0.1; % home field weight
sigma = 0.2; % dominated weight
c1 = dot(A,B); % weight mat 1
d1 = dot(A,C); % weight mat 2
c2 = dot(At,D); % weight mat 1
d2 = dot(At,E); % weight mat 2

w = zeros(1,length(A));
l = zeros(1,length(A));
for i = 1:length(A)
    for j = 1:length(A)
        w(i) = (dw(i) + gamma*c1(i) + sigma*d1(i)) + alpha*At(i,j)*dw(j);
        l(i) = (dl(i) + gamma*c2(i) + sigma*d2(i)) + alpha*A(i,j)*dl(j);
    end
end

s = w -l;

mod_mat = [1:length(A);s];
mod_sorted = sortrows(mod_mat',2,'descend')';
ranked_mod = {};
ranked_w{i} = {};
ranked_l{i} = {};
ranked_s{i} = {};
for i = 1:length(A)
   ranked_mod{i} = teams{mod_sorted(1,i)};
   ranked_w{i} = dw(mod_sorted(1,i));
   ranked_l{i} = dl(mod_sorted(1,i));
   ranked_s{i} = s(mod_sorted(1,i));
end

%% Ranked Teams
scores = [dw;iw;s];
ranked = [ranked1;ranked2;ranked_mod;ranked_w;ranked_l;ranked_s];
