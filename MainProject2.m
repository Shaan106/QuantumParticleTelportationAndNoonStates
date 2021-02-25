%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Investigates tunnelling between 3 wells

clear all
tic

global H
hbar = 1;

%%%%%%%%%%%%%%%%%%
%% Parameters	%%
%%%%%%%%%%%%%%%%%%
N=3;        % Number of atoms
M=3;        % Number of wells/modes - fixed to 2
E(1) = 0;   % Energy of left well
E(2) = 0.5;   % Energy of middle well - This will be modified!
E(3) = 0;    % Energy of right well
J  = 0.1;   % Tunnelling constant
U  = 0;    % Interaction constant (-ve means attractive)

%*********** Load Matrices ***********************
% Loads matrices and if they haven't been created it makes them first
data=sprintf('BHM_data_N=%i_M=%i/Hamiltonian',N,M);
dir=sprintf('BHM_data_N=%i_M=%i',N,M);

if exist(dir,'dir')
  load(data);
else
  ham = func_ham_maker_BHM(N,M);
  load(data);
end

L=length(basis);

%%%%%%%%%%%%%% Hamiltonian %%%%%%%%%%%%%%%%%%%%		
HET = zeros(L,L);	
for m=1:M
  HET = HET + E(m)*HE(:,:,m);
endfor

if N>1
  H = J*HJ + U*HU + HET;
else
  H = J*HJ + HET;
end

[vec, val] = eig(H);


%%%%%%%%%%% Time Evolution %%%%%%%%%%%%%%%%%%%%%%%%
%********** Initial state ************************
psi = zeros(L,1);
psi(1)=1;

%*********** Evolution Parameters *****************
tmax = 300;
dt = 0.1;
t = 0:dt:tmax;
Lt = length(t);
psiE = vec'*psi;

Exp_L = zeros(Lt,1);
Exp_R = zeros(Lt,1);
Exp_C = zeros(Lt,1);


for n=1:Lt
  psiEt = exp(-i*diag(val)*t(n)/hbar).*psiE;
  psit = vec*psiEt;
  Exp_L(n)=psit'*diag(basis(:,1))*psit;
  Exp_R(n)=psit'*diag(basis(:,3))*psit;
  Exp_C(n)=psit'*diag(basis(:,2))*psit;
  
endfor
figure(1)
plot(t,Exp_L,'r',t,Exp_C,'g',t,Exp_R,'b')
#figure(2)
#image([Exp_L,Exp_C,Exp_R]*50)







toc
