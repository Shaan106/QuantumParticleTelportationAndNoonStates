function ham = func_ham_maker_BHM(N,M)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ham_maker_simple.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description
%----------------------------------------------------------------------------
% Creates Hamiltonian matrix for atoms confinded to a Bose-Hubbard model
% trapping potential with arbitrary wells and atoms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs Matrices
%----------------------------------------------------------------------------
% basis - All possible configurations of the environment and system
% INTERACTION TERMS
% HU - Interaction strength is constant in all wells

% ENERGY OF LEVELS IN ENVIRONMENT 
% HE - Energy in each well can be different

% TUNNELING BETWEEN WELLS IN ENVIRONMENT
% HJ - Tunnelling strength between all wells is constant

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% N - Number of atoms in the environment
% M - Number of energy levels (modes) in the environment

% Creates the directory which the matrices will be saved to
data=sprintf('BHM_data_N=%i_M=%i',N,M);
mkdir(data);

% Creates a matrix holding all configurations of atoms
% [Note: I have changed so it works in freemat]

basis = create_basis_BHM(N,M);
L=length(basis);

%toc

%%%%%%%%%%%%%%%%%%%%%%
%%	Diagonal Terms	%%
%%%%%%%%%%%%%%%%%%%%%%
HU=sparse(zeros(L,L));      % Initiate interaction Hamiltonian
HE = zeros(L,L,M);          % Initiate energy of the well Hamiltonian
for m = 1:M
  HU = HU + diag(basis(:,m).*(basis(:,m)-1));
  HE(:,:,m) = diag(basis(:,m));
endfor

%%%%%%%%%%%%%%%%%%%%%%
% OFF DIAGONAL TERMS %
%%%%%%%%%%%%%%%%%%%%%%
HJ=sparse(zeros(L,L));              % Initiate tunnelling Hamiltonian

for j1 = 1:L
	vl = ones(L,1)*basis(j1,:);
  for m=1:M-1
    vr = basis;
    vr(:,m) = vr(:,m)+ones(L,1);
    vr(:,m+1) = vr(:,m+1)-ones(L,1);
    ind = find(sum(abs(vr-vl),2)==0);
    if length(ind)>0.1
      HJ(ind,j1) = HJ(ind,j1) - sqrt((basis(j1,m+1)+1)*basis(j1,m));
    end
  end
end
HJ = HJ+HJ';
%figure(1)
%image(abs(HJ*100))
%colorbar


%full(HJ + HU + HE(:,:,1))

ham = true;

if N>1
  data=sprintf('BHM_data_N=%i_M=%i/Hamiltonian.mat',N,M);
  save(data,'basis','HU','HE','HJ')
else
  data=sprintf('BHM_data_N=%i_M=%i/Hamiltonian.mat',N,M);
  save(data,'basis','HE','HJ')
endif



