function [ J_opt, u_opt_ind ] = ValueIteration( P, G )
%VALUEITERATION Value iteration
%   Solve a stochastic shortest path problem by value iteration.
%
%   [J_opt, u_opt_ind] = ValueIteration(P, G) computes the optimal cost and
%   the optimal control input for each state of the state space.
%
%   Input arguments:
%
%       P:
%           A (MN x MN x L) matrix containing the transition probabilities
%           between all states in the state space for all attainable
%           control inputs. The entry P(i, j, l) represents the transition
%           probability from state i to state j if control input l is
%           applied.
%
%       G:
%           A (MN x L) matrix containing the stage costs of all states in
%           the state space for all attainable control inputs. The entry
%           G(i, l) represents the cost if we are in state i and apply
%           control input l.
%
%   Output arguments:
%
%       J_opt:
%       	A (1 x MN) matrix containing the optimal cost-to-go for each
%       	element of the state space.
%
%       u_opt_ind:
%       	A (1 x MN) matrix containing the indices of the optimal control
%       	inputs for each element of the state space.

% put your code here
%constants
MN = size(G,1);
L  = size(G,2);

J_opt = max(G,[],2);
u_opt_ind = ones(1,MN);
Jk = max(G,[],2)';
Jk1 = zeros(1,MN);
is_done = 0;

iter = 0;
while (is_done == 0)
	iter = iter + 1
	u_old = u_opt_ind;
	for i=1:MN

		tempJ = 0;
		Jk1(i) = 0;
		for u=1:L
			tempJ = G(i,u);
			for j=1:MN
				tempJ = tempJ + P(i,j,u)*Jk(j);
			end
			% display([iter,i,u,tempJ])
			if(i==8 && P(8,9,u) ~= 0)
				display([8,9,u,P(8,9,u),Jk(9),tempJ])
			end
			if(tempJ>Jk1(i))
				Jk1(i) = tempJ;
				u_opt_ind(i) = u;
				%display(1)
			end
		end
	end
	if((Jk - Jk1) < 0.00001)
		J_opt = Jk1;
		is_done = 1;
	end
	Jk = Jk1;
end
display(iter, 'total iterations')
end

