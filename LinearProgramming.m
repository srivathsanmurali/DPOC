function [ J_opt, u_opt_ind ] = LinearProgramming( P, G )
%LINEARPROGRAMMING Value iteration
%   Solve a stochastic shortest path problem by linear programming.
%
%   [J_opt, u_opt_ind] = LinearProgramming(P, G) computes the optimal cost
%   and the optimal control input for each state of the state space.
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

MN = size(G,1);
L  = size(G,2);

u_opt_ind 	= ones(1,MN);
J_opt 		= zeros(1,MN);

A = zeros(MN*L,MN);
B = zeros(MN*L,1);

for i=1:MN
	for u=1:L
		for j=1:MN
			a = ((i-1) * MN) + u;	
			if(i~=j)
				A(a,j) = P(i,j,u);
			else
				A(a,j) = P(i,j,u) - 1;
			end
			B(a) = G(i,u);
		end
	end
end

f = ones(MN,1) * -1;
lb = zeros(MN,1);
J_opt = linprog(f,A,B,[],[],lb);


J_opt = J_opt';

for i=1:MN
	QBest = 0;
	for u=1:L
		Qsa = G(i,u) + P(i,:,u)*J_opt';
		if(Qsa > QBest)
			u_opt_ind(i) = u;
			QBest = Qsa;
		end
	end
end

end

