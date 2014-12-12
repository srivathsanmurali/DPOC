function [ J_opt, u_opt_ind ] = PolicyIteration( P, G )
%POLICYITERATION Value iteration
%   Solve a stochastic shortest path problem by policy iteration.
%
%   [J_opt, u_opt_ind] = PolicyIteration(P, G) computes the optimal cost and
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

MN = size(G,1);
L  = size(G,2);

u_opt_ind 	= ones(1,MN);
J_opt 		= zeros(1,MN);
done 	= 0;

PI_iter = 0;
max_iter = 1000;
quiet = true;
while done == 0
	% Main loop 
	% 	- ends when there is no change in the u_opt_ind;
	
	done = 0;
	% Left to one if no change occurs

	J_opt = getCost(u_opt_ind);
	% Step 1
	% solving the equations using VI method

	% Step 2
	% Improving the policy
	u_ind = u_opt_ind;
	for i=1:MN
		QBest = J_opt(i);
		for u=1:L
			Qsa = G(i,u) + P(i,:,u)*J_opt';
			if(Qsa > QBest)
				u_ind(i) = u;
				QBest = Qsa;
			end
		end
	end
	if(u_ind == u_opt_ind)
		done =1;
		if(~quiet)
			disp('Stopped because PI converged');
		end
	end
	u_opt_ind = u_ind;

	PI_iter = PI_iter + 1;
	
	if(PI_iter == max_iter)
		done = 1;
		if(~quiet)
			disp('Stoped because too many iterations');
		end
	end
	if(~quiet)
		display(PI_iter,'No of P iterations done');
	end
end
if(~quiet)
	display(PI_iter,'Total number of PI_iter');
end
% mx = max(J_opt);
% mn = min(J_opt);
% J_opt = (J_opt - mn) * (1000/(mx-mn));

	function [J] = getCost(u_ind)
		% getCost
		% 	- J (1XMN)
		% 		- Output matrix of Cost
		% 	- u_ind (1,MN)
		% 		- input matrix of policy

		Jk = zeros(1,MN);
		Jk1 = zeros(1,MN);

		T = 100;
		iter_done = 0;
		t = 1;
		while (~iter_done)
			for i=1:MN
				Jk1(i) = G(i,u_ind(i)) + P(i,:,u_ind(i))*Jk';
			end

			t = t +1;
			if(t == T)
				iter_done = 1;
				J = Jk;
			elseif((Jk -Jk1) == 0)
				iter_done = 1;
				J = Jk;
			else				
				Jk = Jk1;
			end
		end
	end
end

