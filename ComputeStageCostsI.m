function G = ComputeStageCostsI( stateSpace, controlSpace, disturbanceSpace, mazeSize, walls, targetCell )
%COMPUTESTAGECOSTSI Compute stage costs.
% 	Compute the stage costs for all states in the state space for all
%   attainable control inputs.
%
%   G = ComputeStageCostsI(stateSpace, controlSpace, disturbanceSpace,
%   mazeSize, walls, targetCell) computes the stage costs for all states in
%   the state space for all attainable control inputs.
%
%   Input arguments:
%
%       stateSpace:
%           A (MN x 2) matrix, where the i-th row represents the i-th
%           element of the state space. Note that the state space also
%           contains the target cell, in order to simplify state indexing.
%
%       controlSpace:
%           A (L x 2) matrix, where the l-th row represents the l-th
%           element of the control space.
%
%       disturbanceSpace:
%           A (S x 3) matrix 'disturbanceSpace', where the first two
%           columns of each row represent an element of the disturbance
%           space, and the third column represents the corresponding
%           probability.
%
%       mazeSize:
%           A (1 x 2) matrix containing the width and the height of the
%           maze in number of cells.
%
%   	walls:
%          	A (2 x 2K) matrix containing the K wall segments, where the start
%        	and end point of the k-th segment are stored in column 2k-1
%         	and 2k, respectively.
%
%    	targetCell:
%          	A (2 x 1) matrix describing the position of the target cell in
%         	the maze.
%
%   Output arguments:
%
%       G:
%           A (MN x L) matrix containing the stage costs of all states in
%           the state space for all attainable control inputs. The entry
%           G(i, l) represents the cost if we are in state i and apply
%           control input l.

% put your code here

% Algo
%   if possible path
%       if takes two steps g = 2
%       if takes one steps g = 1
%       if takes zero steps g = 0
%% Constants
MN  = size(stateSpace, 1);
M   = mazeSize(1); %Vertical
N   = mazeSize(2); %Horizontal
L   = size(controlSpace,1);

MoveMatrix = getMoveMatrix();
controlCost = sum(abs(controlSpace),2);
G = zeros(MN,L);
% Need to add the changes caused due to walls.
for i=1:MN
    [controlSpaceNew] = getPossibleMoves(i);
    for l=1:L
        if(controlSpaceNew(l) == 0)
            continue;
        end
        
        G(i,l) = controlCost(l);
    end
end
%% Get Possible Moves
    function [L_new] = getPossibleMoves(i)
        L_new = zeros(L,1);
        x = i/M;
        y = mod(i,M);
        
        for ll=1:L
            u = controlSpace(ll,2) + (controlSpace(ll,1) * M);
            x = i + u;
            if(x<=MN && x>=1 && MoveMatrix(i,x) == 0)            
                L_new(ll) = 1;
            end
        end
		diagn = [M+1,M-1,-M+1,-M-1];
		diagnL = [12,10,4,2];
		dnodes = [M,1;M,-1;-M,1;-M,-1];
		for dg=1:4
			d = i + diagn(dg);
			n1 = i + dnodes(dg,1);
			n2 = i + dnodes(dg,2);

			if(d<1 || d>MN || n1<1 || n1>MN || n2<1 || n2>MN)
				L_new(diagnL(dg)) = 0;
			else
				if(MoveMatrix(i,n1) == 0 && MoveMatrix(i,n2) == 0 ...
					&& MoveMatrix(d,n1) == 0 && MoveMatrix(d,n2) == 0 )
					L_new(diagnL(dg)) = 1;
				else
					L_new(diagnL(dg)) = 0;				
				end
			end

		end		

        if(y==0)
            L_new(4) = 0;
            L_new(12) = 0;
        end
        if(y==1)
            L_new(10) = 0;
            L_new(2) = 0; 
        end
        if(x==1)
            L_new(2) = 0;
            L_new(4) = 0;
        end
        if(x==M-1)
            L_new(10) = 0;
            L_new(12) = 0;
        end
    end

    function [MM] = getMoveMatrix()
        MM = zeros(MN,MN);
        %borders
        bdrs = [1];%,M+1,-M+1]
        for i=M:M:MN
            for j=1:length(bdrs)
                x = i + bdrs(j);
                if(x>=1 && x<=MN)
                    MM(i,x) = 1;
                end
            end
        end
        for i=1:M:MN
            for j=1:length(bdrs)
                x = i - bdrs(j);
                if(x>=1 && x<=MN)
                    MM(i,x) = 1;
                end
            end
        end
        WL = size(walls,2)/2;
        for wl=1:WL
            from_x = walls(1,2*wl-1);
            from_y = walls(2,2*wl-1);
            to_x = walls(1,2*wl);
            to_y = walls(2,2*wl);

            if(from_x == to_x) %vertical wall
                from    = max(to_y,from_y) + ((from_x-1)*M);
                to      = max(to_y,from_y) + ((from_x)*M);
                while (to<=MN)
                    while (from>=1)
                        MM(from,to) = 1;
                        MM(to,from) = 1;
                        from = from - M;
                    end
                    from    = max(to_y,from_y) + ((from_x-1)*M);
                    to = to + M;
                end
            end
            if(from_y == to_y) %horizontal wall
                from    = (from_y)  + (min(to_x,from_x)*M);
                to      = (from_y+1)+ (min(to_x,from_x)*M);
                to_stop = (ceil(to/M) * M);
                from_stop = ( floor(from/M) * M) + 1;
                while (to<=to_stop)
                    while (from>=from_stop)
                        MM(from,to) = 1;
                        MM(to,from) = 1;
                        from = from - 1;
                    end
                    from    = (from_y)  + (min(to_x,from_x)*M);
                    to = to + 1;
                end
            end
            
        end 
    end
end

