function P = ComputeTransitionProbabilitiesII( stateSpace, controlSpace, disturbanceSpace, mazeSize, walls, targetCell, holes, resetCell )
%COMPUTETRANSITIONPROBABILITIESII Compute transition probabilities.
% 	Compute the transition probabilities between all states in the state
%   space for all attainable control inputs.
%
%   P = ComputeTransitionProbabilitiesII(stateSpace, controlSpace,
%   disturbanceSpace, mazeSize, walls, targetCell, holes, resetCell)
%   computes the transition probabilities between all states in the state
%   space for all attainable control inputs.
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
%    	holes:
%         	A (2 x H) matrix containg the H holes of the maze. Each column
%         	represents the position of a hole.
%
%   	resetCell:
%         	A (2 x 1) matrix describing the position of the reset cell in
%           the maze.
%
%   Output arguments:
%
%       P:
%           A (MN x MN x L) matrix containing the transition probabilities
%           between all states in the state space for all attainable
%           control inputs. The entry P(i, j, l) represents the transition
%           probability from state i to state j if control input l is
%           applied.

% put your code here
%% Constants

MN  = size(stateSpace, 1);
M   = mazeSize(1); %Vertical
N   = mazeSize(2); %Horizontal
S   = size(disturbanceSpace,1);
L   = size(controlSpace,1);
p_pc= 1/L; %probability of each control space
p_pc_norm = 0;
P = [];
target 		= targetCell(2) + ((targetCell(1)-1)*M);
rst 		= resetCell(2) + ((resetCell(1)-1) *M);
Holes = holes(2,:) + ((holes(1,:) -1) * M);

Walls = getWalls();
% Wall matrix
%     - 1 when allowed to move to the node
%     - 0 when not allowed to move to the node
MoveMatrix = getMoveMatrix();
% Move matrix
%     - 1 when allowed to move to the node
%     - 0 when not allowed

% display(MoveMatrix(9,:),'MM(9,:)')
P = zeros(MN,MN,L);

% Need to add the changes caused due to walls.
for i=1:MN
    [controlSpaceNew,p_pc1] = getPossibleMoves(i);
    % Possible Moves
    %     1 - If allowed
    %     0 - Not allowed
    for l=1:L
        if(controlSpaceNew(l) == 0)
            continue;
        end
        % control movement % disturbance movement
        x_c = controlSpace(l,1);
        y_c = controlSpace(l,2);
        u = y_c + (x_c*M);
        %dynamics of the system
        x = i + u;
        if(x<=MN && x>= 1 && MoveMatrix(i,x) == 1)

        	got_hole = find(Holes == x);
        	if(isempty(got_hole))
	            P(i,x,l) = p_pc1(l);
	            
	            for s=1:S
	                x_d = disturbanceSpace(s,1);
	                y_d = disturbanceSpace(s,2);
	                d = y_c + (x_c*M);
	                xd = x+d;
	                if(xd <= MN && xd >= 1 && MoveMatrix(x,xd) == 1)
	                    got_hole = find(Holes == xd);
	                    if(isempty(got_hole))
	                    	P(i,xd,l) = p_pc1(l) * disturbanceSpace(s,3);
	                    else
	                    	P(i,rst,l) = p_pc1(l) * disturbanceSpace(s,3);
	                    end
	                end
	            end
	        else
	        	P(i,x,l) = 0;
	        	P(i,rst,l) = p_pc1(l);
	        end
        end 
    end
end
P(target,:,:) = 0;
P(target,target,:) = 1;

%% Get Possible Moves
    function [L_new,p_pc_m] = getPossibleMoves(i)
        L_new = zeros(L,1);
        p_pc_m = zeros(L,1);
        
        for ll=1:L
            u = controlSpace(ll,2) + (controlSpace(ll,1) * M);
            x = i + u;
            if(x<=MN && x>=1 && MoveMatrix(i,x) == 1)            
                L_new(ll) = 1;
                p_pc_m(ll) = p_pc;
            else
                p_pc_m(ll) = 0;
                % p_pc_m(7) = p_pc_m(ll) + p_pc;
            end                
        end

        % p = 1/length(find(L_new == 1));
        % p_pc_m = p_pc_m  * p;
        
    end

    function [MM] = getMoveMatrix()
        MM = zeros(MN,MN);
        for i = 1:MN
            for u= 1:L
                x_c = controlSpace(u,1);
                y_c = controlSpace(u,2);
                d_c = x_c * y_c;

                x_i = i/M;
                y_i = mod(i,M);

                u = y_c + (x_c*M);
                x = i + u;
                if(x>=1 && x<=MN)
                    if(d_c == 0 && Walls(i,x) == 1)

                        MM(i,x) = 1;
                    else
                        if(y_i == 0 && y_c == 1)
                            MM(i,x) = 0;
                        elseif(y_i == 1 && y_c == -1)
                            MM(i,x) = 0;
                        elseif(x_i == 1 && x_c == -1)
                            MM(i,x) = 0;
                        elseif(x_i == N-1 && x_c == 1)
                            MM(i,x) = 0;
                        else
                            n1 = i + y_c;
                            n2 = i + (x_c*M);
                            if(Walls(i,n1)== 1 && Walls(i,n2) == 1 ...
                                && Walls(x,n1) == 1 && Walls(x,n2) == 1)

                                MM(i,x) = 1;
                            else
                                MM(i,x) = 0;
                            end
                        end
                    end

                end
            end
        end
    end




    function [Walls] = getWalls()
        Walls = ones(MN,MN);
        %borders
        bdrs = [1];
        for i=M:M:MN
            for j=1:length(bdrs)
                x = i + bdrs(j);
                if(x>=1 && x<=MN)
                    
                    Walls(i,x) = 0;
                    Walls(i-1,x) = 0;
                    Walls(i,x+1) = 0;
                end

            end
        end
        for i=1:M:MN
            for j=1:length(bdrs)
                x = i - bdrs(j);
                if(x>=1 && x<=MN)
                    Walls(i,x) = 0;
                    Walls(i+1,x) = 0;
                    Walls(i,x-1) = 0;
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

                        Walls(from,to) = 0;
                        Walls(to,from) = 0;
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
                        Walls(from,to) = 0;
                        Walls(to,from) = 0;
                        from = from - 1;
                    end
                    from    = (from_y)  + (min(to_x,from_x)*M);
                    to = to + 1;
                end
            end
            
        end 

    end
end

