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
target = targetCell(2) + ((targetCell(1)-1)*M);

Walls = getWalls();
% Wall matrix
%     - 1 when allowed to move to the node
%     - 0 when not allowed to move to the node
MoveMatrix = getMoveMatrix();
% Move matrix
%     - 1 when allowed to move to the node
%     - 0 when not allowed

controlCost = sum(abs(controlSpace),2);
G = ones(MN,L) * -10;
% Need to add the changes caused due to walls.
for i=1:MN
    [controlSpaceNew] = getPossibleMoves(i);
    for l=1:L
        x_c = controlSpace(l,1);
        y_c = controlSpace(l,2);
        c = y_c + (x_c * M);
        
        if(controlSpaceNew(l) == 0)
            G(i,l) = -10;
        else
            G(i,l) = 1 * (c~=0);
        	x = i + c;
        	if(x == target)        		
        		G(i,l) = 1000000;
        	end
        end
    end
end
G(target,:) = 0;

%% Get Possible Moves
    function [L_new] = getPossibleMoves(i)
        L_new = zeros(L,1);
        x = i/M;
        y = mod(i,M);
        
        for ll=1:L
            u = controlSpace(ll,2) + (controlSpace(ll,1) * M);
            x = i + u;
            if(x<=MN && x>=1 && MoveMatrix(i,x) == 1)            
                L_new(ll) = 1;
            end
        end
        
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
                        elseif(x_i == M-1 && x_c == 1)
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
        bdrs = [1];%,M+1,-M+1]
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

