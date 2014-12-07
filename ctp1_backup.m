function P = ComputeTransitionProbabilitiesI( stateSpace, ...
    controlSpace, disturbanceSpace, mazeSize, walls, targetCell )
%COMPUTETRANSITIONPROBABILITIESI Compute transition probabilities.
% 	Compute the transition probabilities between all states in the state
%   space for all attainable control inputs.
%
%   P = ComputeTransitionProbabilitiesI(stateSpace, controlSpace,
%   disturbanceSpace, mazeSize, walls, targetCell) computes the transition
%   probabilities between all states in the state space for all attainable
%   control inputs.
%
%   Input arguments:
%
%       stateSpace:
%           A (MN x 2) matrix, where the i-th row represents the i-th
%           element of the state space. Note that the state space also
%           contains the target cell, in order to simplify state indexing.
%
%       controlSpace:
%           A (L x 2) matrix, where the l-th row represents the l-ththe
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
%       P:
%           A (MN x MN x L) matrix containing the transition probabilities
%           between all states in the state space for all attainable
%           control inputs. The entry P(i, j, l) represents the transition
%           probability from state i to state j if control input l is
%           applied.

% put your code here
%% Constants
MN  = size(stateSpace, 1);
M   = mazeSize(1);
N   = mazeSize(2);
S   = size(disturbanceSpace,1);
L   = size(controlSpace,1);
p_pc= 1/L; %probability of each control space
p_pc_norm = 0;
Walls = convertWallMatrix();
Walled = getWalledMatrix();
%% Calculating the transition probability tatble
P = zeros(MN,MN,L);
% Need to add the changes caused due to walls.
for i=1:MN
    [controlSpaceNew,p_pc] = getPossibleMoves(i);
    for l=1:L
        if(controlSpaceNew(l) == 0)
            continue;
        end
        
        for s=1:S
            % control movement % disturbance movement
            x_c = controlSpace(l,2) + disturbanceSpace(s,2);
            y_c = controlSpace(l,1) + disturbanceSpace(s,1);
            uw = y_c + (x_c*M);
            %dynamics of the system
            x = i + uw;
            if( x>MN || x<1)
                continue;
                % if this increases the probability of possible 
                % control spaces then add to them here
            end
            p_pc_norm = p_pc_norm +1;
            P(i,x,l) = P(i,x,l) + p_pc*disturbanceSpace(s,3);
        end
    end
end

%% Get Possible Moves
    function [L_new,p_pc] = getPossibleMoves(i)
        L_new = zeros(L,1);
        for ll=1:L
            u = controlSpace(ll,2) + (controlSpace(ll,1) * M);
            x = i + u;
            if(x<=MN && x>=1)
                if(controlSpace(ll,2) ~= 0 && controlSpace(ll,1) ~=0)
                    if(Walled(x) == 0)
                        temp = controlSpace(ll,2);
                        if(temp >= 1 && temp <=MN && Walls(i,temp) ==0)
                            temp = controlSpace(ll,1) * M;
                            if(temp >= 1 && temp <=MN && Walls(i,temp) ==0)
                                L_new(ll) = 1;
                            end
                        end
                    else
                        L_new(ll) = 0;
                    end
                else
                    if(Walls(i,x) == 0)
                        L_new(ll) = 1;
                    else
                        L_new(ll) = 0;
                    end
                end
            else
                L_new(ll) = 0;
            end
        end
        p_pc = (1/length(L_new==1));
    end
%% Chk if node has wall
    function [Walled] = getWalledMatrix()
        sides = [1,-1,M,-M];
        Walled = zeros(MN,1);
        for i = 1:MN
            for cur = 1:length(sides)
                x_cur = i+sides(cur);
                if(x_cur<=MN && x_cur>=1 && Walls(i,x_cur) == 1)
                    Walled(i) = 1;
                end
            end
        end
    end
%% change Wall Matrix
    function [Walls] = convertWallMatrix()
        WL = size(walls,2)/2;
        Walls = zeros(MN,MN);
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
                        Walls(from,to) = 1;
                        Walls(to,from) = 1;
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
                        Walls(from,to) = 1;
                        Walls(to,from) = 1;
                        from = from - 1;
                    end
                    from    = (from_y)  + (min(to_x,from_x)*M);
                    to = to + 1;
                end
            end
            
        end        
    end
end