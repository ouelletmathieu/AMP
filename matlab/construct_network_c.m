function [Xu, fV] = construct_network_c(Xs, dXs, Xui, conn, a)
% Function for constructing networks given topology and initial guess
% both for infinitesimal or finite motions
%
% Inputs
% Xs:       d x n       matrix of specified node positions
% dXs:      d x n x z   matrix of z specified node motions
% Xui:      d x k       matrix of initial positions of unspecified nodes
% conn:     s x 2       connection matrix from node i to node j
% a:        1 x 1       1 for finite, 0 for infinitesimal
%
% Outputs
% Xu: 2dxk matrix of unspecified node positions and motions/final positions
% fV: 1xk vector of absolute value of conic at Xu (ideally 0)


% Initial values
d = size(Xs,1);
k = size(Xui,2);
z = size(dXs, 3);

% Iterate across unspecified nodes and optimize position
Xu = zeros(2*d,k);
fV = zeros(1,k);

% Optimization Options
options = optimset('TolFun', 1e-6, 'TolX', 1e-6, 'MaxFunEvals',10000, 'Display', 'off');

% Iterate through each unspecified node initial guess
for i = 1:k
    % Find all specified nodes connected to unspecified node i
    sInds = conn(find(conn(:,2)==(i+size(Xs,2))),1);
    
    % Collect cost functions for solution spaces generated by each motion
    f = cell(1,z);
    for j = 1:z
        [Q, W, v] = construct_conic(Xs(:,sInds), dXs(:,sInds,j), a);
        % Convert conic matrices to common spatial coordinates
        if(z > 1)
            P = [W(1:d,:), v(1:d);...
                 zeros(1,d), 1]^-1;
            Q = P'*Q*P;
        end
        f{j} = @(x) ([x 1] * Q * [x 1]')^2;
    end
    m = size(Q,1)-1;
    % Cost function to find where all conics vanish
    fF = @(x) sum(cellfun(@(F)F(x), f));
    % Convert initial condition to coordinate representation
    if(z == 1)
        cP1 = W(1:d,:)\(Xui(:,i) - v(1:d));
        x0 = cP1(1:m)';
    else
        x0 = Xui(:,i)';
    end
    [XP, fVal] = fminsearch(fF, x0, options);
    if(z > 1)
        Xu(1:d,i) = XP';
    else
        Xu(:,i) = [W(1:(2*d),:) v(1:(2*d))] * [XP 1]';
    end
    fV(i) = sqrt(fVal);
end

end