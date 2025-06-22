% Function to, by a vector of Als (X) and their Y value (Y), adjust
% them to a parabola and take the vertex as the true front of the 
% cell

function [Al_X,Al_Y] = parab(X,Y,side)
    
    % If the vector of Als is empty or any is non a number, Al_X 
    % and Al_y will be 0, the frame is not taken into account
    if (isempty(X) || any(isnan(X)))
        Al_X = 0;
        Al_Y = 0;

    % In case the length of the vector of Als is equal or lower 
    % than 3 then we can not adjust to a parabola, so we take the 
    % max value as the cell front
    elseif (length(X) <= 3)
        [~, idx] = max(X(:,1));
        Al_X = X(idx);
        Al_Y = Y(idx);

    else
        
        % The vectors are adjusted to a parabola and we set up an 
        % epsilon so that if the factor "a" is lower than epsilon 
        % then it is almost a straight line
        p = polyfit(Y, X, 2);
        a = p(1);
        b = p(2);
        
        epsilon = 1e-3;
        
        % Defining of the vertex of the parabola
        Y_vertex = -b / (2 * a);
        X_vertex = polyval(p, Y_vertex);
        
        % If the factor "a" is lower than epsilon, the max value of
        % the vector of Als will be the cell front
        if abs(a) < epsilon
            [~, idx] = max(X(:,1));
            Al = X(idx);
            y = Y(idx);

        % If the factor "a" is greater than 0, the parabola is 
        % oriented backwards, so the frame is not taken into 
        % account, only if calculating front
        elseif side == "front" & a > 0
            Al = 0;
            y = 0;

        % If the factor "a" is lower than 0, the parabola is 
        % oriented forward, so the frame is not taken into account,
        % only if calculating back front
        elseif side == "back" & a < 0
            Al = 0;
            y = 0;

        % If the vertex is a number and the Y of the vertex is 
        % between the channel's edges, the vertex is the front of 
        % the cell
        elseif ((~isnan(X_vertex) && ~isnan(Y_vertex)) && (Y_vertex > min(Y) && Y_vertex < max(Y)))
            Al = X_vertex;
            y = Y_vertex;

        % If non of the above are true, the front of the cell is 
        % the max of the Als in the input vector X or the back
        % front is the min of the Als
        else
            if (side == "front")
                [~, idx] = max(X(:,1));
                Al = X(idx);
                y = Y(idx);
            elseif (side == "back")
                [~, idx] = min(X(:,1));
                Al = X(idx);
                y = Y(idx);
            end
            
        end
        
        % The output is the Al_max and Y_max we got after the 
        % function
        Al_X = Al;
        Al_Y = y;
    end
end