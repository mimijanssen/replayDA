function cmap = bluewhitered(n)

if nargin < 1
    n = 256; % default size
end

% Define colors
bottom = [0 0 1];   % blue
middle = [1 1 1];   % white
top    = [1 0 0];   % red

% Interpolate
cmap = [linspace(bottom(1),middle(1),n/2)', ...
        linspace(bottom(2),middle(2),n/2)', ...
        linspace(bottom(3),middle(3),n/2)';
        linspace(middle(1),top(1),n/2)', ...
        linspace(middle(2),top(2),n/2)', ...
        linspace(middle(3),top(3),n/2)'];
end