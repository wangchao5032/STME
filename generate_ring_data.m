function [D] = genRingData(cx, cy, r1, r2, grid_size, rand_move_par, label)
[x, y] = meshgrid(-r2:grid_size:r2, -r2:grid_size:r2);
sz = size(x);
data_len = sz(1)*sz(2);
x = reshape(x, [data_len 1]) + (2*rand(data_len, 1)-1) * rand_move_par*grid_size/2;
y = reshape(y, [data_len 1]) + (2*rand(data_len, 1)-1) * rand_move_par*grid_size/2;
r = sqrt(x.*x + y.*y);
invalid_pts = (r < r1 | r> r2);
x(invalid_pts) = [];
y(invalid_pts) = [];

x = x + cx;
y = y + cy;

data_len = length(x);

if data_len < 10
    return;
end

D = zeros(data_len, 5);
D(:, 2) = x;
D(:, 3) = y;
D(:, 4) = label;

end
