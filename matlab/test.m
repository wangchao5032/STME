


data = xlsread('data\output-python.xlsx');
types = unique(data(:, 5));
figure(2); clf(2); hold on;
for ti = 1:length(types)
    d = data((data(:, 5) == types(ti)), :);
    fprintf('%d - %d\n', types(ti), length(d))
    plot(d(:, 2), d(:, 3), '.')
end