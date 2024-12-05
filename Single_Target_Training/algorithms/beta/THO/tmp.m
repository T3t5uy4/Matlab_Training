fitness = [5, 3, 8, 1];
[fitness, idx] = sort(fitness);

disp('排序后的 fitness:');
disp(fitness); % 输出： [1, 3, 5, 8]

disp('原数组索引:');
disp(idx); % 输出： [4, 2, 1, 3]
