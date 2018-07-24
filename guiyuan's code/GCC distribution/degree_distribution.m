function dist = degree_distribution(B)
degree = full(sum(B));
dist = tabulate(degree);