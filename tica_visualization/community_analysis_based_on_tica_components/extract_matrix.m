%%%%%%%%%%%get the tica importance matrices
file='atom1_atom2_from1_tic1_square.contribution';
a = importdata(file);  %max index:1583
index_keep = importdata('indexes_keep_from1');

matrix = zeros(1600);
for j = 1:size(a, 1)
	matrix(a(j, 1), a(j, 2)) = a(j, 3);
end
matrix = matrix+matrix';

sub_matrix = matrix(index_keep, index_keep);
dlmwrite('tic1_importance_matrix.mtx', sub_matrix, 'delimiter', ' ', 'precision', '%20.16f');

