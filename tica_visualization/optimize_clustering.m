b = importdata('ward/lumping_into_4states.txt');

sub_matrix = importdata('tic1_importance_matrix.mtx');
n_state = max(b)+1;
matrix1 = zeros(n_state);
for j = 0:max(b)
	index1 = find(b == j);
	for k = 0:max(b)
		index2 = find(b == k);
		matrix1(j+1, k+1) = sum(sum(sub_matrix(index1, index2)));
	end
end

initial_lump = b;
initial_matrix = matrix1;
initial_value = sum(sum(diag(matrix1))); %we want to minimize this value

best_lump = b;
best_value = initial_value;
best_matrix = b;
temp_lump = b;
%begin to randomly permute and optimize
for trial = 1:10000
	trial
	pick = randi([1 size(sub_matrix, 1)]);
	bian = randi([0 max(b)]);
	new_matrix = zeros(n_state);
	temp_lump(pick) = bian;
	for j = 0:max(b)
		index1 = find(temp_lump == j);
		for k = 0:max(b)
			index2 = find(temp_lump == k);
			new_matrix(j+1, k+1) = sum(sum(sub_matrix(index1, index2)));
		end
	end
	if (sum(sum(diag(new_matrix)))<best_value)
		disp('accept');
		best_lump(pick) = bian;
		best_value = sum(sum(diag(new_matrix)));
		best_matrix = new_matrix;
	else
		temp_lump(pick) = best_lump(pick);
	end
end

best_value
best_matrix

dlmwrite('optimized_lump.lump', best_lump);
