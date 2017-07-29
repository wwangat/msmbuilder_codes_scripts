resultdir = './';

old_new_relation=importdata('old_new_start_from_zero')+1; %now we start from zero
data=importdata('main_tic1_ext_tic_micro_assignment_without_trim_start_from_zero.txt');
assignment=data(:,4)+1; %now we start from 1, old index
phi = data(:,1);
psi = data(:,2);
nMacro=5; %we care about first five modes

method1='transpose';
microtpm=importdata(strcat('kcenters_microstate_', method1, '_transmat_.txt'));  %it is row normalized,using the new index to calculate
[right_vectors1, left_vectors1] = plot_2Devector_v2(resultdir, method1, microtpm, old_new_relation, phi, psi, assignment, nMacro);

method2='mle';
microtpm=importdata(strcat('kcenters_microstate_', method2, '_transmat_.txt'));  %it is row normalized,using the new index to calculate
[right_vectors2, left_vectors2] = plot_2Devector_v2(resultdir, method2, microtpm, old_new_relation, phi, psi, assignment, nMacro);

%now we calculate the overlap matrix for these two vectors

Y_matrix = zeros(nMacro, nMacro);
for j = 1:nMacro
    for k = 1:nMacro
	Y_matrix(j,k) = dot(right_vectors1(:, j), right_vectors2(:, k))/(sqrt(dot(right_vectors1(:, j), right_vectors1(:, j)))*sqrt(dot(right_vectors2(:, k), right_vectors2(:, k))));
    end
end

abs(Y_matrix)


Y_matrix_left = zeros(nMacro, nMacro);
for j = 1:nMacro
    for k = 1:nMacro
	Y_matrix_left(j,k) = dot(left_vectors1(:, j), left_vectors2(:, k))/(sqrt(dot(left_vectors1(:, j), left_vectors1(:, j)))*sqrt(dot(left_vectors2(:, k), left_vectors2(:, k))));
    end
end

abs(Y_matrix_left)

