b = importdata('/home/wang/gwang-PolII/gwang_hpc2_data/round2-analysis/20181119_projected_tICA_analysis_basedon_cac_removepoortrajs/cac/analysis_bootstrap_results/tic_projections/evaluate_tics/ward/lumping_into_3states.txt');
matrix1 = zeros(3);
for j = 0:2
	index1 = find(b == j);
	for k = 0:2
		index2 = find(b == k);
		matrix1(j+1, k+1) = sum(sum(sub_matrix(index1, index2)));
	end
end

matrix1
