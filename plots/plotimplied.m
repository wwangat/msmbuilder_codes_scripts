function  [data_to_plot] = plotimplied(resultdir, trajMicro, traj_num, traj_len, nMicro, starting, interval, ending) %for trajMicro, should start from 1
	%when sample is not enough, eigenvalus can be negative, thus implied timescale would look very weird
	disp('You should judge whether Assignment is starting from 1, if not, revise your input data');
	traj = trajMicro;
	timescale = starting:interval:ending;
	num = length(timescale);
	if nMicro>10
		num_eigs = 10;
	else
		num_eigs = nMicro-1;
	end
	data_to_plot = zeros(num_eigs, num);
	for time = 1:num
		tau = timescale(time);
		CountMatrix = zeros(nMicro, nMicro);
		for j = 1:traj_num
	 		trajectory = traj{j};
	 		for k = 1:traj_len(j)-tau
	 			old_index = trajectory(k);
	 			new_index = trajectory(k+tau); %start from 1
                if (old_index >=0 & new_index >=0)
	 			    CountMatrix(new_index, old_index) = CountMatrix(new_index, old_index)+1;
	 		    end
            end
	 	end
		transmatrix=(CountMatrix+transpose(CountMatrix))/2.0;
		%%original transition count matrix
		%Need to check whether there are zero transition microstates and remove them in both CountMatrix and original assignment file
		index = find(sum(transmatrix, 1) == 0);
		if length(index)~=0
		    transmatrix(index, :) = [];
		    transmatrix(:, index) = [];
		end
		sizeMicro = length(transmatrix);
		%calculate transition probability matrix, at this moment, col normalize
		TPM = zeros(sizeMicro, sizeMicro);
		for j = 1:sizeMicro
			colsum = sum(transmatrix(:, j));
			TPM(:, j) = transmatrix(:, j)/colsum;
		end
		%column normalized
		disp(strcat('for lagtime ', num2str(tau/10),' ps, TPM is of size ', num2str(sizeMicro)))
		evalues = sort(real(eig(TPM)),'descend');
%		disp(strcat('eigenvalue at timestep', num2str(tau)));
		evalues = evalues(2:end);
%		evalues(1:5)
        if (sizeMicro<10)		
            num_eigs = sizeMicro-1;
        end
		data_to_plot(1:num_eigs, time) = -tau./log(evalues(1:num_eigs));
	end
	%now begin to plot implied timescale
	disp('begin to plot implied timescale');
	for j = 1:num_eigs
		semilogy(timescale*1e-1, data_to_plot(j, :)*1e-1,'o-'); %for alanine dipeptide, unit:ps
		hold on;
	end
	hold off;
	xlabel('\tau (ps)');
	ylabel('\tau_i');
%	ylim([0 1000]);
	title('ImpliedTimescale');
end
