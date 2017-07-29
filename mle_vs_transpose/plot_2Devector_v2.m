function [right_vectors, left_vectors] = plot_2Devector_v2(resultdir, method, microtpm, old_new_relation, phi, psi, assignment, nMacro)

nMicro = length(microtpm);

[v d w] = eig(microtpm);
[eigenvalues, sortindex] = sort(diag(d), 'descend');
right_vectors = v(:, sortindex);
left_vectors = w(:, sortindex);

disp('data collected');

RdBu = importdata('/home/wang/RdBu');

RdBu = RdBu.data;

for j = 1:nMacro
    vector = left_vectors(:, j)/max(abs(left_vectors(:, j)));
    if vector(1)<0
        vector = -vector;
    end
    colorvector = floor(1+(vector+1)*127/2.0);
    for k = 1:nMicro  %trimed number
        index = find(assignment == old_new_relation(k,1));  %old index
        colorww = RdBu(colorvector(k), :);
        plot(phi(index), psi(index), '.', 'color', colorww);
        %        plot(phi(1:100), psi(1:100), 'color', RdBu(image_data_left(1:100, 4), :));
        hold on;
    end
    hold off;
    set(gcf,'PaperPositionMode','auto');
    colormap(RdBu);
%    axis([-4 2 -4 2]);
    print(strcat(method, '_left-eigenvector_', num2str(j)), '-dpng');
end

for j = 1:nMacro
    vector = right_vectors(:, j)/max(abs(right_vectors(:, j)));
    if vector(1)<0
        vector = -vector;
    end
    colorvector = floor(1+(vector+1)*127/2.0);
    for k = 1:nMicro
        index = find(assignment == old_new_relation(k,1));
        colorww = RdBu(colorvector(k), :);
        plot(phi(index), psi(index), '.', 'color', colorww);
        %        plot(phi(1:100), psi(1:100), 'color', RdBu(image_data_left(1:100, 4), :));
        hold on;
    end
    hold off;
    set(gcf,'PaperPositionMode','auto');
    colormap(RdBu);
%    axis([-4 2 -4 2]);  %this range will vary for different systems
    print(strcat(method, '_right-eigenvector_', num2str(j)), '-dpng');
end

