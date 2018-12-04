%%%%%%%%%%%evaluate when the tica eigenvectors become invariant%%%%%%%%%%%%%%%%%%%%
covariance_matrix = importdata('tica_covariance_matrix_corr280.txt'); % calculated using the lagtime=10 data
%comp1 = 1;
%comp2 = 2;
ev1 = importdata('eigvector_corr280.txt');   %reference
temp = [];
for lagtime = 20:20:280

filename = strcat('eigvector_corr', num2str(lagtime), '.txt');
ev2 = importdata(filename);

angle = zeros(10);
angle1 = zeros(10);
inner = zeros(10);
for comp1 = 1:10
    for comp2 = 1:10
	angle(comp1, comp2) = 180/pi*acos(abs(ev1(:, comp1)'*covariance_matrix*ev2(:, comp2))/(sqrt(ev1(:, comp1)'*covariance_matrix*ev1(:, comp1))*sqrt(ev2(:, comp2)'*covariance_matrix*ev2(:, comp2))));
	inner(comp1, comp2) = abs(ev1(:, comp1)'*covariance_matrix*ev2(:, comp2))/(sqrt(ev1(:, comp1)'*covariance_matrix*ev1(:, comp1))*sqrt(ev2(:, comp2)'*covariance_matrix*ev2(:, comp2)));
	angle1(comp1, comp2) = real(180/pi*acos(abs(ev1(:, comp1)'*covariance_matrix*ev2(:, comp2))/(sqrt(ev1(:, comp1)'*covariance_matrix*ev1(:, comp1))*sqrt(ev2(:, comp2)'*covariance_matrix*ev2(:, comp2)))));
    end
end

temp = [temp;inner(1,1) inner(2,2) inner(3,3) inner(4,4) inner(5,5)];
end


