resultdir = './';

nMicro=36;
phi = importdata('../phi_angles.txt');
psi = importdata('../psi_angles.txt');
assignment = importdata('MicroAssignment.txt');

TPM=importdata('microstate_TPM_transpose_at_1ps.txt');
[u,v,w]=eig(TPM);
[sortvalue,sortindex]=sort(diag(v),'descend');
u=u(:,sortindex);
w=w(:,sortindex);
v=v(:,sortindex);

nMacro=4;
%right_vectors = importdata(strcat(prefix, '_membership.txt'));
right_vectors=u;

prefix='rightev';
disp('data collected');

RdBu = importdata('/sda2/wang/RdBu');

RdBu = RdBu.data;

for j = 1:nMacro
    vector = right_vectors(:, j)/max(abs(right_vectors(:,j)));  %for the eigenvector, which lies in [-1,1]
    colorvector = floor(1+(vector+1)*127/2.0); %for the eigenvector
%    colorvector = floor(63*vector+64); %only for the membership function,since it is between 0 to 1
    nkk=1;
    for k=1:nMicro
        index = find(assignment == k);
        if(length(index)~=0)
            colorww = RdBu(colorvector(nkk), :);
            plot(phi(index(1:5:end)), psi(index(1:5:end)), '.', 'color', colorww);
            hold on;
	    nkk=nkk+1;
        end
    end
    hold off;
    axis([min(phi) max(phi) min(psi) max(psi)]);
%    axis([-180 180 -180 180]);
    set(gcf,'PaperPositionMode','auto');
    print(strcat(resultdir, '/',prefix, '_state_', num2str(j)), '-dpng');
end

prefix='leftev';
left_vectors=w;
for j = 1:nMacro
    vector = left_vectors(:, j)/max(abs(left_vectors(:,j)));  %for the eigenvector, which lies in [-1,1]
    colorvector = floor(1+(vector+1)*127/2.0); %for the eigenvector
%    colorvector = floor(63*vector+64); %only for the membership function,since it is between 0 to 1
    nkk=1;
    for k=1:nMicro
        index = find(assignment == k);
        if(length(index)~=0)
            colorww = RdBu(colorvector(nkk), :);
            plot(phi(index(1:5:end)), psi(index(1:5:end)), '.', 'color', colorww);
            hold on;
	    nkk=nkk+1;
        end
    end
    hold off;
    axis([min(phi) max(phi) min(psi) max(psi)]);
%    axis([-180 180 -180 180]);
    set(gcf,'PaperPositionMode','auto');
    print(strcat(resultdir, '/',prefix, '_state_', num2str(j)), '-dpng');
end

