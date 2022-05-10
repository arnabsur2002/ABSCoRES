% Solar Mapping -- NY and NREL Coordinates
% Arnab Sur
% March 14, 2022. 


Solar_NY = readmatrix('solarGenny22f.csv', 'Range','F2:G6');
Solar_NREL_v1 = readmatrix('solar_meta.csv', 'Range','C2:D315');
Solar_NREL(:,1) = Solar_NREL_v1(:,2);
Solar_NREL(:,2) = Solar_NREL_v1(:,1);

Solar_NY_filt = unique(Solar_NY,"rows");
Solar_NREL_filt = unique(Solar_NREL,"rows");
Mapping_Matrix = [];
for i = 1:length(Solar_NY_filt(:,1))
    for j = 1:length(Solar_NREL_filt(:,1))
        Mapping_Matrix(i,j) = distance('gc',[Solar_NY_filt(i,:)], [Solar_NREL_filt(j,:)]);
    end
end

I = [];
M = [];
for k = 1:length(Solar_NY_filt(:,1))
    M(k) = min(Mapping_Matrix(k,:));
    I(k) = find(Mapping_Matrix(k,:)== M(k));
end 

A = [];
for l = 1: length(Solar_NY_filt(:,1))
    A(l,:) = [Solar_NY_filt(l,:) Solar_NREL_filt(I(l),:)];
end

writematrix(A,'Solar_Mapping_NY_NREL.csv');




