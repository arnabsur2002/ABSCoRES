% Wind Mapping -- NY and NREL Coordinates
% Arnab Sur
% March 14, 2022. 


Wind_NY = readmatrix('windGenny22.csv', 'Range','F2:G28');
Wind_NREL_v1 = readmatrix('wind_meta.csv', 'Range','C2:D81');
Wind_NREL(:,1) = Wind_NREL_v1(:,2);
Wind_NREL(:,2) = Wind_NREL_v1(:,1);

Wind_NY_filt = unique(Wind_NY,"rows");
Wind_NREL_filt = unique(Wind_NREL,"rows");
Mapping_Matrix = [];
for i = 1:length(Wind_NY_filt(:,1))
    for j = 1:length(Wind_NREL_filt(:,1))
        Mapping_Matrix(i,j) = distance('gc',[Wind_NY_filt(i,:)], [Wind_NREL_filt(j,:)]);
    end
end

I = [];
M = [];
for k = 1:length(Wind_NY_filt(:,1))
    M(k) = min(Mapping_Matrix(k,:));
    I(k) = find(Mapping_Matrix(k,:)== M(k));
end 

A = [];
for l = 1: length(Wind_NY_filt(:,1))
    A(l,:) = [Wind_NY_filt(l,:) Wind_NREL_filt(I(l),:)];
end

writematrix(A,'Wind_Mapping_NY_NREL.csv');


