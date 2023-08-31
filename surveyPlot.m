% IMPORT DATA
clear; close all; set(0,'DefaultFigureWindowStyle','docked');

casing{1} = importSurveyReport('ET-2022(h)_Survey_Report_Run07_8.5in @4549.77m TD.xlsx');
casing{2} = importSurveyReport('ET-2024(h)_Geodetic Survey Report__Run05_TD@4639m.xlsx');
nCasings = size(casing,2);
locationGrid = [2464335 5852410 0];
meanStage1 = [0 0 0];

clusters{1} = importClustersNew("ET2022h clusters.csv");
clusters{2} = importClustersNew("ET2024h clusters.csv");

figure; hold on;
f1 = gca();
xlabel(f1,'x [m]'); ylabel(f1,'y [m]'); zlabel(f1,'depth [m]');
set(f1,'ZDir','reverse');
for i = 1:nCasings
    EW = casing{i}.EWm(~isnan(casing{i}.EWm));
    NS = casing{i}.NSm(~isnan(casing{i}.NSm));
    TVD = casing{i}.TVDm(~isnan(casing{i}.TVDm));
    plot3(f1,EW+locationGrid(1)-meanStage1(1),NS+locationGrid(2)-meanStage1(2),TVD,'o')
end
hold on
c = 1;
for i = 1:size(clusters,2)
    EW = casing{i}.EWm(~isnan(casing{i}.EWm));
    NS = casing{i}.NSm(~isnan(casing{i}.NSm));
    TVD = casing{i}.TVDm(~isnan(casing{i}.TVDm));
    MD = casing{i}.MDm(~isnan(casing{i}.TVDm));
    for n = 0:2
        cTVD = mean(clusters{1},2);
        meanTVD = mean(cTVD(10*n+1:10*(n+1)));
        I = find(abs(MD-meanTVD) == min(abs(MD-meanTVD)));
        casingPoints(c,:) = [EW(I) NS(I) TVD(I)];
        c = c + 1;
    end
end
plot3(f1,casingPoints(:,1)+locationGrid(1)-meanStage1(1),casingPoints(:,2)+locationGrid(2)-meanStage1(2),casingPoints(:,3),'+','Color',[0 0 0])
for i = 1:size(casingPoints,1)
    text(casingPoints(i,1)+locationGrid(1)-meanStage1(1), casingPoints(i,2)+locationGrid(2)-meanStage1(2), casingPoints(i,3)+200, num2str(i), 'FontSize', 10);
end
for i = [1 2 4 5]
    resultados.dStages(i) = norm(casingPoints(i,:)-casingPoints(i+1,:));
end
for i = 1:size(casingPoints,1)/2
    resultados.dCasings(i) = norm(casingPoints(i,:)-casingPoints(i+3,:));
end