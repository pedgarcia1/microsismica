% IMPORT DATA
clear; close all; set(0,'DefaultFigureWindowStyle','docked');

casing{1} = importSurveyReport('ET-2022(h)_Survey_Report_Run07_8.5in @4549.77m TD.xlsx');
casing{2} = importSurveyReport('ET-2024(h)_Geodetic Survey Report__Run05_TD@4639m.xlsx');
nCasings = size(casing,2);
locationGrid = [0 0 0];
meanStage1 = [0 0 0];

clusters{1} = importClustersNew("ET2022h clusters.csv");
clusters{2} = importClustersNew("ET2024h clusters.csv");

figure; hold on; view(-10,15);
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

vec = casingPoints(3,:) - casingPoints(1,:);
vec = vec/norm(vec);
nvec = cross(vec,[0 0 1]); % vector perpendicular al casing
hold on
quiver3(casingPoints(1,1),casingPoints(1,2),casingPoints(1,3),400*vec(1),400*vec(2),400*vec(3));
quiver3(casingPoints(1,1),casingPoints(1,2),casingPoints(1,3),400*0,400*0,400*1);

for i = 1:3
    quiver3(casingPoints(i,1),casingPoints(i,2),casingPoints(i,3),500*nvec(1),500*vec(2),500*vec(3));
end

aVec = acosd(dot(vec,nvec)/(norm(vec)*norm(nvec)));
atan2d(norm(cross(vec,nvec)),dot(vec,nvec))

resultados.dZ = mean(casingPoints(1:3,3)) - mean(casingPoints(4:6,3));

resultados.dProj = [];
for i = 1:3
    PA = casingPoints(i+3,:)-casingPoints(3,:);
    resultados.dProj(i) = norm(PA - dot(PA,vec)*vec );
end

axis tight 
axis square
