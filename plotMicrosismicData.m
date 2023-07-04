%% MICROSISMIC DATA PLOT, ET2021h
% IMPORT DATA
clear; close all; set(0,'DefaultFigureWindowStyle','docked');

pozo = '2022'; fprintf("Well ET%sh \n",pozo);

bonardaTopZ = 2700; % m depth  
deltaZ = 400;

switch pozo
    case '2021'
        data = importMicroData("2021h_final-events_Campo-Inchauspe_Argentina2_849m-Zdatum.xlsx", "2021_final_all_times", [2, Inf]);
        casingData = importSurveyReport2021("ET.xp-2021 Pilot_Survey_Run10_8.5in_3330.07m.xlsx", "Sheet1", [24, Inf]);
        locationGrid = [2464325.24 5852410 0]; %2021 N 5852410.000 m, E 2464325.240 m

        % vector con las stages elegidas
        stages = [1:31]; % plotea todas
        % stages = [1 6 7 8 13 14]; % plotea solo 1 3 y 5
    case '2022'
        data = importMicroData("G:\Mi unidad\Proyecto Fracking\Microsismica\2022h_final-events_Campo-Inchauspe_Argentina2_849m-Zdatum.xlsx", "2022_final_all-times", [2, Inf]);
        casingData = importSurveyReport('ET-2022(h)_Survey_Report_Run07_8.5in @4549.77m TD.xlsx');
        locationGrid = [2464335 5852410 0];

        % vector con las stages elegidas
        stages = [1:15]; % plotea todas
        % stages = [1 6 7 8 13 14]; % plotea solo 1 3 y 5
        
%         % fault event classification
%         % 1
%         log0 = data.Z > bonardaTopZ+51;
%         % 2
%         log1 = (data.Z > bonardaTopZ & data.Z < bonardaTopZ + 51) & (data.STRIKE < 68 | data.STRIKE > 105);
%         log2 = ((data.Z > bonardaTopZ & data.Z < bonardaTopZ + 51) & (data.STRIKE > 68 & data.STRIKE < 105)) & ( (data.X > 2465010 & data.Y < 5850840) | data.X > 2465100 | ( data.X > 2464925 & data.Y > 5851705) );
%         % 3
%         log3 = (data.Z > bonardaTopZ - 41 & data.Z < bonardaTopZ) & ( data.X < 2463960 | data.X > 2465090 | ( data.X > 2464930 & data.Y < 5850760) | data.Y > 5851900 );
% 
%         log = ~(log0 | log1 | log2 | log3);
%         data = data(log,:);

        
end
[data,~,meanStage1] = moveXYZ(data,locationGrid); % traslada sistema de coordenadas al centro del stage 1

%% ELEGIR STAGES PARA PLOTEAR

switchPlot = 'Mw/RAKE'; % 'U' 'D/S' 'both' 'Mw'
plotFilterPLanes = false; 
plotPlanes = false;
%% INITIALIZATION
% SOLO PLOTEA RESERVOIR EVENTS
data = data(data.Fault == 0,:);

Colors = colormap(jet(length(stages)));
i = 1; legendText    = cell(1,3*length(stages));
magFactor = 200; planeRange = 500; plotHandles = [];

f1 = gca();
axis equal tight
plot3(f1,locationGrid(1)-meanStage1(1),locationGrid(2)-meanStage1(2),locationGrid(3),'*','markersize',14,'markerfacecolor',Colors(i,:));
hold on
zlim([2500 3500])

TVD = casingData.TVDm(160:end-10);
TVD = TVD(~isnan(TVD));
profundidadCasingPromedio = mean(TVD);
% [data] = filterZdata(profundidadCasingPromedio+deltaZ,profundidadCasingPromedio-deltaZ,data,f1,planeRange,plotFilterPLanes);
% [data] = filterZdata(bonardaTopZ+deltaZ,bonardaTopZ,data,f1,planeRange,plotFilterPLanes);


for stage_i = stages
    X = data.X(data.STAGE == stage_i); Y = data.Y(data.STAGE == stage_i); Z = data.Z(data.STAGE == stage_i);
    UX = data.UX(data.STAGE == stage_i); UY = data.UY(data.STAGE == stage_i); UZ = data.UZ(data.STAGE == stage_i);
    DIP = data.DIP(data.STAGE == stage_i); STRIKE = data.STRIKE(data.STAGE == stage_i); RAKE = data.RAKE(data.STAGE == stage_i);
    Mw = data.Mw(data.STAGE == stage_i);
    
    switch switchPlot
        case 'U'
            auxHandle = quiver3(f1,X,Y,Z,UX,UY,UZ,'Color',Colors(i,:));
            title(f1,"UX UY UZ plot");
        case 'D/S'
            auxHandle = DipRakePlane(DIP,STRIKE,X,Y,Z,i,Colors,f1);
            title(f1,"DIP STRIKE plot");
        case 'both'
            auxHandle = quiver3(f1,X,Y,Z,UX,UY,UZ,'Color',Colors(i,:));
            DipRakePlane(DIP,STRIKE,X,Y,Z,i,Colors,f1);
        case 'Mw_old'
            % no usar, no es igual a como plotea el informe
            auxHandle = scatter3(f1,X,Y,Z,[],abs(Mw),'filled');
            colorbar(f1,"eastoutside");
            title(f1,"Mw plot");
        case 'Mw/Stage'
            auxHandle = scatter3(f1,X,Y,Z,abs(Mw),Colors(i,:));
            title(f1,"Events are sized by magnitude and colored by stage");
        case 'Mw/RAKE'
            auxHandle = scatter3(f1,X,Y,Z,abs(Mw),abs(RAKE));
%             ticks = [0, 45, 90, 135, 180]; colors = [1 0 0; 1 0.5 0; 0 0.5 1; 1 0.5 0; 1 0 0];
%             colormap(colors);
            colorbar(f1,"eastoutside","Ticks",[0 45 90 135 180]);
            title(f1,"Events are sized by magnitude and colored by stage");
    end
    plotHandles = [plotHandles,auxHandle];
    hold on
    axis equal tight
    
    %% plot ployfit
    [n,V,p] = affine_fit([X Y Z]);
    
    [S1,S2] = meshgrid([-planeRange 0 planeRange]);

%   generate the pont coordinates
    planeX = p(1)+[S1(:) S2(:)]*V(1,:)';
    planeY = p(2)+[S1(:) S2(:)]*V(2,:)';
    planeZ = p(3)+[S1(:) S2(:)]*V(3,:)';

%   plot the plane
if plotPlanes
    surf(f1,reshape(planeX,3,3),reshape(planeY,3,3),reshape(planeZ,3,3),'facecolor',Colors(i,:),'facealpha',0.5);
%     plot plane normal
    quiver3(f1,p(1),p(2),p(3),magFactor*n(1,:),magFactor*n(2,:),magFactor*n(3,:),'Color',Colors(i,:));
end

    % plot avg point for each stage
    plot3(f1,p(1),p(2),p(3),'o','markersize',12,'markerfacecolor',Colors(i,:));

%     C_proj = [X Y Z] - dot([X Y Z] - p.*ones(size([X Y Z])), n'.*ones(size([X Y Z]))) * n';

    %% project points to plane
    C_proj = zeros(length(X),3);
    for ii = 1:length(X)
        C_proj(ii,:) = [X(ii) Y(ii) Z(ii)] - dot([X(ii) Y(ii) Z(ii)] - p, n) * n';
        %         plot3(C_proj(ii,1),C_proj(ii,2),C_proj(ii,3),'o','markersize',12,'markerfacecolor',Colors(i,:));
    end

    %% histogram/heat map plot
    figure('Name',sprintf('Stage %d',stage_i))
    histPlot = histogram2(C_proj(:,1),C_proj(:,2),[20 20],'DisplayStyle','tile','ShowEmptyBins','on');
    xlabel('x'); ylabel('y');
    colormap('default')
    colorbar

    % Density plot
    % es parecido al plot del histograma pero no tan claro, se puede
    % descomentar y ver
    %     figure
    %     f = DataDensityPlot(C_proj(:,1),C_proj(:,2),20);

    fprintf('Stage %d: media en X %.2f SDx %.2f media en Y %.2f SDy %.2f \n',...
        stage_i,mean(C_proj(:,1)),std(C_proj(:,1)),mean(C_proj(:,2)),std(C_proj(:,2)));

    i = i + 1;
end

%% CASING PLOT
EW = casingData.EWm(~isnan(casingData.EWm));
NS = casingData.NSm(~isnan(casingData.NSm));
TVD = casingData.TVDm(~isnan(casingData.TVDm));
plot3(f1,EW+locationGrid(1)-meanStage1(1),NS+locationGrid(2)-meanStage1(2),TVD,'o')

fnd = find(casingData.MDm == 4420,72);
plot3(f1,casingData.EWm(fnd)+locationGrid(1)-meanStage1(1),casingData.NSm(fnd)+locationGrid(2)-meanStage1(2),casingData.TVDm(fnd),'*','Color',[1 0 0]);
casingMatrix = [casingData.EWm+locationGrid(1)-meanStage1(1),casingData.NSm+locationGrid(2)-meanStage1(2),casingData.TVDm];

% plotCluster(f1,casingData,casingMatrix,stages,Colors);
count = 1;
for i = stages
    stageName{count} = sprintf('Stage %d',i);
    count = count + 1;
end

legend(f1,plotHandles,stageName,'location','bestoutside');

%%
% legend(f1,legendText)
xlabel(f1,'x [m]'); ylabel(f1,'y [m]'); zlabel(f1,'depth [m]');
set(f1,'ZDir','reverse')
axes(f1)

%% END
function [data,locationGrid,out] = moveXYZ(data,locationGrid)
    
    stage_i = 1;
    X = data.X(data.STAGE == stage_i); Y = data.Y(data.STAGE == stage_i);

    data.X = data.X - mean(X);
    data.Y = data.Y - mean(Y);

    locationGrid(1:2) = locationGrid(1:2) - [mean(X) mean(Y)];

    out = [mean(X),mean(Y)];
end

function [n,V,p] = affine_fit(X)
    %Computes the plane that fits best (lest square of the normal distance
    %to the plane) a set of sample points.
    %INPUTS:
    %
    %X: a N by 3 matrix where each line is a sample point
    %
    %OUTPUTS:
    %
    %n : a unit (column) vector normal to the plane
    %V : a 3 by 2 matrix. The columns of V form an orthonormal basis of the
    %plane
    %p : a point belonging to the plane
    %
    %NB: this code actually works in any dimension (2,3,4,...)
    %Author: Adrien Leygue
    %Date: August 30 2013
    
    %the mean of the samples belongs to the plane
    p = mean(X,1);
    
    %The samples are reduced:
    R = bsxfun(@minus,X,p);
    %Computation of the principal directions if the samples cloud
    [V,~] = eig(R'*R);
    %Extract the output from the eigenvectors
    n = V(:,1);
    V = V(:,2:end);
end

function [ f ] = DataDensityPlot( x, y, levels )
%DATADENSITYPLOT Plot the data density 
%   Makes a contour map of data density
%   x, y - data x and y coordinates
%   levels - number of contours to show
%
% By Malcolm Mclean
%
    map = dataDensity(x, y, 256, 256);
    map = map - min(min(map));
    map = floor(map ./ max(max(map)) * (levels-1));
    f = figure();
    
    image(map);
    colormap(jet(levels));
    set(gca, 'XTick', [1 256]);
    set(gca, 'XTickLabel', [min(x) max(x)]);
    set(gca, 'YTick', [1 256]);
    set(gca, 'YTickLabel', [min(y) max(y)]);
%     uiwait;
end

function [ dmap ] = dataDensity( x, y, width, height, limits, fudge )
%DATADENSITY Get a data density image of data 
%   x, y - two vectors of equal length giving scatterplot x, y co-ords
%   width, height - dimensions of the data density plot, in pixels
%   limits - [xmin xmax ymin ymax] - defaults to data max/min
%   fudge - the amount of smear, defaults to size of pixel diagonal
%
% By Malcolm McLean
%
    if(nargin == 4)
        limits(1) = min(x);
        limits(2) = max(x);
        limits(3) = min(y);
        limits(4) = max(y);
    end
    deltax = (limits(2) - limits(1)) / width;
    deltay = (limits(4) - limits(3)) / height;
    if(nargin < 6)
        fudge = sqrt(deltax^2 + deltay^2);
    end
    dmap = zeros(height, width);
    for ii = 0: height - 1
        yi = limits(3) + ii * deltay + deltay/2;
        for jj = 0 : width - 1
            xi = limits(1) + jj * deltax + deltax/2;
            dd = 0;
            for kk = 1: length(x)
                dist2 = (x(kk) - xi)^2 + (y(kk) - yi)^2;
                dd = dd + 1 / ( dist2 + fudge); 
            end
            dmap(ii+1,jj+1) = dd;
        end
    end
            
end

function DipRakePlane(DIP,STRIKE,X,Y,Z,i_color,Colors,fig)

    strikeS2 = zeros(size(STRIKE));
    for i = 1:length(DIP)
        if STRIKE(i)>180
            strikeS2(i) = STRIKE(i) - 180;
        else
            strikeS2(i) = STRIKE(i);
        end
    end
    % senoStrike = sind(strikeS2)
    vStrikeS2 = [sind(strikeS2),cosd(strikeS2), zeros(length(strikeS2),1)];

    vDipS2 = [cosd(DIP), zeros(size(DIP,1),1), sind(DIP)];
    normal = zeros(size(STRIKE,1),3);

    for i = 1:size(normal,1)
        n_aux = cross(vDipS2(i,:),vStrikeS2(i,:));
        normal(i,:) = n_aux/norm(n_aux);
    end
    
    quiver3(fig,X,Y,Z,normal(:,1),normal(:,2),normal(:,3),'Color',Colors(i_color,:));
    hold on

end

function [data] = filterZdata(zMax,zMin,data,f1,planeRange,plotPlanesZ,Colors)

    log = (data.Z > zMax | data.Z < zMin);

    data(log,:) = [];

    if plotPlanesZ
    n = 4;
    x = linspace(-planeRange, planeRange, n); y = linspace(-planeRange, planeRange, n);
    [X,Y] = meshgrid(x,y);
    Z1 = zMax*ones(size(X));
    Z2 = zMin*ones(size(X));

    surf(f1,X,Y,Z1,'FaceAlpha',0.3);
    surf(f1,X,Y,Z2,'FaceAlpha',0.3);
    end
end

function plotCluster(f1,casingData,casingMatrix,stages,Colors)

clusters = importClusters("Clusters.csv", [1, Inf]);
casDataMDm = casingData.MDm(~isnan(casingData.MDm));
count = 1;
plotHandles = [];
for i = stages

    closest_vals = interp1(casDataMDm, casDataMDm, clusters(i,:), 'nearest', 'extrap');
    idx = ismember(casingData.MDm, closest_vals);

    stageName{count} = sprintf('Stage %d',i);
    aux = plot3(f1,casingMatrix(idx,1),casingMatrix(idx,2),casingMatrix(idx,3),'*','Color',Colors(count,:),'DisplayName',stageName{count});
    plotHandles = [plotHandles,aux];
    count = count + 1;

end

legend(f1,plotHandles,stageName,'location','bestoutside');

% for i = 1:length(stageName)
% 
% entry = findobj(hL.EntryContainer.Children, 'DisplayName', stageName{i});
% set(entry, 'IconDisplayStyle', '*', 'Color', Colors(i,:));
% 
% end


end