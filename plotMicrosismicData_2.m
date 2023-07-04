%% MICROSISMIC DATA PLOT, ET2022h
% IMPORT DATA
clear; close all; set(0,'DefaultFigureWindowStyle','docked');

data = importMicroData("G:\Mi unidad\Proyecto Fracking\Microsismica\2022h_final-events_Campo-Inchauspe_Argentina2_849m-Zdatum.xlsx", "2022_final_all-times", [2, Inf]);
data = moveXYZ(data); % traslada sistema de coordenadas al centro del stage 1

%% ELEGIR STAGES PARA PLOTEAR
% vector con las stages elegidas
stages = [1 8 14];
switchPlot = 'U'; % 'U' 'D/S' 'both'

%% INITIALIZATION
process_stages = [6 7 8 14]; % seba dijo 13 OK
Colors = colormap(jet(length(stages)));
i = 1; legendText = cell(1,3*length(stages));
magFactor = 200; planeRange = 500;

%% statistical parameters
num_clusters = 1*ones(1,15);
num_clusters(6) = 2; k_lof(6) = 30; lof_threshold(6) = 1.3; % para 6
num_clusters(7) = 1; k_lof(7) = 30; lof_threshold(7) = 2.1; % para 7
% hace un solo cluster aunque elija mas de uno 
num_clusters(8) = 5; k_lof(8) = 30; lof_threshold(8) = 1.3; % para 8
num_clusters(13) = 1; k_lof(13) = 30; lof_threshold(13) = 1; % para 13
num_clusters(14) = 3; k_lof(14) = 30; lof_threshold(14) = 1.5; % para 14

f1 = gca();
axis equal tight
for stage_i = stages
    X = data.X(data.STAGE == stage_i); Y = data.Y(data.STAGE == stage_i); Z = data.Z(data.STAGE == stage_i);
    UX = data.UX(data.STAGE == stage_i); UY = data.UY(data.STAGE == stage_i); UZ = data.UZ(data.STAGE == stage_i);
    DIP = data.DIP(data.STAGE == stage_i); STRIKE = data.STRIKE(data.STAGE == stage_i);

    %% process data
    if any(stage_i == process_stages)
        [cluster_labels, outlier_labels] = kmeans_clustering_with_outliers([X Y Z], num_clusters(stage_i), k_lof(stage_i), lof_threshold(stage_i),stage_i);
    else
        outlier_labels = zeros(length(X),1);
    end

    X = X(~outlier_labels); Y = Y(~outlier_labels); Z = Z(~outlier_labels);
    UX = UX(~outlier_labels); UY = UY(~outlier_labels); UZ = UZ(~outlier_labels);
    DIP = DIP(~outlier_labels); STRIKE = STRIKE(~outlier_labels);

    switch switchPlot
        case 'U'
            quiver3(f1,X,Y,Z,UX,UY,UZ,'Color',Colors(i,:));
            title(f1,"UX UY UZ plot");
        case 'D/S'
            DipRakePlane(DIP,STRIKE,X,Y,Z,i,Colors,f1);
            title(f1,"DIP STIKE plot");
        case 'both'
            quiver3(f1,X,Y,Z,UX,UY,UZ,'Color',Colors(i,:));
            DipRakePlane(DIP,STRIKE,X,Y,Z,i,Colors,f1);
    end
    hold on
    axis equal tight

    %% legend
    switch switchPlot
        case {'U','D/S'}
            legendText{4*i-3} = sprintf('Stage %d',stage_i);
            legendText{4*i-2} = ''; legendText{4*i-1} = ''; legendText{4*i} = '';
        case 'both'
            legendText{5*i-4} = sprintf('Stage %d',stage_i);
            legendText{5*i-3} = ''; legendText{5*i-2} = ''; legendText{5*i-1} = ''; legendText{5*i} = '';
    end

    %% plot ployfit
    [n,V,p] = affine_fit([X Y Z]);

    [S1,S2] = meshgrid([-planeRange 0 planeRange]);

    %   generate the pont coordinates
    planeX = p(1)+[S1(:) S2(:)]*V(1,:)';
    planeY = p(2)+[S1(:) S2(:)]*V(2,:)';
    planeZ = p(3)+[S1(:) S2(:)]*V(3,:)';

    axes(f1);
    hold on;
    %   plot the plane
    surf(f1,reshape(planeX,3,3),reshape(planeY,3,3),reshape(planeZ,3,3),'facecolor',Colors(i,:),'facealpha',0.5);
    % plot plane normal
    quiver3(f1,p(1),p(2),p(3),magFactor*n(1,:),magFactor*n(2,:),magFactor*n(3,:),'Color',Colors(i,:));
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
    figure('Name',sprintf('Hist: Stage %d',stage_i));
    f2 = gca();
    histPlot = histogram2(f2,C_proj(:,1),C_proj(:,2),[20 20],'DisplayStyle','tile','ShowEmptyBins','on');
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
legend(f1,legendText)
xlabel(f1,'x [m]'); ylabel(f1,'y [m]'); zlabel(f1,'depth [m]');
set(f1,'ZDir','reverse')
axes(f1);

function data = moveXYZ(data)

stage_i = 1;
X = data.X(data.STAGE == stage_i); Y = data.Y(data.STAGE == stage_i);

data.X = data.X - mean(X);
data.Y = data.Y - mean(Y);

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

