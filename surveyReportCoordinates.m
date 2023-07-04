%% COORDINATE CHANGE
clear; close all;

data = importSurveyReport("G:\Mi unidad\Proyecto Fracking\Microsismica\ET-2022(h)_Survey_Report_Run07_8.5in @4549.77m TD.xlsx", "Sheet1", [24, Inf]);

[X,Y] = libinski2xy(data.TVDm,data.Incl,data.AzimGrid);

figure
plot3(X,Y,data.TVDm,'o')
f1 = gca(); set(f1,'ZDir','reverse');
zlim([0 3000]);
axis equal tight
figure
plot3(data.EWm,data.NSm,data.TVDm,'o')
f2 = gca(); set(f2,'ZDir','reverse');
% zlim([0 3000]);
axis equal tight


function [X, Y] = libinski2xy(depth, inclination, azimuth)
% Converts Libinski coordinates (depth, inclination, azimuth) to X-Y coordinates (east is X, north is Y)
% Inputs:
% depth - array of depths (in feet)
% inclination - array of inclinations (in degrees)
% azimuth - array of azimuths (in degrees)
% Outputs:
% X - array of X-coordinates (in feet)
% Y - array of Y-coordinates (in feet)

% Convert angles to radians
inc_rad = deg2rad(inclination);
azi_rad = deg2rad(azimuth);

% Compute X and Y coordinates
X = depth .* sin(azi_rad) .* sin(inc_rad);
Y = depth .* cos(azi_rad) .* sin(inc_rad);

end
