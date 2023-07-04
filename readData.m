clc
clearvars

Table = readtable('2022h_final-events_Campo-Inchauspe_Argentina2_849m-Zdatum.xlsx');
FilteredTable = Table(Table.STAGE == 3,:);

strikeS2 = zeros(size(FilteredTable.STRIKE,1),1);
for i = 1:size(FilteredTable.STRIKE,1)
    if Table.STRIKE(i)>180
        strikeS2(i) = FilteredTable.STRIKE(i) - 180;
    else
        strikeS2(i) = FilteredTable.STRIKE(i);       
    end
end
% senoStrike = sind(strikeS2)
vStrikeS2 = [ sind(strikeS2),cosd(strikeS2), zeros(size(strikeS2,1),1)];

vDipS2 = [cosd(FilteredTable.DIP), zeros(size(FilteredTable.DIP,1),1), sind(FilteredTable.DIP)];
normal = zeros(size(FilteredTable.STRIKE,1),3);

for i = 1:size(normal,1)
    normal(i,:) = cross(vDipS2(i,:),vStrikeS2(i,:));
end

figure 
quiver3(FilteredTable.X-2464000,FilteredTable.Y-5854000,FilteredTable.Z,normal(:,1),normal(:,2),normal(:,3),1)
hold on
plot3(FilteredTable.X-2464000,FilteredTable.Y-5854000,FilteredTable.Z,'ro')
set(gca, 'ZDir','reverse')

figure 
plot3(Table.X,Table.Y,Table.Z,'ro')