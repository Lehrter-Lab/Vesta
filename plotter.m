%% CMikolaitis @ USA/DISL, 2025
warning('off','MATLAB:table:ModifiedAndSavedVarnames');

%% Parameters
list = ["USGS*","node*","NOAA*"];
idToNode = table;
idToNode.Validation = [01304200; 01304562; 01304650; 8510560];
idToNode.Model = [48505; 13488; 42773; 80345];
idToNode.Name = ["Orient";"Peconic";"Shelter";"Montauk"];
baseTime = datetime(2022,1,1);

%% File parsing
allVarNames = {'Time', 'Temperature', 'Salinity', 'Elevation'};
vTemp = struct(); mTemp = struct();
vSal  = struct(); mSal  = struct();
vElev = struct(); mElev = struct();

for j = 1:length(list)
    folder = dir(list(j));
    for i = 1:length(folder)
        filename = folder(i).name;
        parts = split(filename, "_");
        % Extract site ID and match to the correct name
        if contains(filename, "USGS") || contains(filename, "NOAA")
            siteID = str2double(parts{2});
            row = find(idToNode.Validation == siteID);
        elseif contains(filename, "node")
            siteID = str2double(parts{2});
            row = find(idToNode.Model == siteID);
        end
        siteName = idToNode.Name{row};

        % Read and process the file
        if contains(filename, 'USGS')
            t = readtable(filename, "FileType", "text", 'Delimiter', '\t');
            time = t.datetime;
            temp = t{:,7};
            sal  = t{:,9};
            elev = t{:,5};
            vTemp.(char(siteName+"_Temperature")) = temp;
            vTemp.(char(siteName+"_Time")) = time;
            vSal.(char(siteName+"_Salinity")) = sal;
            vSal.(char(siteName+"_Time")) = time;
            vElev.(char(siteName+"_Elevation")) = elev*0.3048;
            vElev.(char(siteName+"_Time")) = time;
    
        elseif contains(filename, 'node')
            tRaw = readtable(filename, "FileType", "text", 'Delimiter', '\t');
            if contains(filename, 'temperature')
                t = tRaw(tRaw.vgrid_layer == 32, :);  % Filter for top layer
                time = baseTime + seconds(t.time);
                mTemp.(char(siteName+"_Temperature")) = t{:,4};
                mTemp.(char(siteName+"_Time")) = time;
            elseif contains(filename, 'salinity')
                t = tRaw(tRaw.vgrid_layer == 32, :);  % Filter for top layer
                time = baseTime + seconds(t.time);
                mSal.(char(siteName+"_Salinity")) = t{:,4};
                mSal.(char(siteName+"_Time")) = time;
            elseif contains(filename, 'elevation')
                t = tRaw(tRaw.vgrid_layer == 1, :);
                time = baseTime + seconds(t.time);
                mElev.(char(siteName+"_Elevation")) = t{:,4}*-1;
                mElev.(char(siteName+"_Time")) = time;
            end
        elseif contains(filename, 'NOAA')
            t = readmatrix(filename);
            time = datetime(t(:,1), t(:,2), t(:,3), t(:,4), t(:,5), t(:,6));
            if contains(filename, 'temperature')
                vTemp.(char(siteName+"_Temperature")) = t(:,7);
                vTemp.(char(siteName+"_Time")) = time;
            elseif contains(filename, 'elevation')
                vElev.(char(siteName+"_Elevation")) = t(:,7)*0.3048;
                vElev.(char(siteName+"_Time")) = time;
            end
        end
    end
end
%% Sanity
clearvars -except m* v* idToNode
%% Plot
t1 = datetime(2022, 1, 1); t2 = datetime(2022, 12, 31);
plotTileComparison(vTemp, mTemp, idToNode.Name, 'Temperature', '°C', [t1 t2]);
plotTileComparison(vSal, mSal, idToNode.Name, 'Salinity', 'PSU', [t1 t2]);
plotTileComparison(vElev, mElev, idToNode.Name, 'Elevation', 'm', [t1 t2]);
%% Plot function
function plotTileComparison(vStruct, mStruct, siteNames, variableLabel, yLabel, timeRange)
% Parameters:
%   vStruct       - struct containing validation data (e.g., vTemp, vSal)
%   mStruct       - struct containing model data (e.g., mTemp, mSal)
%   siteNames     - cell array of site names (e.g., idToNode.Name)
%   variableLabel - string like 'Temperature', 'Salinity', or 'Elevation'
%   yLabel        - string label for Y-axis (e.g., '°C', 'PSU', 'm')
%   timeRange     - optional 1x2 datetime array for x-axis limits (e.g., [startTime, endTime])

    if nargin < 6
        timeRange = [];  % No limit by default
    end
    nSites = numel(siteNames);
    figure('Name', variableLabel + " Comparison", 'Position', [100, 100, 1200, 600]);
    tiledlayout('flow');
    for i = 1:nSites
        site = siteNames{i};
        fieldVal = site + "_" + variableLabel;
        fieldTime = site + "_Time";
        hasV = isfield(vStruct, fieldVal) && isfield(vStruct, fieldTime);
        hasM = isfield(mStruct, fieldVal) && isfield(mStruct, fieldTime);
        if hasV || hasM
            nexttile;
            hold on;
            if hasV
                plot(vStruct.(fieldTime), vStruct.(fieldVal), 'b-', 'DisplayName', 'Validation');
            end
            if hasM
                plot(mStruct.(fieldTime), mStruct.(fieldVal), 'r--', 'DisplayName', 'Model');
            end
            title(site + " " + variableLabel);
            xlabel("Time");
            ylabel(yLabel);
            legend('Location', 'best');
            grid on;
            if ~isempty(timeRange)
                xlim(timeRange);
            end
        end
    end
end
