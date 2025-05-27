%% CMikolaitis @ USA/DISL, 2025
warning('off','MATLAB:table:ModifiedAndSavedVarnames');
warning('off','MATLAB:print:ContentTypeImageSuggested');
%% Parameters
path = "./Source/";
list = ["USGS*","node*","NOAA*"];
% Mapping Table
idToNode            = table;
idToNode.Name       = ["Orient";"Peconic";"Shelter";"Montauk"];
idToNode.Validation = [01304200; 01304562; 01304650; 8510560];
idToNode.Model      = [48821; 13490; 43039; 80775]; % [48505; 13488; 42773; 80345]
% Offset for differing datums
NGVD29toNAVD88 = -0.95*0.3048; % Based on 2025 Peconic River data
NAVD88toMSL    = -0.101;       % Based on 1983 Epoch Montauk NOAA data
% Validation Window
baseYear = 2021;
baseTime = datetime(baseYear,1,1,'TimeZone','UTC');
t1 = datetime(baseTime,'TimeZone','America/New_York'); 
t2 = datetime(baseYear, 12, 31,'TimeZone','America/New_York');
times = timerange(t1,t2);
% Variable info
varInfo = table;
varInfo.Name = ["Temperature";"Salinity";"Elevation"];
varInfo.Unit = ["°C";"PSU";"m"];
% Figure outputs
saveDir = "./";
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end
%% File parsing
vTemp = struct(); mTemp = struct();
vSal  = struct(); mSal  = struct();
vElev = struct(); mElev = struct();

boolM = string(idToNode.Model);
boolV = string(idToNode.Validation);
for j = 1:length(list)
    folder = dir(path+list(j));
    for i = 1:length(folder)
        filename = path+folder(i).name;
        parts = split(filename, "_");
        if ~contains(filename,boolM) && ~contains(filename,boolV)
            continue
        end
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
        if contains(filename, 'USGS') % NGVD29
            t = readtable(filename, "FileType", "text", 'Delimiter', '\t');
            time = datetime(t.datetime,'TimeZone','America/New_York');
            temp = t{:,7};
            sal  = t{:,9};
            elev = t{:,5};
            vTemp.(char(siteName+"_Temperature")) = temp;
            vTemp.(char(siteName+"_Time")) = time;
            vSal.(char(siteName+"_Salinity")) = sal;
            vSal.(char(siteName+"_Time")) = time;
            vElev.(char(siteName+"_Elevation")) = elev*0.3048+NGVD29toNAVD88+NAVD88toMSL;
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
        elseif contains(filename, 'NOAA') % MSL
            t = readmatrix(filename);
            time = datetime(t(:,1), t(:,2), t(:,3), t(:,4), t(:,5), t(:,6),'TimeZone','America/New_York');
            if contains(filename, 'temperature')
                vTemp.(char(siteName+"_Temperature")) = t(:,7);
                vTemp.(char(siteName+"_Time")) = time;
            elseif contains(filename, 'elevation')
                vElev.(char(siteName+"_Elevation")) = t(:,7);
                vElev.(char(siteName+"_Time")) = time;
            end
        end
    end
end
%% Sanity
clearvars -except m* v* idToNode base* times varInfo saveDir
%% Generate time filtered paired datasets
pairTemp = caller(vTemp,mTemp,idToNode.Name,'Temperature',times,baseTime);
pairSal  = caller(vSal, mSal, idToNode.Name,'Salinity',   times,baseTime);
pairElev = caller(vElev,mElev,idToNode.Name,'Elevation',  times,baseTime);
%% Get elevation bias
bias = elevBias(pairElev);
biasCorrection = bias.Montauk;
%% Plot
%plotTileComparison(vTemp,mTemp,idToNode.Name,varInfo.Name(1),varInfo.Unit(1),times,biasCorrection,saveDir);
%plotTileComparison(vSal, mSal, idToNode.Name,varInfo.Name(2),varInfo.Unit(2),times,biasCorrection,saveDir);
%plotTileComparison(vElev,mElev,idToNode.Name,varInfo.Name(3),varInfo.Unit(3),times,biasCorrection,saveDir);
%% Master Func
function pairTable = caller(vStruct,mStruct,siteNames,variableLabel,times,baseTime)
    for i = 1:numel(siteNames)
        site          = siteNames{i};
        fieldVal      = site + "_" + variableLabel;
        fieldTime     = site + "_Time";
        vTable        = table;
        mTable        = table;
        try
            vTable.(site) = vStruct.(fieldVal);
            vTable.Time   = vStruct.(fieldTime);
            Validation    = table2timetable(vTable,'RowTimes','Time');
        catch
            warning([site ' is missing ' variableLabel ' data'])
            Validation = timetable('RowTimes',baseTime);
        end
        mTable.(site) = mStruct.(fieldVal);
        mTable.Time   = mStruct.(fieldTime);
        Model         = table2timetable(mTable,'RowTimes','Time');
        pairCurrent   = synchronize(Validation,Model);
        pairCurrent   = pairCurrent(times,:); 
        if i == 1
            pairTable = pairCurrent;
        else
            pairTable = synchronize(pairTable,pairCurrent);
        end
    end
end
%% Bias Function
function biasTable = elevBias(pairElev)
columns = pairElev.Properties.VariableNames;
    for i = 1:(width(pairElev)/2)
        site = extractBefore(columns{i*2-1},"_");
        mData = pairElev.(columns{i*2});
        vData = pairElev.(columns{i*2-1});
        biasTable.(site) = mean(mData,"omitnan")-mean(vData,"omitnan");
    end
end
%% Plot function
function plotTileComparison(pairTable, siteNames, yLabel, timeRange, bias, saveDir)
    fig = figure('Name', variableLabel + " Comparison", 'Position', [100, 100, 1200, 600]);
    tiledlayout('flow');
    nexttile;
    hold on;
    scatter(vStruct.(fieldTime), vStruct.(fieldVal),1,...
        'b','filled','DisplayName', 'Validation');
    if contains(variableLabel,'Elevation')
        temp = mStruct.(fieldVal)-bias;
    else
        temp = mStruct.(fieldVal);
    end
    scatter(mStruct.(fieldTime), temp,1,...
        'r','filled','DisplayName','Model');
    title(site + " " + variableLabel);
    xlabel("Time");
    ylabel(yLabel);
    legend('Location', 'best');
    grid on;
    if ~isempty(timeRange)
        xlim(timeRange);
    end
    filename = fullfile(saveDir, "Comparison_"+variableLabel+".pdf");
    exportgraphics(fig, filename, 'ContentType', 'vector');
end
%% Kling-Gupta
function [kge,r,alpha,beta] = klinggupta(model,validation)
model(isnan(validation)) = NaN;
varin  = [model,validation];
mStd   = std(model,"omitnan");
vStd   = std(validation,"omitnan");
mMean  = mean(modelled,"omitnan");
vMean  = mean(validation,"omitnan");
r      = corrcoef(varin,'rows','pairwise'); 
r      = r(1,2);
alpha  = mStd/vStd;   % variability of prediction errors
beta   = mMean/vMean; % bias
%KGE timeseries 
kge    = 1-sqrt(((r-1)^2) + ((alpha-1)^2) + ((beta-1)^2));
end
