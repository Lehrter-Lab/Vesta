% Parameters
inputDir = './outputs/';
varList = {'temperature'};
filePatterns = {'temperature_*.nc'};
nodeList = [13488, 42773, 48505, 80345];  % Pull 4 nodes
timeVar = 'time';
nodeCount = length(nodeList);

% Build file maps for each variable
fileMaps = dictionary;
fileSuffixes = dictionary;
for v = 1:length(varList)
    var = varList{v};
    files = dir(fullfile(inputDir, filePatterns{v}));
    names = {files.name};
    suffixes = regexp(names, '(?<=_)\d+', 'match');
    % Sort suffixes numerically
    numericSuffixes = cellfun(@(x) str2double(x{1}), suffixes);
    [~, sortedIdx] = sort(numericSuffixes);
    sortedSuffixes = suffixes(sortedIdx);
    % Build file map with sorted suffixes
    suffixMap = dictionary;
    for i = 1:length(sortedSuffixes)
        suffixMap(sortedSuffixes{i}) = fullfile(inputDir, names{i});
    end
    fileMaps(var) = {suffixMap};
    fileSuffixes(var) = {sortedSuffixes};
end

% Access suffix list for master variable
varSizes = cellfun(@numel, fileSuffixes.values);
assert(all(varSizes == varSizes(1)), 'Not all cell arrays have the same size.');
masterVar = varList{1};
masterSuffixes = fileSuffixes(masterVar); masterSuffixes = masterSuffixes{1};

% Initialize output files
fids = dictionary;
for i = 1:nodeCount
    nodeID = nodeList(i);
    outputFilename = sprintf('node_%d_%s.tsv', nodeID, strjoin(varList, '_'));
    fid = fopen(outputFilename, 'w');
    fprintf(fid, 'node\tvgrid_layer\ttime');
    for v = 1:length(varList)
        fprintf(fid, '\t%s', varList{v});
    end
    fprintf(fid, '\n');
    fids(nodeList(i)) = fid;
end

%% Loop through files by time step (suffix)
for s = 1:length(masterSuffixes)
    suffix = masterSuffixes{s};
    masterMap = fileMaps(masterVar); masterMap = masterMap{1};
    masterFile = char(masterMap(suffix));
    timeVals = ncread(masterFile, timeVar);
    numTime = length(timeVals);
    if strcmp(masterVar, 'elevation')
        numLayers = 1;
    else
        info = ncinfo(masterFile, masterVar);
        numLayers = info.Size(1);
    end
    varDataAll = dictionary;
    % Precompute node index mapping once & allocate
    [~, nodePos] = ismember(nodeList, nodeList);
    fullBlock = NaN(numLayers, nodeCount, numTime);
    %% Computationally slow part 
    for v = 1:length(varList)
        var = varList{v};
        varMap = fileMaps(var); varMap = varMap{1};
    
        if isKey(varMap, suffix)
            varFile = char(varMap(suffix));
            % Fetch file info
            fileInfo = ncinfo(varFile, var);
            allNodes = max(fileInfo.Size);
            % Only process valid nodes that exist in the file
            validNodes = nodeList <= allNodes;
            targetNodeIndices = nodeList(validNodes);
            % Get the required nodes one at a time and stuff in dataBlock
            for numNodes = 1:length(targetNodeIndices)
                try
                    if strcmp(var, 'elevation')
                        rawData = ncread(varFile, var, [targetNodeIndices(numNodes),1],[1,numTime]);
                        dataBlock = reshape(rawData, [1, 1, numTime]);
                    else
                        dataBlock = ncread(varFile, var, [1, targetNodeIndices(numNodes), 1],[numLayers,1,numTime]);
                    end
                    fullBlock(:, numNodes, :) = dataBlock;
                catch e
                    % Handle any ncread errors
                    fprintf('Error reading %s for time step %s: %s\n', var, char(suffix), e.message);
                end
            end
            varDataAll(var) = {fullBlock};
        else
            % If variable file is missing, fill with NaN
            varDataAll(var) = {NaN(numLayers, nodeCount, numTime)};
            fprintf('Missing %s file for time step %s — filled with NaN\n', var, char(suffix));
        end
    end
    %% Write each node's data
    for n = 1:nodeCount
        nodeID = nodeList(n);
        lines = strings(numTime * numLayers, 1);  % Preallocate buffer
        lineIdx = 1;
        for t = 1:numTime
            for l = 1:numLayers
                parts = sprintf('%d\t%d\t%.6f', nodeID, l, timeVals(t));
                for v = 1:length(varList)
                    block = varDataAll(varList{v}); block = block{1};
                    parts = [parts, sprintf('\t%.6f', block(l, n, t))];
                end
                lines(lineIdx) = parts;
                lineIdx = lineIdx + 1;
            end
        end
    
        fid = fids(nodeID);
        fprintf(fid, '%s\n', lines{:});  % Write entire buffer in one call
    end
    disp(string(masterSuffixes{s}));
end

% Close files
for i = 1:nodeCount
    fclose(fids(nodeList(i)));
end
