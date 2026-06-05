function Mesh = parseMSHV4(filename)
    fid = fopen(filename, 'rb', 'l');
    if fid == -1, error('파일을 열 수 없습니다.'); end
    
    Mesh = struct();
    Mesh.SkippedSections = {}; 
    Mesh.Version = 0;
    Mesh.isBinary = false;
    Mesh.sizeT = 'uint64'; % MSH 4.1 기본 8바이트
    
    Mesh.EntityToPhysMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
    Mesh.EntityToPartitionMap = containers.Map('KeyType', 'char', 'ValueType', 'double');
    
    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line), break; end
        header = strtrim(line);
        
        if startsWith(header, '$') && ~startsWith(header, '$End')
            switch header
                case '$MeshFormat'
                    fmt = fgetl(fid);
                    v_info = sscanf(fmt, '%f %d %d');
                    Mesh.Version = v_info(1);
                    Mesh.isBinary = (v_info(2) == 1);
                    dataSize = v_info(3);
                    if dataSize == 4, Mesh.sizeT = 'uint32'; end
                    
                    if Mesh.Version < 4.0
                        fclose(fid); error('이 파일은 V2.2 포맷입니다. parseMSHV2를 사용하세요.');
                    end
                    
                    if Mesh.isBinary
                        posBeforeOne = ftell(fid); 
                        one = fread(fid, 1, 'int32'); 
                        if one ~= 1
                            if one == 0x01000000
                                disp('빅 엔디안 감지됨: 모드 전환');
                                fclose(fid); fid = fopen(filename, 'rb', 'b'); fseek(fid, posBeforeOne, 'bof');
                                if fread(fid, 1, 'int32') ~= 1, fclose(fid); error('엔디안 복구 실패'); end
                            else
                                fclose(fid); error('손상된 파일');
                            end
                        end
                        fgetl(fid); 
                    end
                    skipToEnd(fid, '$EndMeshFormat');
                    
                case '$PhysicalNames'
                    Mesh = parsePhysicalNames(fid, Mesh);
                    skipToEnd(fid, '$EndPhysicalNames');
                    
                case '$Entities'
                    Mesh = parseEntities(fid, Mesh);
                    skipToEnd(fid, '$EndEntities');
                    
                case '$PartitionedEntities'
                    Mesh = parsePartitionedEntities(fid, Mesh);
                    skipToEnd(fid, '$EndPartitionedEntities');
                    
                case '$Nodes'
                    Mesh = parseNodesV4(fid, Mesh);
                    skipToEnd(fid, '$EndNodes');
                    
                case '$Elements'
                    Mesh = parseElementsV4(fid, Mesh);
                    skipToEnd(fid, '$EndElements');
                    
                otherwise
                    sectionName = header(2:end);
                    Mesh.SkippedSections{end+1} = sectionName;
                    skipToEnd(fid, ['$End', sectionName]);
            end
        end
    end
    fclose(fid);
    
    if isfield(Mesh, 'EntityToPhysMap'), Mesh = rmfield(Mesh, 'EntityToPhysMap'); end
    
    % =========================================================================
    % 🚀 [궁극의 최적화: 데이터 다이어트 및 Region 매핑] 
    % (V2와 동일한 매커니즘 사용, 생략 없이 V2의 최적화 블록을 그대로 붙여넣습니다)
    % =========================================================================
    % 0. Nodes 최적화 (행 번호 = 노드 ID)
    if isfield(Mesh, 'Nodes') && ~isempty(Mesh.Nodes)
        [~, unique_node_idx] = unique(Mesh.Nodes(:, 1), 'sorted');
        unique_nodes = Mesh.Nodes(unique_node_idx, :);
        max_node_id = unique_nodes(end, 1);
        optimized_nodes = zeros(max_node_id, 3);
        optimized_nodes(unique_nodes(:, 1), :) = unique_nodes(:, 2:4);
        Mesh.Nodes = optimized_nodes;
    end

    % 1. Elements 최적화 및 Region 분리
    if isfield(Mesh, 'Phys') && isfield(Mesh, 'Mat')
        Mesh.Region = struct();
        Mesh.ElementPartitions = cell(length(Mesh.Elements),1);
        temp_ele_phys = cell(length(Mesh.Elements),1);
        physNames = Mesh.Phys.NameToID.keys();
        typeToName = containers.Map({1, 2, 3, 4, 5, 15}, {'Line', 'Tri', 'Quad', 'Tet', 'Hex', 'Point'});
        
        CleanElements = cell(length(Mesh.Elements), 1);
        RawToUniqueMaps = cell(length(Mesh.Elements), 1);
        
        for eType = 1:length(Mesh.Elements)
            if ~isempty(Mesh.Elements{eType})
                raw_data = Mesh.Elements{eType};
                elm_ids = raw_data(:, 1);
                [~, unique_idx, raw_to_unique_map] = unique(elm_ids, 'stable');
                RawToUniqueMaps{eType} = raw_to_unique_map;
                
                CleanElements{eType} = raw_data(unique_idx, [1, 5:end]);
                part_col = raw_data(unique_idx, 4);
                if any(part_col > 0), Mesh.ElementPartitions{eType} = [raw_data(unique_idx, 1), part_col]; end
                phys_col = raw_data(unique_idx, 2);
                if any(phys_col), temp_ele_phys{eType} = [raw_data(unique_idx, 1), phys_col]; end
            end
        end
        
        valid_phys_cells = temp_ele_phys(~cellfun('isempty', temp_ele_phys));
        if ~isempty(valid_phys_cells)
            Mesh.ele_physTag = sortrows(vertcat(valid_phys_cells{:}), 1);
        else
            Mesh.ele_physTag = [];
        end
        
        for i = 1:length(physNames)
            pName = physNames{i}; pID = Mesh.Phys.NameToID(pName); vName = matlab.lang.makeValidName(pName); 
            foundTypes = []; tempIdxMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
            for eType = 1:length(Mesh.Elements)
                if ~isempty(Mesh.Elements{eType})
                    matched_raw_idx = find(Mesh.Elements{eType}(:, 2) == pID);
                    if ~isempty(matched_raw_idx)
                        true_idx = unique(RawToUniqueMaps{eType}(matched_raw_idx), 'stable'); 
                        foundTypes(end+1) = eType; tempIdxMap(eType) = true_idx;
                    end
                end
            end
            if isempty(foundTypes), Mesh.Region.(vName) = [];
            elseif length(foundTypes) == 1, Mesh.Region.(vName) = tempIdxMap(foundTypes(1));
            else
                for j = 1:length(foundTypes)
                    eT = foundTypes(j); if isKey(typeToName, eT), fName = typeToName(eT); else, fName = sprintf('Type_%d', eT); end
                    Mesh.Region.(vName).(fName) = tempIdxMap(eT);
                end
            end
        end
        Mesh.Elements = CleanElements;
    end
end

% =========================================================================

%% ParseEntities
function Mesh = parseEntities(fid, Mesh)
    if ~Mesh.isBinary
        % ---------------------------------------------------------
        % [ASCII 모드]
        % ---------------------------------------------------------
        counts = [];
        while isempty(counts) && ~feof(fid)
            line = strtrim(fgetl(fid));
            counts = sscanf(line, '%lu'); 
        end
        
        if length(counts) < 4
            counts = [counts; zeros(4-length(counts), 1)];
        end

        dims = [0, 1, 2, 3]; % Point, Curve, Surface, Volume
        
        for d = 1:4
            for i = 1:counts(d)
                line_data = sscanf(strtrim(fgetl(fid)), '%f');
                if isempty(line_data), continue; end
                
                tag = line_data(1);
                
                % [💡수정된 부분] 정확한 피지컬 태그 개수(numPhysicalTags)의 위치 지정
                % Point는 [Tag, X, Y, Z] 4개 뒤인 5번째
                % 나머지는 [Tag, minX, minY, minZ, maxX, maxY, maxZ] 7개 뒤인 8번째
                offset = 5; 
                if dims(d) > 0, offset = 8; end
                
                numPhysTags = line_data(offset);
                
                if numPhysTags > 0
                    physTags = line_data(offset+1 : offset+numPhysTags);
                    key = sprintf('%d_%d', dims(d), tag);
                    Mesh.EntityToPhysMap(key) = physTags;
                end
            end
        end
        
    else
        % ---------------------------------------------------------
        % [Binary 모드] (기존과 동일 - 바이너리는 명세대로 정확히 작동함)
        % ---------------------------------------------------------
        counts = double(fread(fid, 4, Mesh.sizeT));
        if isempty(counts), counts = zeros(4,1); end
        
        dims = [0, 1, 2, 3];
        
        for d = 1:4
            for i = 1:counts(d)
                tag = fread(fid, 1, 'int32');
                
                if dims(d) == 0
                    fread(fid, 3, 'double'); 
                else
                    fread(fid, 6, 'double'); 
                end
                
                numPhysTags = fread(fid, 1, Mesh.sizeT);
                if numPhysTags > 0
                    physTags = fread(fid, numPhysTags, 'int32');
                    key = sprintf('%d_%d', dims(d), tag);
                    Mesh.EntityToPhysMap(key) = physTags;
                end
                
                if dims(d) > 0
                    numBounds = fread(fid, 1, Mesh.sizeT);
                    if numBounds > 0
                        fread(fid, numBounds, 'int32'); 
                    end
                end
            end
        end
    end
end
%% parsePartitionedEntities
function Mesh = parsePartitionedEntities(fid, Mesh)
    if ~Mesh.isBinary
        % [ASCII 모드]
        fscanf(fid, '%lu', 1); % numPartitions 스킵
        numGhosts = fscanf(fid, '%lu', 1);
        fgetl(fid); 
        
        for i = 1:numGhosts
            fgetl(fid); % Ghost 정보는 임시 스킵
        end
        
        line = strtrim(fgetl(fid));
        counts = sscanf(line, '%lu');
        if length(counts) < 4, counts = [counts; zeros(4-length(counts), 1)]; end
        dims = [0, 1, 2, 3];
        
        for d = 1:4
            for i = 1:counts(d)
                line = strtrim(fgetl(fid));
                if isempty(line), continue; end
                data = sscanf(line, '%f');
                
                tag = data(1);
                numParts = data(4);
                
                if numParts > 0
                    partTag = data(5); % 첫 번째 할당 파티션 번호
                    key = sprintf('%d_%d', dims(d), tag);
                    Mesh.EntityToPartitionMap(key) = partTag;
                end
            end
        end
    else
        % [Binary 모드]
        fread(fid, 1, Mesh.sizeT); % numPartitions 스킵
        numGhosts = fread(fid, 1, Mesh.sizeT);
        
        if numGhosts > 0
            fread(fid, numGhosts * 2, 'int32');
        end
        
        counts = double(fread(fid, 4, Mesh.sizeT));
        if isempty(counts), counts = zeros(4,1); end
        dims = [0, 1, 2, 3];
        
        for d = 1:4
            for i = 1:counts(d)
                tag_info = fread(fid, 3, 'int32'); 
                tag = tag_info(1);
                
                numParts = fread(fid, 1, Mesh.sizeT);
                if numParts > 0
                    partTags = fread(fid, numParts, 'int32');
                    key = sprintf('%d_%d', dims(d), tag);
                    Mesh.EntityToPartitionMap(key) = partTags(1); 
                end
                
                % Coordinate / Bounds 데이터 스킵 처리
                if dims(d) == 0
                    fread(fid, 3, 'double');
                else
                    fread(fid, 6, 'double');
                end
                
                numPhysTags = fread(fid, 1, Mesh.sizeT);
                if numPhysTags > 0, fread(fid, numPhysTags, 'int32'); end
                
                if dims(d) > 0
                    numBounds = fread(fid, 1, Mesh.sizeT);
                    if numBounds > 0, fread(fid, numBounds, 'int32'); end
                end
            end
        end
    end
end
%% ParseNodeV4 (완벽 수정본)
function Mesh = parseNodesV4(fid, Mesh)
    if ~Mesh.isBinary
        % ---------------------------------------------------------
        % [ASCII 모드] 
        % ---------------------------------------------------------
        header = fscanf(fid, '%f', 4);
        numEntityBlocks = header(1);
        numNodesTotal = header(2);
        
        nodeIDs = zeros(numNodesTotal, 1);
        nodeCoords = zeros(numNodesTotal, 3);
        currIdx = 1;
        
        for b = 1:numEntityBlocks
            blockHeader = fscanf(fid, '%f', 4);
            if isempty(blockHeader), break; end
            
            entityDim = blockHeader(1);
            parametric = blockHeader(3);
            numNodesInBlock = blockHeader(4);
            
            if numNodesInBlock == 0, continue; end
            
            % 노드 ID 배열 읽기
            blockIDs = fscanf(fid, '%f', numNodesInBlock);
            
            % 💡 [수정 1] XYZ 3차원 좌표 (3 * N)만 정확하게 분리해서 읽기
            rawCoords = fscanf(fid, '%f', 3 * numNodesInBlock);
            blockCoords = reshape(rawCoords, 3, [])';
            
            % 💡 [수정 2] 파라매트릭 좌표(U, V, W)가 뒤에 붙어있다면, 따로 읽어서 휴지통에 버리기
            if parametric >= 1 && entityDim > 0
                fscanf(fid, '%f', entityDim * numNodesInBlock);
            end
            
            idxRange = currIdx : currIdx + numNodesInBlock - 1;
            nodeIDs(idxRange) = blockIDs;
            nodeCoords(idxRange, :) = blockCoords;
            
            % 💡 [수정 3] 루프 인덱스 업데이트 누락 복구!
            currIdx = currIdx + numNodesInBlock; 
        end
        
    else
        % ---------------------------------------------------------
        % [Binary 모드] 
        % ---------------------------------------------------------
        header = fread(fid, 4, Mesh.sizeT);
        numEntityBlocks = header(1);
        numNodesTotal = header(2);
        
        nodeIDs = zeros(numNodesTotal, 1);
        nodeCoords = zeros(numNodesTotal, 3);
        currIdx = 1;
        
        for b = 1:numEntityBlocks
            blockHeaderInt = fread(fid, 3, 'int32');
            entityDim = blockHeaderInt(1);
            parametric = blockHeaderInt(3);
            
            numNodesInBlock_raw = fread(fid, 1, Mesh.sizeT);
            if isempty(numNodesInBlock_raw) || numNodesInBlock_raw == 0
                continue; 
            end
            numNodesInBlock = double(numNodesInBlock_raw);
            
            % 노드 ID 배열 읽기
            blockIDs = fread(fid, numNodesInBlock, Mesh.sizeT);
            
            % 💡 [수정 1] XYZ 좌표만 한 번에 정확히 읽어오기 (totalCols 조합 삭제)
            rawCoords = fread(fid, [3, numNodesInBlock], 'double')';
            blockCoords = rawCoords(:, 1:3);
            
            % 💡 [수정 2] 파라매트릭 좌표가 있으면 fseek으로 C++처럼 건너뛰기 (double은 8바이트)
            if parametric >= 1 && entityDim > 0
                fseek(fid, entityDim * numNodesInBlock * 8, 'cof'); 
            end
            
            idxRange = currIdx : currIdx + numNodesInBlock - 1;
            nodeIDs(idxRange) = blockIDs;
            nodeCoords(idxRange, :) = blockCoords;
            
            currIdx = currIdx + numNodesInBlock;
        end
    end
    
    % 최종 결과 통합: [ID, X, Y, Z]
    Mesh.Nodes = [nodeIDs, nodeCoords];
end

%% ParseElementsV4
function Mesh = parseElementsV4(fid, Mesh)
    Mesh.Elements = cell(93, 1); 
    tempBlocks = cell(93, 1);
    
    if ~Mesh.isBinary
        % ---------------------------------------------------------
        % [ASCII 모드] 혼용 버그를 막기 위해 전부 fscanf로 읽기
        % ---------------------------------------------------------
        header = fscanf(fid, '%f', 4);
        numEntityBlocks = header(1);
        
        for b = 1:numEntityBlocks
            blockHeader = fscanf(fid, '%f', 4);
            
            if isempty(blockHeader), break; end
            
            entityDim = blockHeader(1);
            entityTag = blockHeader(2);
            elementType = blockHeader(3);
            numElementsInBlock = blockHeader(4);
            
            if numElementsInBlock == 0, continue; end
            
            % ---------------------------------------------------------
            % [수정됨] Physical Tag 매핑 (다중 태그 지원)
            % ---------------------------------------------------------
            key = sprintf('%d_%d', entityDim, entityTag);
            if isKey(Mesh.EntityToPhysMap, key)
                tags = Mesh.EntityToPhysMap(key);
            else
                tags = 0; % 물리 태그가 아예 없는 경우 기본값 0
            end

            if isKey(Mesh.EntityToPartitionMap,key)
                partTag=Mesh.EntityToPartitionMap(key);
            else
                partTag=0;
            end
            
            numNodesPerElem = getNodesCount(elementType);
            
            % 데이터 읽기 (파일 I/O는 무조건 한 번만!)
            rawVals = fscanf(fid, '%f', (1 + numNodesPerElem) * numElementsInBlock);
            blockData = reshape(rawVals, 1 + numNodesPerElem, [])';
            
            elmNum_col = blockData(:, 1);
            nodes_cols = blockData(:, 2:end);
            elemTag_col = repmat(entityTag, numElementsInBlock, 1);
            partTag_col=repmat(partTag,numElementsInBlock,1);
            % 태그 개수만큼 반복하며 요소를 복제하여 할당
            for t = 1:length(tags)
                physTag_col = repmat(tags(t), numElementsInBlock, 1);
                
                newData = [double(elmNum_col), double(physTag_col), double(elemTag_col),double(partTag_col), double(nodes_cols)];
                tempBlocks{elementType}{end+1} = newData;
            end
        end
    else
        % ---------------------------------------------------------
        % [Binary 모드] 초고속 인메모리 스캔
        % ---------------------------------------------------------
        % 헤더 전체 읽기 (size_t)
        header = fread(fid, 4, Mesh.sizeT);
        numEntityBlocks = header(1);
        
        for b = 1:numEntityBlocks
            % 1. 블록 헤더 앞부분 3개 읽기 (int32)
            blockHeaderInt = fread(fid, 3, 'int32');
            entityDim = blockHeaderInt(1);
            entityTag = blockHeaderInt(2);
            elementType = blockHeaderInt(3);
            
            % 2. 블록 헤더 뒷부분 1개 읽기 (size_t)
            numElementsInBlock = fread(fid, 1, Mesh.sizeT);
            
            if numElementsInBlock == 0, continue; end
            
            % ---------------------------------------------------------
            % [수정됨] Physical Tag 매핑 (다중 태그 지원)
            % ---------------------------------------------------------
            key = sprintf('%d_%d', entityDim, entityTag);
            if isKey(Mesh.EntityToPhysMap, key)
                tags = Mesh.EntityToPhysMap(key);
            else
                tags = 0; % 물리 태그가 아예 없는 경우 기본값 0
            end

            if isKey(Mesh.EntityToPartitionMap,key)
                partTag=Mesh.EntityToPartitionMap(key);
            else
                partTag=0;
            end
            
            numNodesPerElem = getNodesCount(elementType);
            blockSize = 1 + numNodesPerElem; % Element ID + Node IDs
            
            % 3. 요소 데이터 한 번에 읽기 (파일 I/O는 무조건 한 번만!)
            rawData = fread(fid, [blockSize, numElementsInBlock], Mesh.sizeT)';
            
            elmNum_col = rawData(:, 1);
            nodes_cols = rawData(:, 2:end);
            elemTag_col = repmat(entityTag, numElementsInBlock, 1);
            partTag_col=repmat(partTag,numElementsInBlock,1);
            % 태그 개수만큼 반복하며 요소를 복제하여 할당
            for t = 1:length(tags)
                physTag_col = repmat(tags(t), numElementsInBlock, 1);
                
                % V2.2와 동일한 포맷으로 병합 (전부 double로 통일)
                newData = [double(elmNum_col), double(physTag_col), double(elemTag_col),double(partTag_col), double(nodes_cols)];
                tempBlocks{elementType}{end+1} = newData;
            end
        end
    end
    
    % 마무리: 임시 블록들을 하나로 병합하여 구조체에 저장
    for i = 1:93
        if ~isempty(tempBlocks{i})
            Mesh.Elements{i} = vertcat(tempBlocks{i}{:});
        end
    end
end

function n = getNodesCount(elmType)
    % Gmsh MSH 요소 타입별 노드 개수 (V2.2 & V4.1 완벽 호환)
    % 인덱스가 곧 elmType이 되도록 100번까지 넉넉히 할당
    counts = zeros(100, 1);
    
    % =========================================================================
    % 1~19번 (기존 V2.2 및 기본 요소)
    % =========================================================================
    counts(1)  = 2;   % 1: 2-node line.
    counts(2)  = 3;   % 2: 3-node triangle.
    counts(3)  = 4;   % 3: 4-node quadrangle.
    counts(4)  = 4;   % 4: 4-node tetrahedron.
    counts(5)  = 8;   % 5: 8-node hexahedron.
    counts(6)  = 6;   % 6: 6-node prism.
    counts(7)  = 5;   % 7: 5-node pyramid.
    counts(8)  = 3;   % 8: 3-node second order line (2 nodes associated with the vertices and 1 with the edge).
    counts(9)  = 6;   % 9: 6-node second order triangle (3 nodes associated with the vertices and 3 with the edges).
    counts(10) = 9;   % 10: 9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face).
    counts(11) = 10;  % 11: 10-node second order tetrahedron (4 nodes associated with the vertices and 6 with the edges).
    counts(12) = 27;  % 12: 27-node second order hexahedron (8 nodes associated with the vertices, 12 with the edges, 6 with the faces and 1 with the volume).
    counts(13) = 18;  % 13: 18-node second order prism (6 nodes associated with the vertices, 9 with the edges and 3 with the quadrangular faces).
    counts(14) = 14;  % 14: 14-node second order pyramid (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face).
    counts(15) = 1;   % 15: 1-node point.
    counts(16) = 8;   % 16: 8-node second order quadrangle (4 nodes associated with the vertices and 4 with the edges).
    counts(17) = 20;  % 17: 20-node second order hexahedron (8 nodes associated with the vertices and 12 with the edges).
    counts(18) = 15;  % 18: 15-node second order prism (6 nodes associated with the vertices and 9 with the edges).
    counts(19) = 13;  % 19: 13-node second order pyramid (5 nodes associated with the vertices and 8 with the edges).

    % =========================================================================
    % 20~31번 (V4.1 고차 요소 - High-order elements)
    % =========================================================================
    counts(20) = 9;   % 20: 9-node third order incomplete triangle (3 nodes associated with the vertices, 6 with the edges)
    counts(21) = 10;  % 21: 10-node third order triangle (3 nodes associated with the vertices, 6 with the edges, 1 with the face)
    counts(22) = 12;  % 22: 12-node fourth order incomplete triangle (3 nodes associated with the vertices, 9 with the edges)
    counts(23) = 15;  % 23: 15-node fourth order triangle (3 nodes associated with the vertices, 9 with the edges, 3 with the face)
    counts(24) = 15;  % 24: 15-node fifth order incomplete triangle (3 nodes associated with the vertices, 12 with the edges)
    counts(25) = 21;  % 25: 21-node fifth order complete triangle (3 nodes associated with the vertices, 12 with the edges, 6 with the face)
    counts(26) = 4;   % 26: 4-node third order edge (2 nodes associated with the vertices, 2 internal to the edge)
    counts(27) = 5;   % 27: 5-node fourth order edge (2 nodes associated with the vertices, 3 internal to the edge)
    counts(28) = 6;   % 28: 6-node fifth order edge (2 nodes associated with the vertices, 4 internal to the edge)
    counts(29) = 20;  % 29: 20-node third order tetrahedron (4 nodes associated with the vertices, 12 with the edges, 4 with the faces)
    counts(30) = 35;  % 30: 35-node fourth order tetrahedron (4 nodes associated with the vertices, 18 with the edges, 12 with the faces, 1 in the volume)
    counts(31) = 56;  % 31: 56-node fifth order tetrahedron (4 nodes associated with the vertices, 24 with the edges, 24 with the faces, 4 in the volume)

    % =========================================================================
    % 92~93번 (V4.1 고차 헥사헤드론)
    % =========================================================================
    counts(92) = 64;  % 92: 64-node third order hexahedron (8 nodes associated with the vertices, 24 with the edges, 24 with the faces, 8 in the volume)
    counts(93) = 125; % 93: 125-node fourth order hexahedron (8 nodes associated with the vertices, 36 with the edges, 54 with the faces, 27 in the volume)

    % 유효성 검사
    if elmType < 1 || elmType > length(counts) || counts(elmType) == 0
        error('지원하지 않거나 알 수 없는 요소 타입입니다: %d', elmType);
    end
    
    n = counts(elmType);
end


%% SkipToEnd
function skipToEnd(fid, endTag)
    chunkSize = 5242880; % 5MB 청크 단위 (빠른 I/O)
    
    tagBytes = uint8(endTag);
    tagLen = length(tagBytes);
    
    newlineByte = uint8(10); % \n (LF)
    crByte      = uint8(13); % \r (CR)
    spaceByte   = uint8(32); % 공백 문자
    
    while ~feof(fid)
        startPos = ftell(fid);
        % 인코딩 없이 순수 바이트로 초고속 읽기
        chunk = fread(fid, chunkSize, '*uint8')';
        
        idx = strfind(chunk, tagBytes);
        
        % 태그를 찾았다면, 그것이 '진짜'인지 검증합니다.
        for i = 1:length(idx)
            matchIdx = idx(i);
            
            % ----------------------------------------------------
            % [검증 1] 태그 앞부분 검사 (새로운 줄에서 시작했는가?)
            % ----------------------------------------------------
            isValidStart = false;
            if matchIdx == 1
                % 청크의 맨 첫 글자에 걸린 경우, 파일 포인터를 뒤로 돌려 이전 글자 확인
                if startPos == 0
                    isValidStart = true; % 파일의 맨 처음이면 인정
                else
                    fseek(fid, startPos - 1, 'bof');
                    prevChar = fread(fid, 1, '*uint8');
                    if prevChar == newlineByte || prevChar == crByte
                        isValidStart = true;
                    end
                    fseek(fid, startPos, 'bof'); % 원상 복구
                end
            else
                % 청크 중간에 있다면 바로 앞 인덱스 확인
                prevChar = chunk(matchIdx - 1);
                if prevChar == newlineByte || prevChar == crByte
                    isValidStart = true;
                end
            end
            
            % ----------------------------------------------------
            % [검증 2] 태그 뒷부분 검사 (태그 직후에 개행이나 공백이 오는가?)
            % ----------------------------------------------------
            isValidEnd = false;
            endIdx = matchIdx + tagLen; 
            
            if endIdx <= length(chunk)
                nextChar = chunk(endIdx);
                % 태그 뒤에 이상한 바이너리 값이 아닌 제어 문자가 오는지 확인
                if nextChar == newlineByte || nextChar == crByte || nextChar == spaceByte
                    isValidEnd = true;
                end
            else
                % 청크 끝부분에 잘려서 뒤를 알 수 없다면, 
                % 일단 무시하고 다음 루프에서 (fseek 후진 후) 다시 검사함
                continue; 
            end
            
            % ----------------------------------------------------
            % [최종 확인] 앞뒤가 모두 완벽한 진짜 태그라면!
            % ----------------------------------------------------
            if isValidStart && isValidEnd
                % 태그 뒤에 이어지는 진짜 개행 문자(\n) 위치 탐색
                subChunk = chunk(matchIdx:end);
                newlineIdx = strfind(subChunk, newlineByte);
                
                if ~isempty(newlineIdx)
                    % 개행 문자 바로 다음 위치로 파일 포인터 세팅
                    offset = matchIdx - 1 + newlineIdx(1);
                else
                    offset = matchIdx - 1 + tagLen;
                end
                
                fseek(fid, startPos + offset, 'bof');
                return; % 찾았으므로 함수 완벽 종료
            end
        end
        
        % 태그를 못 찾았거나, 찾았는데 바이너리 '가짜' 태그였다면?
        % 청크 경계면 잘림을 대비해 태그 길이만큼 살짝 후진 후 다음 5MB 스캔
        if ~feof(fid)
            fseek(fid, -tagLen, 'cof');
        end
    end
end


%% ParsePhysicalNames
function Mesh = parsePhysicalNames(fid, Mesh)
    % 빈 줄에 대비해 fscanf로 깔끔하게 숫자만 읽어옵니다.
    numNames = fscanf(fid, '%d', 1);
    fgetl(fid); % 숫자 뒤의 줄바꿈 소모
    
    % 양방향 고속 검색 및 차원(dim) 보존을 위한 Map (딕셔너리) 생성
    Mesh.Phys = struct();
    Mesh.Phys.NameToID  = containers.Map('KeyType', 'char', 'ValueType', 'double');
    Mesh.Phys.IDToName  = containers.Map('KeyType', 'double', 'ValueType', 'char');
    Mesh.Phys.NameToDim = containers.Map('KeyType', 'char', 'ValueType', 'double'); % 차원 정보 보존
    
    % 솔버 스크립트에서 직관적으로 쓸 수 있는 매핑 구조체 (자동완성 지원)
    Mesh.Mat = struct(); 
    
    currIdx = 1;
    while currIdx <= numNames
        line = strtrim(fgetl(fid));
        if isempty(line), continue; end
        
        data = textscan(line, '%d %d %q'); 
        
        dim = double(data{1});
        id = double(data{2});
        name = char(data{3});
        
        % 1. Map (딕셔너리)에 데이터 꼼꼼히 저장
        Mesh.Phys.NameToID(name) = id;
        Mesh.Phys.IDToName(id) = name;
        Mesh.Phys.NameToDim(name) = dim;
        
        % 2. 점(.) 연산자로 바로 ID에 접근하기 위한 구조체 저장
        validName = matlab.lang.makeValidName(name);
        Mesh.Mat.(validName) = id;
        
        currIdx = currIdx + 1;
    end
end