%% initialization
function Mesh = parseMSH(filename)
    fid = fopen(filename, 'rb', 'l');
    if fid == -1, error('파일을 열 수 없습니다.'); end
    
    Mesh = struct();
    Mesh.SkippedSections = {}; 
    Mesh.Version = 0;
    Mesh.isBinary = false;
    Mesh.sizeT = 'uint64'; % MSH 4.1에서 사용될 size_t 타입 (기본 8바이트)
    
    % V4를 위한 Entity -> Physical Tag 매핑 컨테이너
    Mesh.EntityToPhysMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
    Mesh.EntityToPartitionMap=containers.Map('KeyType','char','ValueType','double');

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
                    dataSize = v_info(3); % 4.1에서는 sizeof(size_t)를 의미
                    
                    if dataSize == 4, Mesh.sizeT = 'uint32'; end
                    
                    if Mesh.Version ~= 2.2 && Mesh.Version < 4.0
                        fclose(fid); error('지원하지 않는 버전입니다: %.1f', Mesh.Version);
                    end
                    
                    if Mesh.isBinary
                        posBeforeOne = ftell(fid); 
                        one = fread(fid, 1, 'int32'); 
                        
                        if one ~= 1
                            % 애초에 'l'로 열었는데 1이 아니라면, 이 파일은 무조건 'b' (빅 엔디안)입니다.
                            if one == 0x01000000
                                disp('빅 엔디안(Big-Endian) 파일 감지됨: 모드를 전환합니다.');
                                fclose(fid);
                                
                                % 묻지도 따지지도 않고 'b' 옵션으로 다시 열기
                                fid = fopen(filename, 'rb', 'b');
                                fseek(fid, posBeforeOne, 'bof');
                                
                                one_retry = fread(fid, 1, 'int32');
                                if one_retry ~= 1
                                    fclose(fid); error('엔디안 자동 복구에 실패했습니다.');
                                end
                            else
                                fclose(fid); error('손상된 파일입니다. (값: %d)', one);
                            end
                        end
                        fgetl(fid); % 줄바꿈 소모
                    end
                    skipToEnd(fid, '$EndMeshFormat');

                case '$PhysicalNames'
                    Mesh = parsePhysicalNames(fid, Mesh);
                    skipToEnd(fid, '$EndPhysicalNames');
                    
                case '$Entities'
                    if Mesh.Version >= 4.0
                        Mesh = parseEntities(fid, Mesh);
                    end
                    skipToEnd(fid, '$EndEntities');
                case '$PartitionedEntities'
                    if Mesh.Version>=4.0
                        Mesh=parsePartitionedEntities(fid,Mesh);
                    end
                    skipToEnd(fid,'$EndPartitionedEntities');

                case '$Nodes'
                    if Mesh.Version >= 4.0
                        Mesh = parseNodesV4(fid, Mesh);
                    else
                        Mesh = parseNodes(fid, Mesh, Mesh.isBinary);
                    end
                    skipToEnd(fid, '$EndNodes');

                case '$Elements'
                    if Mesh.Version >= 4.0
                        Mesh = parseElementsV4(fid, Mesh);
                    else
                        Mesh = parseElements(fid, Mesh, Mesh.isBinary);
                    end
                    skipToEnd(fid, '$EndElements');

                otherwise
                    sectionName = header(2:end);
                    Mesh.SkippedSections{end+1} = sectionName;
                    skipToEnd(fid, ['$End', sectionName]);
            end
        end
    end

    fclose(fid);

    % 임시 매핑 데이터 제거
    if isfield(Mesh, 'EntityToPhysMap')
        Mesh = rmfield(Mesh, 'EntityToPhysMap');
    end
    
    % =========================================================================
    % 🚀 [궁극의 최적화: 중복 제거, 다중 태그 매핑, 데이터 다이어트]
    % =========================================================================
    if isfield(Mesh, 'Phys') && isfield(Mesh, 'Mat')
        Mesh.Region = struct();
        Mesh.ElementPartitions=cell(length(Mesh.Elements),1);
        temp_ele_phys=cell(length(Mesh.Elements),1);

        physNames = Mesh.Phys.NameToID.keys();
        typeToName = containers.Map({1, 2, 3, 4, 5, 15}, {'Line', 'Tri', 'Quad', 'Tet', 'Hex', 'Point'});
        % ---------------------------------------------------------------------
        % 단계 0. 노드 데이터(Mesh.Nodes) 최적화: 1열(ID) 제거, 행=노드번호
        % ---------------------------------------------------------------------
        if isfield(Mesh, 'Nodes') && ~isempty(Mesh.Nodes)
            % 1. 노드 ID(1열) 기준으로 고유값 인덱스 추출 (자동으로 오름차순 정렬됨)
            [~, unique_node_idx] = unique(Mesh.Nodes(:, 1), 'sorted');
            unique_nodes = Mesh.Nodes(unique_node_idx, :);
            
            % 2. 1행이 1번 노드, 2행이 2번 노드가 되도록 보장 (비어있는 노드 ID 방어)
            max_node_id = unique_nodes(end, 1);
            optimized_nodes = zeros(max_node_id, 3); % X, Y, Z (3열)
            
            % 3. 고유하고 정렬된 좌표 데이터를 해당 노드 번호 행에 직접 꽂아넣기
            optimized_nodes(unique_nodes(:, 1), :) = unique_nodes(:, 2:4);
            
            % 4. 기존의 뚱뚱한 Nodes 배열을 '순수 좌표 배열'로 영구 교체!
            Mesh.Nodes = optimized_nodes;
        end
        % ---------------------------------------------------------------------
        % 단계 1. 각 요소 타입별로 중복 제거 및 데이터 다이어트
        % ---------------------------------------------------------------------
        CleanElements = cell(length(Mesh.Elements), 1);
        RawToUniqueMaps = cell(length(Mesh.Elements), 1);
        
        for eType = 1:length(Mesh.Elements)
            if ~isempty(Mesh.Elements{eType})
                raw_data = Mesh.Elements{eType};
                elm_ids = raw_data(:, 1);
                
                % 고유 요소 번호를 기준으로 중복 제거 (stable로 파일 순서 유지)
                [~, unique_idx, raw_to_unique_map] = unique(elm_ids, 'stable');
                RawToUniqueMaps{eType} = raw_to_unique_map;
                
                % ✂️ [핵심: 2열과 3열 잘라내기]
                % raw_data 구조: [1:ID, 2:PhysTag, 3:ElemTag, 4:PartTag, 5:Node1, 6:Node2 ...]
                % 여기서 2열과 3열을 건너뛰고, 1열과 4열부터 끝까지만 가져옵니다!
                CleanElements{eType} = raw_data(unique_idx, [1, 5:end]);

                part_col=raw_data(unique_idx,4);
                if any(part_col>0)
                    Mesh.ElementPartitions{eType}=[raw_data(unique_idx,1),part_col];
                end
                phys_col=raw_data(unique_idx,2);
                if any(phys_col)
                    temp_ele_phys{eType}=[raw_data(unique_idx,1),phys_col];
                end
            end
        end
        
        valid_phys_cells=temp_ele_phys(~cellfun('isempty',temp_ele_phys));
        if ~isempty(valid_phys_cells)
            Mesh.ele_physTag=vertcat(valid_phys_cells{:});
            Mesh.ele_physTag=sortrows(Mesh.ele_physTag,1);
        else
            Mesh.ele_physTag=[];
        end
        % ---------------------------------------------------------------------
        % 단계 2. Region 정보 매핑 (순수 인덱스 기준)
        % ---------------------------------------------------------------------
        for i = 1:length(physNames)
            pName = physNames{i};
            pID = Mesh.Phys.NameToID(pName);
            vName = matlab.lang.makeValidName(pName); 
            
            foundTypes = [];
            tempIdxMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
            
            for eType = 1:length(Mesh.Elements)
                if ~isempty(Mesh.Elements{eType})
                    raw_data = Mesh.Elements{eType}; % 원본 뚱뚱한 데이터에서 검색
                    matched_raw_idx = find(raw_data(:, 2) == pID);
                    
                    if ~isempty(matched_raw_idx)
                        % [마법의 순간 ✨] 원본 행 번호를 -> 순수 행 번호로 일괄 번역
                        true_idx = RawToUniqueMaps{eType}(matched_raw_idx);
                        
                        % 중복 방지 (안전장치)
                        true_idx = unique(true_idx, 'stable'); 
                        
                        foundTypes(end+1) = eType;
                        tempIdxMap(eType) = true_idx;
                    end
                end
            end
            
            % -----------------------------------------------------------------
            % 단계 3. 직관적 할당 (Auto-Unwrap)
            % -----------------------------------------------------------------
            if isempty(foundTypes)
                Mesh.Region.(vName) = [];
            elseif length(foundTypes) == 1
                Mesh.Region.(vName) = tempIdxMap(foundTypes(1));
            else
                for j = 1:length(foundTypes)
                    eT = foundTypes(j);
                    if isKey(typeToName, eT)
                        fName = typeToName(eT);
                    else
                        fName = sprintf('Type_%d', eT);
                    end
                    Mesh.Region.(vName).(fName) = tempIdxMap(eT);
                end
            end
        end
        
        % ---------------------------------------------------------------------
        % 단계 4. 기존의 무거운 행렬을 '순수 행렬'로 영구 교체!
        % ---------------------------------------------------------------------
        Mesh.Elements = CleanElements;
    end
end % parseMSH 함수 종료

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

%% ParseNodesV2
function Mesh = parseNodes(fid, Mesh, isBinary)
    % 핵심 수정: 바이너리 모드여도 노드 개수는 ASCII로 읽어야 합니다!
    numNodes = fscanf(fid, '%f', 1); % 공백, 엔터 모두 알아서 건너뛰고 숫자 1개만 추출
    fgetl(fid); % 숫자 뒤에 남은 개행문자 소모
    
    if ~isBinary
        % --- ASCII 모드 ---
        data = fscanf(fid, '%f', [4, numNodes])';
        Mesh.Nodes = data;
    else
        % --- Binary 모드 ---
        % 노드 1개당 정확히 28바이트 (ID: 4바이트, X/Y/Z: 각 8바이트)
        % 통째로 바이트(uint8)로 읽어들입니다.
        rawData = fread(fid, [28, numNodes], '*uint8');
        
        % 1~4 바이트: Node ID (int32를 double로 캐스팅)
        ids = double(typecast(reshape(rawData(1:4, :), [], 1), 'int32'));
        
        % 5~28 바이트: X, Y, Z 좌표 (double)
        coords = typecast(reshape(rawData(5:28, :), [], 1), 'double');
        coords = reshape(coords, 3, numNodes)';
        
        % 최종 노드 배열 저장
        Mesh.Nodes = [ids, coords];
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

%% ParseElementsV2
function Mesh = parseElements(fid, Mesh, isBinary)
    Mesh.Elements = cell(19, 1); 
    
    numElementsStr = fgetl(fid);
    numElementsTotal = str2double(strtrim(numElementsStr));
    
    tempBlocks = cell(19, 1); 
    
    % ---------------------------------------------------------
    % 1. 시스템 메모리 기반 동적 임계값 설정
    % ---------------------------------------------------------
    if ispc % Windows 환경인 경우
        [~, sysView] = memory;
        available_RAM = sysView.PhysicalMemory.Available; % 현재 사용 가능한 남은 RAM (바이트)
        
        % 여유 공간의 25%만 파일 로드 임계값으로 설정 (안전 여유분 확보)
        threshold_bytes = available_RAM * 0.25; 
    else
        % Mac/Linux 등 memory 함수 지원이 안 되는 환경의 기본값
        threshold_bytes = 1024 * 1024 * 1024; % 1GB
    end
    
    currentPos = ftell(fid);      
    fseek(fid, 0, 'eof');         
    fileSize_bytes = ftell(fid);  
    fseek(fid, currentPos, 'bof');
    
    loadToMemory = (fileSize_bytes <= threshold_bytes); 
    % ---------------------------------------------------------

    if ~isBinary
        % =========================================================
        % ASCII 모드
        % =========================================================
        if loadToMemory
            % [온메모리] 파일 전체를 한 번에 문자열로 읽고 메모리에서 분리
            rawText = fread(fid, inf, '*char')';
            lines = splitlines(rawText);
            numLines = length(lines);
            lineIdx = 1;
            
            for i = 1:numElementsTotal
                if lineIdx > numLines, break; end
                
                line = strtrim(lines{lineIdx});
                lineIdx = lineIdx + 1;
                
                % 빈 줄 건너뛰기
                while isempty(line) && lineIdx <= numLines
                    line = strtrim(lines{lineIdx});
                    lineIdx = lineIdx + 1;
                end
                if isempty(line), break; end 
                
                val = sscanf(line, '%f')';
                
                elmNum = val(1); elmType = val(2); numTags = val(3);
                
                physTag = 0; elemTag = 0; partTag=0;
                if numTags >= 1, physTag = val(4); end
                if numTags >= 2, elemTag = val(5); end
                if numTags >=4, partTag=val(7);end
                
                nodes = val(3 + numTags + 1 : end);
                tempBlocks{elmType}{end+1} = [elmNum, physTag, elemTag,partTag, nodes];
            end
        else
            % [디스크 I/O] 기존 한 줄씩 읽기 (1GB 초과)
            for i = 1:numElementsTotal
                line = strtrim(fgetl(fid));
                while isempty(line) && ~feof(fid)
                    line = strtrim(fgetl(fid));
                end
                if isempty(line), break; end 
                
                val = sscanf(line, '%f')';
                
                elmNum = val(1); elmType = val(2); numTags = val(3);
                
                physTag = 0; elemTag = 0; partTag=0;
                if numTags >= 1, physTag = val(4); end
                if numTags >= 2, elemTag = val(5); end
                if numTags>=4, partTag=val(7);end
                
                nodes = val(3 + numTags + 1 : end);
                tempBlocks{elmType}{end+1} = [elmNum, physTag, elemTag,partTag, nodes];
            end
        end
        
    else
        % =========================================================
        % Binary 모드
        % =========================================================
        if loadToMemory
            % =========================================================
            % [온메모리] 바이너리 데이터 고속 파싱 (루프 오버헤드 극소화)
            % =========================================================
            allData = fread(fid, inf, '*int32'); 
            totalLen = length(allData);
            
            % 💡 [최적화 1] getNodesCount를 루프 밖에서 미리 배열로 매핑
            % (Gmsh 요소 타입이 대략 100번대까지 있다고 가정하여 넉넉히 150 잡음)
            nodeCountMap = zeros(1, 150);
            for i = 1:150
                try 
                    % 미리 모든 타입의 노드 개수를 계산해 배열에 저장해 둠
                    nodeCountMap(i) = getNodesCount(i); 
                catch
                    % 매핑되지 않은 타입은 무시
                end
            end
            
            % 💡 [최적화 2] tempBlocks 셀 배열 넉넉하게 사전 할당 (Pre-allocation)
            tempBlocks = cell(19, 1);
            tempIdx = ones(19, 1); % 각 요소 타입별 인덱스 추적기
            for i = 1:19
                % 최대 numElementsTotal만큼 들어올 수 있다고 가정하고 미리 방을 만듦
                tempBlocks{i} = cell(numElementsTotal, 1); 
            end
            
            currCount = 0;
            idx = 1; 
            
            while currCount < numElementsTotal && idx <= totalLen
                elmType    = allData(idx);
                numInBlock = allData(idx+1);
                numTags    = allData(idx+2);
                idx = idx + 3; 
                
                % [적용 1] 함수 호출 대신 인덱스 참조 (속도 0초에 수렴)
                numNodes = nodeCountMap(elmType);
                
                blockSize = 1 + numTags + numNodes; 
                blockDataLen = blockSize * numInBlock;
                
                if idx + blockDataLen - 1 > totalLen, break; end
                
                rawData1D = allData(idx : idx + blockDataLen - 1);
                idx = idx + blockDataLen; 
                
                blockMatrix = double(reshape(rawData1D, blockSize, numInBlock))';
                
                elmNum_col  = blockMatrix(:, 1);
                
                % 💡 [최적화 3] zeros() 함수 호출 대신 기존 행렬에 0 곱하기
                physTag_col = elmNum_col * 0; 
                elemTag_col = elmNum_col * 0; 
                partTag_col = elmNum_col * 0;
                if numTags >= 1, physTag_col = blockMatrix(:, 2); end
                if numTags >= 2, elemTag_col = blockMatrix(:, 3); end
                if numTags>=4, partTag_col=blockMatrix(:,5);end
                nodes_cols = blockMatrix(:, 1 + numTags + 1 : end);
                newData = [elmNum_col, physTag_col, elemTag_col,partTag_col, nodes_cols];
                
                % [적용 2] end+1 대신 미리 만들어둔 빈방에 데이터 꽂아넣기
                tempBlocks{elmType}{tempIdx(elmType)} = newData;
                tempIdx(elmType) = tempIdx(elmType) + 1; % 포인터 증가
                
                currCount = currCount + numInBlock;
            end
            
            % 마무리: 임시 블록들을 하나로 병합 (빈 공간 잘라내기)
            for i = 1:19
                if tempIdx(i) > 1
                    % 채워진 부분까지만 추출해서 병합
                    validCells = tempBlocks{i}(1 : tempIdx(i)-1); 
                    Mesh.Elements{i} = vertcat(validCells{:});
                end
            end
        else
            % [디스크 I/O] 기존 블록 단위 파싱 (1GB 초과)
            currCount = 0;
            
            while currCount < numElementsTotal
                header = fread(fid, 3, 'int32');
                if isempty(header) || length(header) < 3, break; end 
                
                elmType = header(1);
                numInBlock = header(2);
                numTags = header(3);
                
                numNodes = getNodesCount(elmType);
                blockSize = 1 + numTags + numNodes; 
                
                rawData = fread(fid, [blockSize, numInBlock], '*int32')';
                blockData = double(rawData);
                
                elmNum_col = blockData(:, 1);
                physTag_col = zeros(numInBlock, 1);
                elemTag_col = zeros(numInBlock, 1);
                
                if numTags >= 1, physTag_col = blockData(:, 2); end
                if numTags >= 2, elemTag_col = blockData(:, 3); end
                
                nodes_cols = blockData(:, 1 + numTags + 1 : end);
                
                newData = [elmNum_col, physTag_col, elemTag_col, nodes_cols];
                tempBlocks{elmType}{end+1} = newData;
                
                currCount = currCount + numInBlock;
            end
        end
    end
    
    % 마무리: 임시 블록들을 하나로 병합
    for i = 1:19
        if ~isempty(tempBlocks{i})
            Mesh.Elements{i} = vertcat(tempBlocks{i}{:});
        end
    end
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

%% GetNodesCount
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