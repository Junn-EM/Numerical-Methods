function Mesh = parseMSHV2(filename)
    fid = fopen(filename, 'rb', 'l');
    if fid == -1, error('파일을 열 수 없습니다.'); end
    
    Mesh = struct();
    Mesh.SkippedSections = {}; 
    Mesh.Version = 0;
    Mesh.isBinary = false;
    
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
                    
                    if Mesh.Version >= 4.0
                        fclose(fid); error('이 파일은 V4.1 포맷입니다. parseMSHV4를 사용하세요.');
                    end
                    
                    if Mesh.isBinary
                        posBeforeOne = ftell(fid); 
                        one = fread(fid, 1, 'int32'); 
                        if one ~= 1
                            if one == 0x01000000
                                disp('빅 엔디안(Big-Endian) 감지됨: 모드를 전환합니다.');
                                fclose(fid);
                                fid = fopen(filename, 'rb', 'b');
                                fseek(fid, posBeforeOne, 'bof');
                                if fread(fid, 1, 'int32') ~= 1, fclose(fid); error('엔디안 복구 실패'); end
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
                    
                case '$Nodes'
                    Mesh = parseNodesV2(fid, Mesh, Mesh.isBinary);
                    skipToEnd(fid, '$EndNodes');
                    
                case '$Elements'
                    Mesh = parseElementsV2(fid, Mesh, Mesh.isBinary);
                    skipToEnd(fid, '$EndElements');
                    
                otherwise
                    sectionName = header(2:end);
                    Mesh.SkippedSections{end+1} = sectionName;
                    skipToEnd(fid, ['$End', sectionName]);
            end
        end
    end
    fclose(fid);
    
    % =========================================================================
    % 🚀 [궁극의 최적화: 데이터 다이어트 및 Region 매핑]
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
                
                % [1:ID, 4:PartTag, 5~:Nodes] 만 남김 (2:Phys, 3:Elem 제거)
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
        
        % Region 매핑
        for i = 1:length(physNames)
            pName = physNames{i};
            pID = Mesh.Phys.NameToID(pName);
            vName = matlab.lang.makeValidName(pName); 
            foundTypes = [];
            tempIdxMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
            
            for eType = 1:length(Mesh.Elements)
                if ~isempty(Mesh.Elements{eType})
                    raw_data = Mesh.Elements{eType};
                    matched_raw_idx = find(raw_data(:, 2) == pID);
                    if ~isempty(matched_raw_idx)
                        true_idx = unique(RawToUniqueMaps{eType}(matched_raw_idx), 'stable'); 
                        foundTypes(end+1) = eType;
                        tempIdxMap(eType) = true_idx;
                    end
                end
            end
            
            if isempty(foundTypes)
                Mesh.Region.(vName) = [];
            elseif length(foundTypes) == 1
                Mesh.Region.(vName) = tempIdxMap(foundTypes(1));
            else
                for j = 1:length(foundTypes)
                    eT = foundTypes(j);
                    if isKey(typeToName, eT), fName = typeToName(eT); else, fName = sprintf('Type_%d', eT); end
                    Mesh.Region.(vName).(fName) = tempIdxMap(eT);
                end
            end
        end
        Mesh.Elements = CleanElements;
    end
end

% =========================================================================
% 로컬 헬퍼 함수들 (parseMSHV2 전용)
% =========================================================================
function Mesh = parseNodesV2(fid, Mesh, isBinary)
    numNodes = fscanf(fid, '%f', 1); fgetl(fid); 
    if ~isBinary
        data = fscanf(fid, '%f', [4, numNodes])';
        Mesh.Nodes = data;
    else
        rawData = fread(fid, [28, numNodes], '*uint8');
        ids = double(typecast(reshape(rawData(1:4, :), [], 1), 'int32'));
        coords = reshape(typecast(reshape(rawData(5:28, :), [], 1), 'double'), 3, numNodes)';
        Mesh.Nodes = [ids, coords];
    end
end

function Mesh = parseElementsV2(fid, Mesh, isBinary)
    % (기존 V2 파서 코드와 완벽히 동일하므로 생략 없이 원본 로직 적용)
    Mesh.Elements = cell(19, 1); 
    numElementsStr = fgetl(fid);
    numElementsTotal = str2double(strtrim(numElementsStr));
    tempBlocks = cell(19, 1); 
    
    if ispc
        [~, sysView] = memory; threshold_bytes = sysView.PhysicalMemory.Available * 0.25; 
    else
        threshold_bytes = 1024 * 1024 * 1024;
    end
    
    currentPos = ftell(fid); fseek(fid, 0, 'eof'); fileSize_bytes = ftell(fid); fseek(fid, currentPos, 'bof');
    loadToMemory = (fileSize_bytes <= threshold_bytes); 
    
    if ~isBinary
        if loadToMemory
            rawText = fread(fid, inf, '*char')'; lines = splitlines(rawText);
            numLines = length(lines); lineIdx = 1;
            for i = 1:numElementsTotal
                if lineIdx > numLines, break; end
                line = strtrim(lines{lineIdx}); lineIdx = lineIdx + 1;
                while isempty(line) && lineIdx <= numLines, line = strtrim(lines{lineIdx}); lineIdx = lineIdx + 1; end
                if isempty(line), break; end 
                val = sscanf(line, '%f')';
                elmNum = val(1); elmType = val(2); numTags = val(3);
                physTag = 0; elemTag = 0; partTag=0;
                if numTags >= 1, physTag = val(4); end
                if numTags >= 2, elemTag = val(5); end
                if numTags >= 4, partTag=val(7); end
                nodes = val(3 + numTags + 1 : end);
                tempBlocks{elmType}{end+1} = [elmNum, physTag, elemTag, partTag, nodes];
            end
        else
            for i = 1:numElementsTotal
                line = strtrim(fgetl(fid));
                while isempty(line) && ~feof(fid), line = strtrim(fgetl(fid)); end
                if isempty(line), break; end 
                val = sscanf(line, '%f')';
                elmNum = val(1); elmType = val(2); numTags = val(3);
                physTag = 0; elemTag = 0; partTag=0;
                if numTags >= 1, physTag = val(4); end
                if numTags >= 2, elemTag = val(5); end
                if numTags >= 4, partTag=val(7); end
                nodes = val(3 + numTags + 1 : end);
                tempBlocks{elmType}{end+1} = [elmNum, physTag, elemTag, partTag, nodes];
            end
        end
    else
        if loadToMemory
            allData = fread(fid, inf, '*int32'); totalLen = length(allData);
            nodeCountMap = zeros(1, 150); for i = 1:150, try nodeCountMap(i) = getNodesCount(i); catch, end, end
            tempIdx = ones(19, 1); for i = 1:19, tempBlocks{i} = cell(numElementsTotal, 1); end
            currCount = 0; idx = 1; 
            while currCount < numElementsTotal && idx <= totalLen
                elmType = allData(idx); numInBlock = allData(idx+1); numTags = allData(idx+2); idx = idx + 3; 
                numNodes = nodeCountMap(elmType);
                blockSize = 1 + numTags + numNodes; blockDataLen = blockSize * numInBlock;
                if idx + blockDataLen - 1 > totalLen, break; end
                rawData1D = allData(idx : idx + blockDataLen - 1); idx = idx + blockDataLen; 
                blockMatrix = double(reshape(rawData1D, blockSize, numInBlock))';
                elmNum_col = blockMatrix(:, 1);
                physTag_col = elmNum_col * 0; elemTag_col = elmNum_col * 0; partTag_col = elmNum_col * 0;
                if numTags >= 1, physTag_col = blockMatrix(:, 2); end
                if numTags >= 2, elemTag_col = blockMatrix(:, 3); end
                if numTags >= 4, partTag_col=blockMatrix(:,5); end
                nodes_cols = blockMatrix(:, 1 + numTags + 1 : end);
                tempBlocks{elmType}{tempIdx(elmType)} = [elmNum_col, physTag_col, elemTag_col, partTag_col, nodes_cols];
                tempIdx(elmType) = tempIdx(elmType) + 1;
                currCount = currCount + numInBlock;
            end
            for i = 1:19, if tempIdx(i) > 1, Mesh.Elements{i} = vertcat(tempBlocks{i}{1:tempIdx(i)-1}); end, end
        else
            % Disk I/O fallback
            currCount = 0;
            while currCount < numElementsTotal
                header = fread(fid, 3, 'int32');
                if isempty(header) || length(header) < 3, break; end 
                elmType = header(1); numInBlock = header(2); numTags = header(3);
                blockSize = 1 + numTags + getNodesCount(elmType); 
                blockData = double(fread(fid, [blockSize, numInBlock], '*int32')');
                elmNum_col = blockData(:, 1); physTag_col = zeros(numInBlock, 1); elemTag_col = zeros(numInBlock, 1);
                if numTags >= 1, physTag_col = blockData(:, 2); end
                if numTags >= 2, elemTag_col = blockData(:, 3); end
                nodes_cols = blockData(:, 1 + numTags + 1 : end);
                tempBlocks{elmType}{end+1} = [elmNum_col, physTag_col, elemTag_col, zeros(numInBlock,1), nodes_cols];
                currCount = currCount + numInBlock;
            end
            for i = 1:19, if ~isempty(tempBlocks{i}), Mesh.Elements{i} = vertcat(tempBlocks{i}{:}); end, end
        end
    end
    if isempty(Mesh.Elements{1}), for i=1:19, if ~isempty(tempBlocks{i}), Mesh.Elements{i}=vertcat(tempBlocks{i}{:}); end, end; end
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