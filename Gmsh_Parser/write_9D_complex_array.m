function write_9D_complex_array(filename, array)
    fid = fopen(filename, 'w');
    if fid == -1
        error('파일을 생성할 수 없습니다: %s', filename);
    end
    
    [num_rows, num_cols] = size(array);
    
    % 배열을 순회하며 한 줄에 9개의 값을 띄어쓰기로 구분하여 출력
    for i = 1:num_rows
        for j = 1:num_cols
            val = array(i, j);
            
            if isreal(val)
                % 실수일 경우
                fprintf(fid, '%.15e ', val);
            else
                % 복소수일 경우 부호 처리
                if imag(val) >= 0
                    fprintf(fid, '%.15e+%.15ei ', real(val), imag(val));
                else
                    fprintf(fid, '%.15e%.15ei ', real(val), imag(val));
                end
            end
        end
        fprintf(fid, '\n'); % 9개 다 찍으면 줄바꿈
    end
    
    fclose(fid);
end