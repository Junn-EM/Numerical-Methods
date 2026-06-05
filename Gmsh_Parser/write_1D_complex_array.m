function write_1D_complex_array(filename, array)
    fid = fopen(filename, 'w');
    if fid == -1
        error('파일을 생성할 수 없습니다: %s', filename);
    end
    
    % 배열을 순회하며 한 줄(1열)씩 출력
    for i = 1:length(array)
        val = array(i);
        
        if isreal(val)
            % 실수일 경우 깔끔하게 소수점 6자리까지 출력
            fprintf(fid, '%.6f\n', val);
        else
            % 복소수일 경우 부호에 맞춰 'a+bi' 또는 'a-bi' 형태로 출력
            if imag(val) >= 0
                fprintf(fid, '%.6f+%.6fi\n', real(val), imag(val));
            else
                fprintf(fid, '%.6f%.6fi\n', real(val), imag(val)); % imag가 음수면 자체 마이너스 부호 출력됨
            end
        end
    end
    
    fclose(fid);
end