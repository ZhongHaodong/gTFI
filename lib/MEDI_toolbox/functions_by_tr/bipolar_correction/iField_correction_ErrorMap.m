%函数功能：校正奇偶回波线形相位差
%方法：ErrorMap是由iField生成的奇偶回波相位差图，利用ErrorMap校正iField
function iField_correction = iField_correction_ErrorMap(iField,ErrorMap)
    matrix_size = size(iField);
    numvox = prod(matrix_size(1:end-1));
    numte = matrix_size(end);
    EchoNum = (1:numte)';
    EchoNum = repmat(EchoNum,[1 numvox]);
    ErrorMap = reshape(ErrorMap,[1 numvox]);
    E = exp(1i *(-1).^EchoNum .*repmat(-ErrorMap,[numte,1]));  
    s0 = permute(iField, [length(matrix_size) 1:length(matrix_size)-1]);
    s0 = reshape(s0,[numte numvox]);
    s0 = s0.*E;
    s0 = reshape(s0,[matrix_size(end),matrix_size(1:3)] );
    iField_correction = permute(s0,[2 3 4 1]);
end