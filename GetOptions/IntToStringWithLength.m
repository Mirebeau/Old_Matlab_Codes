function s = IntToStringWithLength(i,n)
    for k=1:n; 
        if i<10^k; 
            s=sprintf('%s%i',repmat('0',[1,n-k]),i);
            return;
        end;
    end;
    disp(sprintf('IntToStrinWithLength error: %i exceeds 10^%i',i,n));
    s=sprintf('%i',i);
end