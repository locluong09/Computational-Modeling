function out = boundaries(X)
    [row, col] = size(X);
    X_new  = zeros(row+2, col +2);
    for i = 1:row
        for j = 1:col
            X_new(i+1,j+1) = X(i,j);
        end
    end
    out =  X_new;
end