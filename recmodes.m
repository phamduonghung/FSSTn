function imf = recmodes(T,C,d)


[nr N] = size(C);
na = size(T,1);
imf = zeros(nr,N);

for k=1:nr
    for j=1:N
        imf(k,j) = sum(T(max(1,C(k,j)-d):min(na,C(k,j)+d),j));
    end
end

%imf = 2*real(imf);
