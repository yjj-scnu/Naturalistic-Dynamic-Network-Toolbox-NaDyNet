function mat=sf_vec2mat(V,vec)
vec=vec(:);
mat=zeros(V,V);
k=0;
for i=1:(V-1)
    for j=(i+1):V
        k=k+1;
        mat(j,i)=vec(k);
    end
end

temp=ones(V);
IND = find((temp-triu(temp))>0);

vec2 = mat(IND);
if ~isequal(vec,vec2)
    error('Error: vector size does not match, please check')
end
end