function axset(dim1, dim2, dim3)

if nargin == 1
    dim2 = dim1;
    dim3 = dim1;
end

quiver3(0,0,0,dim1,0,0)
quiver3(0,0,0,0,dim2,0)
quiver3(0,0,0,0,0,dim3)

text(dim1,0,0,'x')
text(0,dim2,0,'y')
text(0,0,dim3,'z')

end