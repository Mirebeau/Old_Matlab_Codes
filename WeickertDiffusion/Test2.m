if false
m=[1,2;3,4];
[E,V]=eig(m)

m=[1,2,5;1,2,6];


arrayfun(@(x,y,z)eig([x,y;y,z]),m(:,1),m(:,2),m(:,3),'uniformOutput',false)



[x,y]=meshgrid(1:10,1:10);
x=x'; y=y';
options.a=0;
%st=StructureTensor(x,options);

a=rand(5); a=a'*a;
[v,d]=eig(a);
disp(d); disp(v);
disp(v*d*v');
disp(a);

end


%img=hdf5read('Images/Source/Cos3D_Noisy.hdf5','/ITKImage/0/VoxelData');

options.Weickert_lambda=1;
tensors=zeros(1,1,1,6);
tensors(:,:,:,1)=1;
tensors(:,:,:,3)=2;
tensors(:,:,:,6)=3;
WeickertTensor_3D(tensors,options);

%img2=imfilter(img,[1,1,1],'symmetric');

a=[4,5,6];
a([1,3,1]) = a([1,3,1])+[2,2,3];
TestFun(a);
b=[a;a];


mex -Dchar16_t=uint16_T -output TripletMatrixProduct Toolbox/mex/TripletMatrixProduct.cpp

triplets = [1,2,1;2,2,1;1,2,3];
vector = [1,2;3,4];
times = [1,2];
[result,done]=TripletMatrixProduct(triplets,vector,times);
A=sparse(triplets(1,:),triplets(2,:),triplets(3,:),2,2);
result2=A*vector;