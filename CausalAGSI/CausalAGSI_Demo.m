if false 
n=100;
input.Speed=ones(2*n-1,2*n-1);
input.Seeds=[n;n;0];
input.tol=1e-6;
output=CausalAGSI_Iso2(input);
imagesc(output.Distance);
end;

if true
n=10;
input.Speed=ones(n,n);
input.Seeds=[n/2;n/2;0];
input.tol=1e-6;
output=CausalAGSI_Iso2(input);
imagesc(output.Distance);
end;