ring r=32003,(x,y,z),dp;
poly f=x+y+z;
poly g=1+f+f^2+f^3+f^4;
f;
g;
[f,1];
[f,x^5];
ring r1=32003,(x,y,z),(c,dp);
poly f=fetch(r,f);
poly g=fetch(r,g);
f;
g;
[f,1];
[f,x^5];
// module orders
ring r2=0,x,(ds,c);
matrix m[2][3]=1,2,3,4,5;
print(module(transpose(m)));
"-------";
ring r3=0,x,(dp,c);
matrix m=fetch(r2,m);
print(module(transpose(m)));
"-------";
ring r4=0,x,(Dp,c);
matrix m=fetch(r2,m);
print(module(transpose(m)));
"-------";
ring r5=0,x,(Ds,c);
matrix m=fetch(r2,m);
print(module(transpose(m)));
"-------";
ring r6=0,x,(ls,c);
matrix m=fetch(r2,m);
print(module(transpose(m)));
"-------";
ring r7=0,x,(lp,c);
matrix m=fetch(r2,m);
print(module(transpose(m)));
"-------";
ring r8=0,x,(ds,C);
matrix m=fetch(r2,m);
print(module(transpose(m)));
"-------";
ring r9=0,x,(dp,C);
matrix m=fetch(r2,m);
print(module(transpose(m)));
"-------";
ring r10=0,x,(Dp,C);
matrix m=fetch(r2,m);
print(module(transpose(m)));
"-------";
ring r11=0,x,(Ds,C);
matrix m=fetch(r2,m);
print(module(transpose(m)));
"-------";
ring r12=0,x,(ls,C);
matrix m=fetch(r2,m);
print(module(transpose(m)));
"-------";
ring r13=0,x,(lp,C);
matrix m=fetch(r2,m);
print(module(transpose(m)));

"-------";
ring R2=0,x,(c,ds);
matrix m=fetch(r2,m);
print(module(transpose(m)));
"-------";
ring R3=0,x,(c,dp);
matrix m=fetch(r2,m);
print(module(transpose(m)));
"-------";
ring R4=0,x,(c,Dp);
matrix m=fetch(r2,m);
print(module(transpose(m)));
"-------";
ring R5=0,x,(c,Ds);
matrix m=fetch(r2,m);
print(module(transpose(m)));
"-------";
ring R6=0,x,(c,ls);
matrix m=fetch(r2,m);
print(module(transpose(m)));
"-------";
ring R7=0,x,(c,lp);
matrix m=fetch(r2,m);
print(module(transpose(m)));
"-------";
ring R8=0,x,(C,ds);
matrix m=fetch(r2,m);
print(module(transpose(m)));
"-------";
ring R9=0,x,(C,dp);
matrix m=fetch(r2,m);
print(module(transpose(m)));
"-------";
ring R10=0,x,(C,Dp);
matrix m=fetch(r2,m);
print(module(transpose(m)));
"-------";
ring R11=0,x,(C,Ds);
matrix m=fetch(r2,m);
print(module(transpose(m)));
"-------";
ring R12=0,x,(C,ls);
matrix m=fetch(r2,m);
print(module(transpose(m)));
"-------";
ring R13=0,x,(C,lp);
matrix m=fetch(r2,m);
print(module(transpose(m)));

$
