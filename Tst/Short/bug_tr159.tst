LIB "tst.lib";
tst_init();

// NF was not 0 as expected - not enough searching in kFind...InS
intmat m1[8][8] = 
 2, 2, 2, 4, 1, 1, 3, 3,
 0, 0, 0,-1, 0, 0, 0, 0,
 0, 0, 0, 0,-1, 0, 0, 0,
-1, 0, 0, 0, 0, 0, 0, 0,
 0,-1, 0, 0, 0, 0, 0, 0,
 0, 0,-1, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 0,-1, 0, 0,
 0, 0, 0, 0, 0, 0,-1, 0;
ring R1 = 2, (b_2_1,b_2_2,b_2_3,c_4_8,a_1_0,b_1_1,b_3_4,b_3_5), M(m1);
ideal I;
I[1]=a_1_0^2;
I[2]=a_1_0*b_1_1;
I[3]=b_2_3*a_1_0+b_2_2*a_1_0;
I[4]=b_2_1*b_1_1+b_2_2*a_1_0;
I[5]=b_2_2*b_1_1+b_2_2*a_1_0;
I[6]=a_1_0*b_3_4;
I[7]=a_1_0*b_3_5;
I[8]=b_2_2^2+b_2_1*b_2_3;
I[9]=b_1_1*b_3_4;
I[10]=b_2_2*b_3_4+b_2_1*b_3_5+b_2_1*b_2_2*a_1_0;
I[11]=b_2_3*b_3_4+b_2_2*b_3_5+b_2_1*b_2_2*a_1_0;
I[12]=b_3_4^2+b_2_1*b_2_3^2+b_2_1^2*b_2_3;
I[13]=b_3_4*b_3_5+b_2_2*b_2_3^2+b_2_1*b_2_2*b_2_3;
I[14]=b_3_5^2+b_2_3*b_1_1*b_3_5+b_2_3^3+b_2_1*b_2_3^2+c_4_8*b_1_1^2;
qring Q1 = groebner(I);


intmat m2[5][5] =
 2, 2, 2, 1, 1,
 0,-1,-1, 0, 0,
 0, 0, 0,-1, 0,
-1, 0, 0, 0, 0,
 0,-1, 0, 0, 0;
ring R2 = 2, (@b_2_1,@c_2_2,@c_2_3,@a_1_0,@b_1_1), M(m2);
ideal I;
I[1]=@a_1_0^2;
I[2]=@a_1_0*@b_1_1;
I[3]=@b_2_1*@a_1_0;
I[4]=@b_2_1^2+@c_2_2*@b_1_1^2;
qring Q2 = groebner(I);

def Q = Q1+Q2;
setring Q;
Q;

ideal t = a_1_0, b_1_1, b_2_1, b_2_2, b_2_3, b_3_4, b_3_5, c_4_8;
option(prot);
NF(b_2_1,t); 

tst_status(1);$
