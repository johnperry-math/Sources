//
// test script for syz command
//
pagelength = 10000;
ring r1 = 32003,(x,y,z),(c,ls);
r1;
"-------------------------------";
ideal i1=maxideal(3);
ideal i2=x6,y4,z5,xy,yz,xz;
i2;
syz(i2);
"---------------------------------";
i1;
syz(i1);
"----------------------------------";
vector v1=[x4,y2,z2,2z2];
vector v2=[yz,xy,xz,2xy];
module m=v1,v2;
m;
syz(m);
"--------------------------------";
listvar(all);
kill r1;
$;
