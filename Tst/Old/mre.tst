//
// test script for mres command
//
pagelength = 10000;
option(prot);
ring r1=32003,(x,y,z),dp;
r1;
"-----------------------";
ideal i=x2,y2,z2;
i;
mres(i,0); 
"--------------------------";
module m=[x2,0],[0,y4];
m;
mres(m,0);
"-------------------------------";
ring r2=31991,(t,x,y,z,w),ls;
ideal j=t2x2+tx2y+x2yz,t2y2+ty2z+y2zw,t2z2+tz2w+xz2w,t2w2+txw2+xyw2;
j;
mres(j,0);
"-------------------------";
listvar(all);
kill r1;
$;

