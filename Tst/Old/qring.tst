ring r=32003,(x,y,z),(C,dp);;
ideal q=maxideal(2);
qring qr=q;
ideal i=maxideal(1);
std(i);
i=x3,y2,z;
std(i);
module m=[x],[y],[z],[0,x],[0,y],[0,z],[0,0,x],[0,0,z],[0,0,y];
std(m);
m=[x],[0,y],[0,0,z];
std(m);
setring r;
kill qr;
q=x2-1;
qring qr=q;
ideal i=maxideal(3);
std(i);
module m=[x3,x,1],[y2+x2,x2y];
std(m);
i=0;
std(i);
setring r;
q=0;
qring qq=q;
std(maxideal(1));
$
