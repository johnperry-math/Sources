LIB "tst.lib"; tst_init();
LIB "gmssing.lib";
ring R=0,(x,y),ds;
poly t=x5+x2y2+y5;
def G=gmsring(t,"s");
setring(G);
list l0=gmsnf(gmspoly,0);
print(l0[1]);
list l1=gmsnf(gmspoly,1);
print(l1[1]);
list l=gmsnf(l0[2],1);
print(l[1]);
tst_status(1);$
