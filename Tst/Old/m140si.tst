  ring r=32003,(x,y,z),dp;
  ideal i=x+y,z3+22y;
  write(">save_i",i);
  ring r0=0,(x,y,z),Dp;
  string s="ideal k="+read("save_i")+";";
  s;
  execute s; // define the ideal k
  k;
$
