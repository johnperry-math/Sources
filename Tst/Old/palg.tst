LIB "elim.lib";
ring T=(5,a),(x,y,z,y(1),y(2),y(3)),(dp(3),dp(3));
minpoly=a2+2;
matrix M[24][95]=
(1),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^2+(-1)*y^2+y(1),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*z^2+y(2),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^6+y^6+y(3),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
x+(-2*a)*y,
(2*a)*x,
(2)*y^2,
0,
0,
(-2)*y*z,
(-1)*y^2,
0,
0,
0,
0,
0,
0,
z,
(2*a)*x*z+(-1)*y*z,
(-1*a)*x*z+y*z,
0,
0,
0,
(-1)*y^2*z,
0,
0,
0,
(-1)*x^2+(-1)*y^2+y(1),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*z^2+y(2),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^6+y^6+y(3),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(2)*x,
0,
0,
z,
(-1*a)*y*z,
(1*a)*y^2,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^2+(-1)*y^2+y(1),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*z^2+y(2),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^6+y^6+y(3),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*y,
0,
0,
0,
0,
0,
0,
0,
0,
x,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^2+(-1)*y^2+y(1),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*z^2+y(2),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^6+y^6+y(3),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1),
(-2*a)*y,
0,
(-2*a)*y^2,
0,
(-2)*y^3,
0,
0,
0,
0,
(-1*a)*y,
(-1)*x+(-1*a)*y,
(2*a)*x*y+(-2)*y^2,
0,
0,
(2)*x*y+(-2*a)*y^2,
(2)*x*y^2+(1*a)*y^3,
0,
0,
0,
0,
0,
(-1)*x^2+(-1)*y^2+y(1),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*z^2+y(2),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^6+y^6+y(3),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1),
(-1*a),
x+(-1*a)*y,
(1*a)*y^3,
0,
0,
x+(1*a)*y,
(-2)*x*z+(-1*a)*y*z,
(-2)*x*y+(2*a)*y^2,
0,
(-2)*y^3,
0,
0,
0,
0,
0,
0,
(-1)*x^2+(1*a)*x*y,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^2+(-1)*y^2+y(1),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*z^2+y(2),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^6+y^6+y(3),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(2)*y^2,
0,
0,
0,
0,
0,
(-1),
0,
0,
(-2*a)*y^2,
0,
0,
0,
(-1*a)*x*y^2+(2)*y^3,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^2+(-1)*y^2+y(1),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*z^2+y(2),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^6+y^6+y(3),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-2),
0,
(-1)*y^3,
0,
0,
0,
(2)*y*z,
(2)*y^2,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^2+(-1)*y^2+y(1),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*z^2+y(2),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^6+y^6+y(3),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(1*a),
0,
(-2)*x+(-2*a)*y,
0,
y^2,
0,
(2*a)*y^3,
0,
0,
0,
0,
0,
0,
(1*a)*y^3,
0,
0,
(-1)*x*y^3+(-2*a)*y^4,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^2+(-1)*y^2+y(1),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*z^2+y(2),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^6+y^6+y(3),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(2*a)*y^2,
0,
0,
(-1*a),
(2*a)*z,
0,
0,
y^2,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^2+(-1)*y^2+y(1),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*z^2+y(2),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^6+y^6+y(3),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(2),
0,
0,
0,
(-1*a)*y^2,
0,
0,
0,
0,
(-1),
0,
(2)*x+(1*a)*y,
0,
(2*a)*x*y^2+(2)*y^3,
(-1*a)*x+(2)*y,
(-1*a)*x*y+(-2)*y^2,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^2+(-1)*y^2+y(1),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*z^2+y(2),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^6+y^6+y(3),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-2),
y^2,
0,
0,
(1),
(-2)*z,
0,
0,
(-1*a)*y^2,
0,
y^3,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^2+(-1)*y^2+y(1),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*z^2+y(2),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^6+y^6+y(3),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(2*a),
0,
(2)*y,
0,
(-1*a)*y^2,
0,
0,
0,
0,
0,
0,
(2*a)*y^2,
(-1*a),
0,
(-2)*x*y^2+(1*a)*y^3,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^2+(-1)*y^2+y(1),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*z^2+y(2),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^6+y^6+y(3),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
x+(1*a)*y,
0,
0,
0,
0,
(-1*a),
0,
(2)*y,
0,
(2*a)*x^2+(-1)*x*y+(-1*a)*y^2,
0,
0,
0,
0,
(-2),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^2+(-1)*y^2+y(1),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*z^2+y(2),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^6+y^6+y(3),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-2),
0,
0,
0,
(-2)*y^2,
0,
0,
0,
0,
(2*a),
0,
0,
0,
(1*a)*x+(-2)*y,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^2+(-1)*y^2+y(1),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*z^2+y(2),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^6+y^6+y(3),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
y,
0,
0,
0,
0,
(-2),
0,
0,
0,
(-2)*y^2,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^2+(-1)*y^2+y(1),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*z^2+y(2),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^6+y^6+y(3),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1),
0,
(-2*a)*y,
0,
0,
0,
0,
0,
0,
(-1*a)*y,
0,
0,
x*y+(2*a)*y^2,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^2+(-1)*y^2+y(1),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*z^2+y(2),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^6+y^6+y(3),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-2*a),
0,
0,
0,
0,
0,
0,
(-1),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^2+(-1)*y^2+y(1),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*z^2+y(2),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^6+y^6+y(3),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(1*a),
0,
0,
0,
0,
0,
0,
0,
0,
(-2*a)*x+(-2)*y,
0,
(-2),
(-1)*y^2,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^2+(-1)*y^2+y(1),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*z^2+y(2),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^6+y^6+y(3),
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1),
0,
0,
0,
0,
0,
0,
(1*a),
0,
(-1)*y,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^2+(-1)*y^2+y(1),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*z^2+y(2),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^6+y^6+y(3),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(1*a),
0,
0,
0,
0,
0,
0,
(-2*a),
0,
0,
(2)*x+(-1*a)*y,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^2+(-1)*y^2+y(1),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*z^2+y(2),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^6+y^6+y(3),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(2),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^2+(-1)*y^2+y(1),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*z^2+y(2),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^6+y^6+y(3),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(2),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^2+(-1)*y^2+y(1),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*z^2+y(2),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^6+y^6+y(3),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(1),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^2+(-1)*y^2+y(1),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*z^2+y(2),
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
(-1)*x^6+y^6+y(3)
;
matrix N=std(M);
N=elim(N,1,3);
N;
module(N);

print(N);
$
