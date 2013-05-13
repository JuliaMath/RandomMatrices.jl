#jacobian2by2.jl
# Code 8.1 of Random Eigenvalues by Alan Edelman

#Experiment:    Compute the Jacobian of a 2x2 matrix function
#Comment:       Symbolic tools are not perfect. The author
#               exercised care in choosing the variables

warn("Julia doesn't currently have a good symbolic math package")
warn("This may use features of PyCall and/or SymPy which are not publicly available yet")
using PyCall
pyinitialize("/opt/local/bin/python2")
using SymPy

p,q,r,s, a,b,c,d, t, e1,e2 = Sym(:p,:q,:r,:s, :a,:b,:c,:d, :t, :e1,:e2)
X = SymMatrix([p q; r s])
A = SymMatrix([a b; c d])
jacobian(Y,X) = Sym(Y[:reshape](1,length(Y))[:jacobian](X[:reshape](1,length(X))))
Y = X^2; J = jacobian(Y,X); JAC_square = Sym(J[:det]()[:factor]())
Y = X^3; J = jacobian(Y,X); JAC_cube = Sym(J[:det]()[:factor]())
Y = X[:inv](); J = jacobian(Y,X); JAC_cube = Sym(J[:det]()[:factor]())
Y = A*X; J = jacobian(Y,X); JAC_linear = Sym(J[:det]()[:factor]())
Y = SymMatrix([p q; r/p X[:det]()/p]); J = jacobian(Y,X); JAC_lu = Sym(J[:det]()[:factor]())

x = SymMatrix([p, s, r])
y = SymMatrix([sqrt(p), sqrt(s), r/(sqrt(p)*sqrt(s))])
J = jacobian(y,x)
JAC_DMD = Sym(J[:det]()[:factor]())

x = SymMatrix([p, s])
y = SymMatrix([sqrt(p^2 + s^2), atan(s/p)])
J = jacobian(y,x)
JAC_notrace = Sym(J[:det]()[:factor]())

Q = SymMatrix([cos(t)  -sin(t); sin(t)  cos(t)])
D = SymMatrix([e1 0; 0 e2])
Y = Q*D*Q[:transpose]()
y = SymMatrix([Y[1,1], Y[2,2], Y[1,2]])
x = SymMatrix([t, e1, e2])
J = jacobian(Y,X)
JAC_symeig = Sym(J[:det]()[:factor]())

X = SymMatrix([p s; s r])
Y = Sym(A[:transpose]())*X*A
y = SymMatrix([Y[1,1], Y[2,2], Y[1,2]])
x = SymMatrix([p, r, s])
J = jacobian(Y,X)
JAC_symcong = Sym(J[:det]()[:factor]())

