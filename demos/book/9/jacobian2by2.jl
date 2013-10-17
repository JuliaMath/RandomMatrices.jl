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
jacobian(Y,X) = Sym(Y[:reshape](1,length(Y))[:jacobian](X[:reshape](1,length(X))))

p,q,r,s, a,b,c,d, t, e1,e2 = Sym(:p,:q,:r,:s, :a,:b,:c,:d, :t, :e1,:e2)
X = SymMatrix([p q; r s])
A = SymMatrix([a b; c d])

Y = X^2; J = jacobian(Y,X); JAC_square = Sym(J[:det]()[:factor]()); println("square:\n", J, "\n\n", JAC_square, "\n\n")
Y = X^3; J = jacobian(Y,X); JAC_cube = Sym(J[:det]()[:factor]()); println("cube:\n", J, "\n\n", JAC_cube, "\n\n")
Y = X[:inv](); J = jacobian(Y,X); JAC_inv = Sym(J[:det]()[:factor]()); println("inv:\n", J, "\n\n", JAC_inv, "\n\n")
Y = A*X; J = jacobian(Y,X); JAC_linear = Sym(J[:det]()[:factor]()); println("linear:\n", J, "\n\n", JAC_linear, "\n\n")
Y = SymMatrix([p q; r/p X[:det]()/p]); J = jacobian(Y,X); JAC_lu = Sym(J[:det]()[:factor]()); println("lu:\n", J, "\n\n", JAC_lu, "\n\n")

x = SymMatrix([p, s, r])
y = SymMatrix([sqrt(p), sqrt(s), r/(sqrt(p)*sqrt(s))])
J = jacobian(y,x)
JAC_DMD = Sym(J[:det]()[:factor]())
println("DMD:\n", J, "\n\n", JAC_DMD, "\n\n")

x = SymMatrix([p, s])
y = SymMatrix([sqrt(p^2 + s^2), atan(s/p)])
J = jacobian(y,x)
JAC_notrace = Sym(J[:det]()[:factor]())
println("notrace:\n", J, "\n\n", JAC_notrace, "\n\n")

Q = SymMatrix([cos(t)  -sin(t); sin(t)  cos(t)])
D = SymMatrix([e1 0; 0 e2])
Y = Q*D*Q[:transpose]()
y = SymMatrix([Y[1,1], Y[2,2], Y[1,2]])
x = SymMatrix([t, e1, e2])
J = jacobian(Y,X)
JAC_symeig = Sym(J[:det]()[:factor]())
println("symeig:\n", J, "\n\n", JAC_symeig, "\n\n")

X = SymMatrix([p s; s r])
Y = Sym(A[:transpose]())*X*A
y = SymMatrix([Y[1,1], Y[2,2], Y[1,2]])
x = SymMatrix([p, r, s])
J = jacobian(Y,X)
JAC_symcong = Sym(J[:det]()[:factor]())
println("symcong:\n", J, "\n\n", JAC_symcong, "\n\n")

