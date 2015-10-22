#jacobigrowth2.jl

import Distributions
function jacobigrowth2(m1,m2,p,beta)
   CJ=eye(p)
   SJ=zeros(p, p)
   for k=1:m2
       Ju = sqrt(rand(Chisq(beta),1,p))
       Ju *= sqrt( rand(Chisq(beta*p)) / rand(Chisq(beta*(k+m1-p))) ) / norm(Ju)
       _, _, _, CJ, SJ=svd(CJ,[SJ; Ju])
   end
   diag(CJ)
end
