# Generates joint probability density functions (jpdfs) of eigenvalues
# of finite-dimensional random matrices that are distributed according
# to the classical Gaussian random matrix ensembles
#
# The notation follows closely that in:
#
# Ioana Dumitriu and Alan Edelman, "Matrix Models for Beta Ensembles"
# Journal of Mathematical Physics, vol. 43 no. 11 (2002), pp. 5830--5547
# doi: 10.1063/1.1507823
# arXiv: math-ph/0206043

export VandermondeDeterminant, HermiteJPDF, LaguerreJPDF, JacobiJPDF 

###
# Joint probability density functions
#

#Calculate Vandermonde determinant term
function VandermondeDeterminant{Eigenvalue<:Number}(lambda::Vector{Eigenvalue}, beta::Real)
    n = length(lambda)
    Vandermonde = 1.0
    for j=1:n
        for i=1:j-1
            Vandermonde *= abs(lambda[i] - lambda[j])^beta
        end
    end
    return Vandermonde
end


function HermiteJPDF{Eigenvalue<:Number}(lambda::Vector{Eigenvalue}, beta::Real)
    n = length(lambda)
    #Calculate normalization constant
    c = (2pi)^(-n/2)
    for j=1:n
        c *= gamma(1 + beta/2)/gamma(1 + beta*j/2)
    end

    #Calculate argument of exponential
    Energy = sum(lambda.^2/2)

    return c * VandermondeDeterminant(lambda, beta) * exp(-Energy)
end



#TODO Check m and ns
function LaguerreJPDF{Eigenvalue<:Number}(lambda::Vector{Eigenvalue}, n::Unsigned, beta::Real)
    m = length(lambda)
    #Laguerre parameters
    a = beta*n/2.0
    p = 1.0 + beta*(m-1)/2.0
    #Calculate normalization constant
    c = 2.0^-(m*a)
    for j=1:m
       z = (a - beta*(m-j)/2.0)
       if z < 0 && (int(z) - z) < eps()
           #Pole of gamma function, there is no density here no matter what
           return 0.0
       end
       c *= gamma(1 + beta/2)/(gamma(1 + beta*j/2)*gamma(z))
    end

    #Calculate Laguerre product term
    Prod = prod(lambda.^(a-p))

    #Calculate argument of exponential
    Energy = sum(lambda)/2

    return c * VandermondeDeterminant(lambda, beta) * Prod * exp(-Energy)
end


#TODO Check m and ns
function JacobiJPDF{Eigenvalue<:Number}(lambda::Vector{Eigenvalue}, n1::Unsigned, n2::Unsigned, beta::Real)
    m = length(lambda)
    #Jacobi parameters
    a1 = beta*n1/2.0
    a2 = beta*n2/2.0
    p = 1.0 + beta*(m-1)/2.0
    #Calculate normalization constant
    c = 1.0
    for j=1:m
       z1 = (a1 - beta*(m-j)/2.0)
       if z1 < 0 && (int(z1) - z1) < eps()
           #Pole of gamma function, there is no density here no matter what
           return 0.0
       end
       z2 = (a2 - beta*(m-j)/2.0)
       if z2 < 0 && (int(z2) - z2) < eps()
           #Pole of gamma function, there is no density here no matter what
           return 0.0
       end
       c *= gamma(1 + beta/2)*gamma(a1+a2-beta*(m-j)/2)
       c /= gamma(1 + beta*j/2)*gamma(z1)*gamma(z2)
    end

    #Calculate Laguerre product term
    Prod = prod(lambda.^(a1-p))*prod((1-lambda).^(a2-p))

    #Calculate argument of exponential
    Energy = sum(lambda/2)

    return c * VandermondeDeterminant(lambda, beta) * Prod * exp(-Energy)
end

