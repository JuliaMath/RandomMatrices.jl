import Distributions

#Compute Mauchly's statistic (valid under assumption of multinormality)
#Mauchly, 1940; Kendall and Stuart, 1968
#n is the number of data points (samples)
#This returns both the value of the test statistic and the expected distribution to test against
function MauchlySphericityTestStatistic(PopulationCovariance::AbstractMatrix{T}, SampleCovariance::AbstractMatrix{T}, n::Integer) where {T<:Real}
    p = size(SampleCovariance)[1]
    C = inv(PopulationCovariance) * SampleCovariance
    l = (det(C)/(trace(C)/p)^p)
    x = -n*log(l)
    f = p*(p+1)/2-1
    return x, Chisq(f)
end

#Johnstone's variant of the sphericity test
#n, p>= 10 recommended
#Johnstone (2001)
function JohnstoneSphericityTestStatistic(PopulationCovariance::AbstractMatrix{T}, SampleCovariance::AbstractMatrix{T}, n::Integer) where {T<:Real}
    C = inv(PopulationCovariance) * SampleCovariance
    v = max(eigvals(C))
    mu=(sqrt(n-1)+sqrt(p))^2
    sigma=sqrt(mu)*(1/sqrt(n-1)+1/sqrt(p))^(1/3)
    (n*l-mu)/sigma #To be tested against Tracy-Widom with beta=1
end
