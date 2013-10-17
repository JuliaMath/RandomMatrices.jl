#Calculate numerical derivative using:
#forward first-order finite difference for first element
#backward first-order finite difference for last element
#central first-order finite difference for all intermediate elements
function gradient(y::Vector{Float64}, x::Vector{Float64})
  dy = zeros(length(y))
  n = length(dy)
  dy[1] = (y[2]-y[1])/(x[2]-x[1])
  dy[n] = (y[n]-y[n-1])/(x[n]-x[n-1])
  for i=2:n-1
    dy[i] = (y[i+1]-y[i-1])/(x[i+1]-x[i-1])
  end
  return dy
end

