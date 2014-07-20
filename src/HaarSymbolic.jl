using Catalan

export permutations_in_Sn, compose, cycle_structure, data, part, #Functions working with partitions and permutations
    UniformHaar, expectation, expectedtrace, WeingartenUnitary

#Functions working with partitions and permutations
# TODO Migrate these functions to Catalan

#Iterate over partitions of n in lexicographic order
function part(n::Integer)
    if n==0 produce([]) end
    if n<0 return end
    for p in @task part(n-1)
        p = [p; 1]
        produce(p)
        p = p[1:end-1]
        if length(p) == 1 || (length(p)>1 && p[end]<p[end-1])
            p[end] += 1
            produce(p)
            p[end] -= 1
        end
    end
end

if _HAVE_GSL
function permutations_in_Sn(n::Integer)
    P = permutation_calloc(n)
    while true 
        produce(P)
        try permutation_next(P) catch break end
    end
end

function compose(P::Ptr{gsl_permutation}, Q::Ptr{gsl_permutation})
    #Compose the permutations
    n=convert(Int64, permutation_size(P))
    @assert n==permutation_size(Q)
    x=permutation_alloc(n)
    Pv = [x+1 for x in pointer_to_array(permutation_data(P), (n,))]
    Qv = [x+1 for x in pointer_to_array(permutation_data(Q), (n,))]
    Xv = [Qv[i] for i in Pv]
    for i=1:n
         unsafe_assign(permutation_data(x), Xv[i]-1, i)
    end
    permutation_valid(x)
    x
end

#Compute cycle structure, i.e. lengths of each cycle in the cycle factorization of the permutatation
function cycle_structure(P::Ptr{gsl_permutation})
    n=convert(Int64, permutation_size(P))
    Pcanon = permutation_linear_to_canonical(P)
    PCANON = data(Pcanon)
    HaveNewCycleAtPos = [PCANON[i]>PCANON[i+1] for i=1:n-1]
    cycleindices = Int[1]
    for i=1:n-1 if HaveNewCycleAtPos[i] push!(cycleindices, i+1) end end
    push!(cycleindices, n+1)
    cyclestructure = Int[cycleindices[i+1]-cycleindices[i] for i=1:length(cycleindices)-1]
    @assert sum(cyclestructure) == n
    cyclestructure
end

#Returns a vector of indices (starting from 1 in the Julia convention)
data(P::Ptr{gsl_permutation}) = [convert(Int64, x)+1 for x in
    pointer_to_array(permutation_data(P), (convert(Int64, permutation_size(P)) ,))]

end #_HAVE_GSL

immutable UniformHaar <: ContinuousMatrixDistribution
    beta::Float64
    N::Int
    UniformHaar(beta::Real, N::Integer) = new(float64(beta), int(N))
end
    
# In random matrix theory one often encounters expressions of the form
#
#X = Q * A * Q' * B
#
#where A and B are deterministic matrices with fixed numerical matrix entries
#and Q is a random matrix that does not have explicitly defined matrix
#elements. Instead, one takes an expectation over of expressions of this form
#and this "integrates out" the Qs to produce a numeric result.
#
#expectation(Q, X) #= some number
#
#Here is a function that symbolically calculates the expectation of a product
#of matrices over the symmetric group that Q is uniform Haar over.
#It takes an expression consisting of a product of matrices and replaces it
#with an evaluated symbolic expression which is the expectation.
expectation(N::Integer, MyQ::Symbol, X::Symbol) = Q==X ? 0 : X
expectation(N::Integer, MyQ::Symbol, X::Expr) = expectation(N, MyQ, X, false)
function expectation(N::Integer, MyQ::Symbol, X::Expr, dotrace::Bool)
    if X.head != :call
        error(string("Unexpected type of expression: ", X.head))
    end
    
    n = length(X.args) - 1
    if n < 1 return X end #nothing to do Haar-wise
    
    if X.args[1] != :*
        error("Unexpected expression, only products supported")
    end

    # Parse expression involving products of matrices to extract the
    # positions of Haar matrices and their ctransposes
    Qidx=[] #Indices for Haar matrices
    Qpidx=[] #Indices for ctranspose of Haar matrices
    Others=[]
    for i=1:n
        thingy=X.args[i+1]
        if isa(thingy, Symbol)
            if thingy == MyQ
                Qidx=[Qidx; i]
            else
                Others = [Others; (thingy, i, i+1)]
            end
            println(i, ' ', thingy)
        elseif isa(thingy, Expr)
            if thingy.head==symbol('\'') && length(thingy.args)>=1 #Maybe this is a Q'
                if typeof(thingy.args[1])==Symbol && thingy.args[1] == MyQ
                    println(i, ' ', thingy, "::", MyQ,"'")
                    Qpidx=[Qpidx; i]
                end
            end
        else
            error(string("Unexpected subexpression ", thingy ," of type ", typeof(thingy)))
        end
    end
    if length(Qidx) == length(Qpidx) == 0 return X end #nothing to do Haar-wise
    println(MyQ, " is in places ", Qidx)
    println(MyQ, "' is in places ", Qpidx)

    M = length(Qidx)
    #If there are different Qs and Q's, the answer is a big fat 0
    #TODO return correct size of 0
    if M != length(Qpidx) return 0 end

    ##################################
    # Step 2. Enumerate permutations #
    ##################################
    AllTerms = {}
    for sigma in @task permutations_in_Sn(M)
        for tau in @task permutations_in_Sn(M)
            sigma_inv=permutation_inverse(sigma)
            #Compose the permutations
            perm=compose(sigma_inv, tau)
            
            #Compute cycle structure, i.e. lengths of each
            #cycle in the cycle factorization of the
            #permutatation
            cyclestruct = cycle_structure(perm)
            Qr = Int64[n+1 for n in Qidx]
            Qpr= Int64[Qpidx[n]+1 for n in data(tau)]
            #Change indices to wrap around if this is a
            #trace expression
            if dotrace
                Qr[Qr.==n+1]=1
                Qpr[Qpr.==n+1]=1
            end
            #Consolidate deltas
            Deltas=Dict()
            for i=1:M
                Deltas[Qr[i]] = Qpidx[data(sigma)[i]]
                Deltas[Qidx[i]] =  Qpr[data(tau)[i]]
            end
            #Print deltas
            print("V(", cyclestruct, ") ")
            for k in keys(Deltas)
                print("δ(",k,",",Deltas[k],") ")
            end
            for (Symb, col_idx, row_idx) in Others
                print(Symb,"(",col_idx,",",row_idx,") ")
            end
            println()
            print("= V(", cyclestruct, ") ")
            #Evaluate deltas over the indices of the other
            #matrices
            ReindexedSymbols = []
            for (Symb, col_idx, row_idx) in Others
                new_col, new_row = col_idx, row_idx
                if has(Deltas, col_idx) new_col = Deltas[col_idx] end
                if has(Deltas, row_idx) new_row = Deltas[row_idx] end
                print(Symb,"(",new_col,",",new_row,") ")
                ReindexedSymbols = [ReindexedSymbols; (Symb, new_col, new_row)]
            end
            println()
            #Parse coefficient
            Coefficient= WeingartenUnitary(N, perm)
            println("V(",cyclestruct,") = ",Coefficient)
            #Reconstruct expression
            println("START PARSING")
            println("The term is =", ReindexedSymbols[end])
            Symb, left_idx, right_idx = pop!(ReindexedSymbols)
            Expression={{{Symb}, left_idx, right_idx}}
            while length(ReindexedSymbols) > 0
                pop_idx = expr_idx = do_transpose = is_left = nothing
                for expr_iter in enumerate(Expression)
                    expr_idx, expr_string = expr_iter
                    _, left_idx, right_idx = expr_string
                    for iter in enumerate(ReindexedSymbols)
                        idx, data = iter
                        Symb, col_idx, row_idx = data
                        if row_idx == left_idx
                            pop_idx = idx
                            do_transpose = false
                            is_left = true
                            println("Case A")
                            break
                        elseif col_idx == right_idx
                            pop_idx = idx
                            do_transpose = false
                            is_left = false
                            println("Case B")
                            break
                        elseif row_idx == right_idx
                            pop_idx = idx
                            do_transpose = true
                            is_left = false
                            println("Case C")
                            break
                        elseif col_idx == left_idx
                            pop_idx = idx
                            do_transpose = true
                            is_left = true
                            println("Case D")
                            break
                        end
                    end
                    if pop_idx != nothing break end
                end
                println("Terms left = ", ReindexedSymbols)
                if pop_idx == nothing #Found nothing, start new expression blob
                println("New word: The term is ", ReindexedSymbols[1])
                    Symb, left_idx, right_idx = delete!(ReindexedSymbols, 1)
                    push!(Expression, {{Symb}, left_idx, right_idx})
                else #Found something
                    println("The term is =", ReindexedSymbols[pop_idx])
                    Symb, col_idx, row_idx = delete!(ReindexedSymbols, pop_idx)
                    Term = do_transpose ? Expr(symbol("'"), Symb) : Symb
                    if is_left
                        insert!(Expression[expr_idx][1], 1, Term)
                        Expression[expr_idx][2] = do_transpose ? row_idx : col_idx
                    else
                        push!(Expression[expr_idx][1], Term)
                        Expression[expr_idx][3] = do_transpose ? col_idx : row_idx
                    end
                end
            end
            println("DONE PARSING TREE")

            # Evaluate closed cycles
            NewExpression={}
            for ExprChunk in Expression
                ExprBlob, start_idx, end_idx = ExprChunk
                print(ExprChunk, " => ")
                if start_idx == end_idx #Have a cycle; this is a trace
                    ex= length(ExprBlob)==1 ? :(trace($(ExprBlob[1]))) : 
                        Expr(:call, :trace, Expr(:call, :*, ExprBlob...))
                else #Not a cycle, regular chain of multiplications
                    ex= length(ExprBlob)==1 ? ExprBlob[1] :
                        Expr(:call, :*, ExprBlob...)
                end
                push!(NewExpression, ex)
                println(ex)
            end
            ex = Expr(:call, :*, Coefficient, NewExpression...)
            println("Final expression is: ", ex)
            push!(AllTerms, ex)
        end
    end
    X = length(AllTerms)==1 ? AllTerms[1] : Expr(:call, :+, AllTerms...)
    println("THE ANSWER IS ", X)
    X
end

expectedtrace(N::Integer, MyQ::Symbol, X::Expr)=expectation(N, MyQ, X, true)

if _HAVE_GSL

#Computes the Weingarten function for permutations
function WeingartenUnitary(N::Integer, P::Ptr{gsl_permutation})
    C = cycle_structure(P)
    WeingartenUnitary(N, C)
end
#Computes the Weingarten function for partitions
function WeingartenUnitary(N::Integer, P::Partition)
    println("Computing Weingarten function of ", P)
    n = sum(P)
    m = length(P)
    thesum = 0.0
    for irrep in @task part(n)
        #Character of the partition
        S = character(irrep, P)
        #Character of the identity divided by n!
        T = character_identity(irrep)
        #Denominator f_r(N) of (2.10)
        f = prod([factorial(BigInt(N + P[i] - i)) / factorial(BigInt(N - i))
                for i=1:m])
        thesum += S*T/(factorial(n)*f)
        println(irrep, " ", thesum)
    end
    float64(thesum)
end

end #_HAVE_GSL
