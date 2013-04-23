using GSL

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
    cycleindices = [1]
    for i=1:n-1 if HaveNewCycleAtPos[i] cycleindices = [cycleindices; i+1] end end
    cycleindices = [cycleindices; n+1]
    cyclestructure = [cycleindices[i+1]-cycleindices[i] for i=1:length(cycleindices)-1]
    @assert sum(cyclestructure) == n
    cyclestructure
end

#Returns a vector of indices (starting from 1 in the Julia convention)
data(P::Ptr{gsl_permutation}) = [convert(Int64, x)+1 for x in pointer_to_array(permutation_data(P), (convert(Int64, permutation_size(P)) ,))]


# In random matrix theory one often encounters expressions of the form
#
#X = trace(Q * A * Q' * B)
#
#where A and B are deterministic matrices with fixed numerical matrix entries and Q is a random matrix that does not have explicitly defined matrix elements. Instead, one takes an expectation over of expressions of this form and this "integrates out" the Qs to produce a numeric result.
#
#expectation(X) #= some number
#
#I'm curious how far one can implement, in Julia, a data type for Q and expectation() method for which expressions of this form can be written natively as code. This probably strays into the larger question of how much symbolic algebra can be written and handled in Julia, but I think this would be have real payoffs for making possible extremely elegant manipulations of random matrices that doesn't exist in any known computer language.
#
#Any thoughts on how feasible this is? Perhaps to begin with, it would be nice to have an expectation() that would work on arbitrary strings of products of matrices. Probably this would involve taking an input expression argument Would like to hear some stuff on t would be nice to start with just how one can write  arbitrary products of matrices My initial thoughts are to declare Q::MyRandomMatrix with the type MyRandomMatrix<:AbstractMatrix and to write an expectation(Ex::Expr) function that at first blush traverses the expression Ex, finds a * that takes at least one MyRandomMatrix in its arguments and applies the appropriate substitution.
#
#So... 1) How feasible is this and 2) How do I traverse the expression in the way described in the previous paragraph? The "Metaprogramming" manual is somewhat sparse on this other than saying .args[], but that doesn't quite help

function Expectation(X::Expr)
    if X.head != :call
        error(string("Unexpected type of expression: ", X.head))
    end
    
    n = length(X.args) - 1
    if n < 1 return eval(X) end #nothing to do Haar-wise
    
    if X.args[1] != :*
        error("Unexpected expression, only products supported")
    end

    # Parse expression involving products of matrices to extract the
    # positions of Haar matrices and their ctransposes
    Qidx=[] #Indices for Haar matrices
    Qpidx=[] #Indices for ctranspose of Haar matrices
    Others=[]
    MyQ=None
    for i=1:n
        thingy=X.args[i+1]
        if typeof(thingy)==Symbol
            matrixtype = typeof(eval(thingy))
            if matrixtype == HaarMatrix
                if MyQ==None MyQ=thingy end
                if MyQ == thingy
                    Qidx=[Qidx; i]
                else
                    println(string("only one instance of HaarMatrix supported, skipping the other guy ", thingy))
                end
            else
                Others = [Others; (thingy, i, i+1)]
            end
            println(i, ' ', thingy, "::", typeof(eval(thingy)))
        elseif typeof(thingy)==Expr
            println(i, ' ', thingy, "::Expr")
            if thingy.head==symbol('\'') && length(thingy.args)>=1 #Maybe this is a Q'
                if typeof(thingy.args[1])==Symbol && typeof(eval(thingy.args[1]))==HaarMatrix
                    println("Here is a Qtranspose")
                    Qpidx=[Qpidx; i]
                end
            end
        else
            error(string("Unexpected token ", thingy ," of type ", typeof(thingy)))
        end
    end

    if length(Qidx) == length(Qpidx) == 0 return eval(X) end #nothing to do Haar-wise
    println(MyQ, " is in places ", Qidx)
    println(MyQ, "' is in places ", Qpidx)

    n = length(Qidx)
    #If there are different Qs and Q's, the answer is a big fat 0
    if n != length(Qpidx) return zeros(size(eval(X.args[2]),1)...) end

    ##################################
    # Step 2. Enumerate permutations #
    ##################################
    AllTerms = :(+(A*B))
    delete!(AllTerms.args, 2)
    println("BEGINS WITH ", AllTerms)
    for sigma in @task permutations_in_Sn(n)
        for tau in @task permutations_in_Sn(n)
            sigma_inv=permutation_inverse(sigma)
            #Compose the permutations
            perm=compose(sigma_inv, tau)
            
            #Compute cycle structure, i.e. lengths of each cycle in the cycle
            #factorization of the permutatation
            cyclestruct = cycle_structure(perm)
            Qr = Int64[n+1 for n in Qidx]
            Qpr= Int64[Qpidx[n]+1 for n in data(tau)]
            #Consolidate deltas
            Deltas=Dict()
            for i=1:n
                Deltas[Qr[i]] = Qpr[i]
                Deltas[Qidx[i]] =  Qpidx[data(sigma)[i]]
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
            #Substitute
            ReindexedSymbols = []
            for (Symb, col_idx, row_idx) in Others
                new_col, new_row = col_idx, row_idx
                if has(Deltas, col_idx) new_col = Deltas[col_idx] end
                if has(Deltas, row_idx) new_row = Deltas[row_idx] end
                print(Symb,"(",new_col,",",new_row,") ")
                ReindexedSymbols = [ReindexedSymbols; (Symb, new_col, new_row)]
            end
            println()
            println(ReindexedSymbols)
            println()
            #Reconstruct expression
            println("START PARSING")
            Symb, left_idx, right_idx = delete!(ReindexedSymbols, length(ReindexedSymbols))
            Expression={{{Symb}, left_idx, right_idx}}
            println("E =", Expression)
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
                println("NEW EXPRBLOB: The term is =", ReindexedSymbols[1])
                    Symb, left_idx, right_idx = delete!(ReindexedSymbols, 1)
                    insert!(Expression, length(Expression)+1, {{Symb}, left_idx, right_idx})
                else #Found something
                println("The term is =", ReindexedSymbols[pop_idx])
                    Symb, col_idx, row_idx = delete!(ReindexedSymbols, pop_idx)
                    Term = do_transpose ? Expr(symbol("'"), Symb) : Symb
                        println(Expression[expr_idx])
                        println(Expression[expr_idx][1])
                    if is_left
                        println("insert left")
                        insert!(Expression[expr_idx][1], 1, Term)
                        Expression[expr_idx][2] = do_transpose ? row_idx : col_idx
                    else
                        println("insert right")
                        insert!(Expression[expr_idx][1],
                        length(Expression[expr_idx][1])+1, Term)
                        Expression[expr_idx][3] = do_transpose ? col_idx : row_idx
                    end
                end
                println("E=",Expression)
            end
            println("DONE PARSING TREE")

            # Evaluate closed cycles
            NewExpression=:(+(A*B))
            delete!(NewExpression.args, 2)
            println("HIHIHI", NewExpression.args)
            for ExprBlob in Expression
                if ExprBlob[2]==ExprBlob[3] #Have a cycle; this is a trace
                    println(ExprBlob[1])
                    ex=:(trace(A*B))
                    println(ex.args[2].args)
                    ex.args[2].args={:*; ExprBlob[1]...}
                    println(ex)
                    println(ex.args)
                    insert!(NewExpression.args, length(NewExpression.args)+1, ex)
                else #Not a cycle, regular chain of multiplications
                    ex=:(A*B)
                    ex.args={:*; ExprBlob[1]...}
                    insert!(NewExpression.args, length(NewExpression.args)+1, ex)
                end
            end
            println("Final expression is ", NewExpression)
            insert!(AllTerms.args, length(AllTerms.args)+1, NewExpression)
        end
    end
    println("THE ANSWER IS ", AllTerms)
    eval(AllTerms)
end


#Iterate over partitions of n in lexicographic order
function part(n::Integer)
    if n==0 produce([]) end
    if n<=0 return end
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



###############
# SANDBOX
N=5
A=randn(N,N)
B=randn(N,N)
type HaarMatrix
    beta::Real
end
Q = HaarMatrix(2)


println("Case 1")
println("E(A*B) = ", Expectation(:(A*B)))
println("A*B/N = ", A*B)

println("Case 2")
println("tr(A)*tr(B)/N = ", trace(A)*trace(B)/N)
#println("E(A*Q*B*Q') = ", Expectation(:(A*Q*B*Q')))
println("E.tr(A*Q*B*Q') = ", trace(Expectation(:(A*Q*B*Q'))))

#println("Case 3")
#println(Expectation(:(A*B)) == A*B)
#println("E(A*Q*B*Q'*A*Q*B*Q') = ", Expectation(:(A*Q*B*Q'*A*Q*B*Q')))

