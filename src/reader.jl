function read_deleplanque(filename::String)::DOMPData
    p = n = 0
    D = []
    lambda = []
    open(filename) do f
        while !eof(f)
            line = readline(f)
            tok = split(line, [' ', ':']; keepempty = false)
            if tok[1] == "OPEN_FACILITIES"
                p = parse(Int64, tok[2])
            elseif tok[1] == "DIMENSION"
                n = parse(Int64, tok[2])
                D = zeros(Int64, n, n)
                lambda = zeros(Int64, n)
            elseif tok[1] == "COST_SECTION"
                for i in 1 : n
                    ctok = split(readline(f); keepempty = false)
                    for j in 1 : n
                        D[i, j] = parse(Int64, ctok[j])
                    end
                end
            elseif tok[1] == "LAMBDA_SECTION"
                for i in 1 : n
                    ctok = split(readline(f); keepempty = false)
                    k = parse(Int64, ctok[1])
                    l = parse(Int64, ctok[2])
                    while k > i
                        i += 1
                    end
                    lambda[i] = l
                end
            end
        end
    end
    uD = unique(sort(vec(D)))
    return DOMPData(D, p, lambda, uD)
end

function read_orlib(filename::String)::DOMPData
    nnodes_orlib = 0
    nedges_orlib = 0
    p = 0
    D_orlib = zeros(Int64, 0, 0)
    open(filename) do f
        let
            tok = split(readline(f), [' ', '\t', ',']; keepempty = false)
            nnodes_orlib = parse(Int64, tok[1])
            nedges_orlib = parse(Int64, tok[2])
            D_orlib = 1000000000000 * ones(Int64, nnodes_orlib, nnodes_orlib)
            for n in 1 : nnodes_orlib
                D_orlib[n, n] = 0
            end
            p = parse(Int64, tok[3])
        end
        for e in 1 : nedges_orlib
            tok = split(readline(f), [' ', '\t', ',']; keepempty = false)
            i = parse(Int64, tok[1])
            j = parse(Int64, tok[2])
            v = parse(Int64, tok[3])
            D_orlib[i, j] = D_orlib[j, i] = v
        end
    end
    floyd!(D_orlib)
    @assert(nnodes_orlib % 2 == 0, ("odd number of nodes!"))
    n = round(Int64, nnodes_orlib / 2)
    D = zeros(Int64, n, n)
    for i in 1 : n, j in 1 : n
        D[i, j] = D_orlib[i, j + n]
    end
    uD = unique(sort(vec(D)))
    DOMPData(D, p, rand(1 : 100, n), uD)
end

function read_beasley(filename::String)::DOMPData
    p = n = 0
    D = []
    lambda = []
    i = 0
    open(filename) do f
        line = readline(f)
        # println(line)
        tok = split(line, [' ', ':']; keepempty = false)
        n = parse(Int64, tok[1])
        p = parse(Int64, tok[2])
        D = zeros(Int64, n, n)
        lambda = zeros(Int64, n)
        # println(D)
        # println(n)
        while !eof(f)
            for i in 1 : n
                ctok = split(readline(f); keepempty = false)
                # println(ctok)
                for k in 1 : n
                    l = parse(Int64, ctok[k])
                    # println(l)
                    # println(i)
                    # println(k)
                    D[i, k] = l
                    # println(D[i, k])
                end
            end
              
        end
    end
    uD = unique(sort(vec(D)))
    # println(uD)
    return DOMPData(D, p, lambda, uD)
end

function read_deleplanque_ns(filename::String)::DOMPData
    p = n = 0
    D = []
    lambda = []
    open(filename) do f
        while !eof(f)
            line = readline(f)
            tok = split(line, [' ', ':']; keepempty = false)
            if tok[1] == "OPEN_FACILITIES"
                p = parse(Int64, tok[2])
            elseif tok[1] == "DIMENSION"
                n = parse(Int64, tok[2])
                D = zeros(Int64, n, n)
                lambda = zeros(Int64, n)
            elseif tok[1] == "COST_SECTION"
                for i in 1 : n
                    ctok = split(readline(f); keepempty = false)
                    for j in 1 : n
                        D[i, j] = parse(Int64, ctok[j])
                    end
                end
            elseif tok[1] == "LAMBDA_SECTION"
                
                ctok = split(readline(f); keepempty = false)
                k = parse(Int64, ctok[1])
                for i in 1 : n
                    l = parse(Int64, ctok[i+1])
                    lambda[i] = l
                end
            end
        end
    end
    uD = unique(sort(vec(D)))
    return DOMPData(D, p, lambda, uD)
end

function write_instance(data::DOMPData, filename::String)::Nothing
    name = split(filename, ['.', '/']; keepempty = false)[end - 1]
    nrows = size(data.D, 1)
    ncols = size(data.D, 2)
    open(filename, "w") do f
        println(f, "NAME : $name")
        println(f, "COMMENT : Cherkesly, Contardo & Gruson 2022")
        println(f, "TYPE : DOMP")
        println(f, "DIMENSION : $ncols")
        println(f, "DISPLAY_DATA_TYPE : NO_DISPLAY")
        println(f, "OPEN_FACILITIES : $(data.p)")
        println(f, "COST_SECTION")
        for i in 1 : nrows
            for j in 1 : ncols
                print(f, " $(data.D[i, j])")
            end
            println(f)
        end
        println(f, "LAMBDA_SECTION")
        for i in 1 : nrows
            println(f, "$i $(data.lambda[i])")
        end
    end
    return
end

function floyd!(D::Matrix{Int64})::Matrix{Int64}
    n = size(D, 1)
    for k in 1 : n
        for i in 1 : n
            for j in 1 : n
                D[i, j] = min(D[i, j], D[i, k] + D[k, j])
            end
        end
    end
    return D
end
