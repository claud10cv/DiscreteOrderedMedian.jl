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
    return DOMPData(D, p, lambda)
end

