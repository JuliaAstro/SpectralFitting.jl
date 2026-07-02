function prettyfloat(f)
    if f == 0
        return "0.0"
    elseif f == Inf
        return "Inf"
    elseif f ≥ 1 && f < 1.0e5 && f - trunc(Int, f) < 1.0e-5
        return Printf.@sprintf("%.1f", f)
    else
        return Printf.@sprintf("%#.5g", f)
    end
end

function encapsulate(text)
    # drop everything after last new line
    n = findlast(==('\n'), text)
    s = strip(text[n:end]) == "" ? text[1:n] : text
    lines = split(s, "\n")
    out = map(enumerate(lines)) do (i, line)
        if i == 1
            "┌ " * line
        elseif i == lastindex(lines)
            "└ " * line
        else
            "│ " * line
        end
    end
    return join(out, "\n")
end

function indent(text, n)
    spacer = " "^n
    return replace(text, '\n' => "\n" * spacer)
end
