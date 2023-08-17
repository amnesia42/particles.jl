function my_add!(a;b=1)
    a +=b 
    return nothing
end
function my_ops(a,c;b=1)
    a += c
    c += a
    println(b)
    my_add!(a;b=b)
    println(a,'\t', b,'\t', c)
    return nothing
end

a = 2
c = 1
my_ops(a,c)