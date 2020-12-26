function soma(n1=0, n2=0){ //começa predefinido como 0 para n retornar NAN                  eron flash ->
    return n1 + n2
}

console.log(soma(7))

let v = function(x){
    return x*2
}

console.log(v(5))

function fatorial(n){
    if(n == 1){
        return 1 
    } else {
        return n * fatorial (n - 1)
    }
}

