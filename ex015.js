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

function multmatiz(m1, m2){
    let mr 
    for (let i = 0;   i < 4; i++){
        for(let j = 0; j < 4; j++){
            for(let z = 0; z < 4; z++)
                mr[i][j] += (m1[i][z] * m2[j][z])
        }        
    }
    return mr
}



