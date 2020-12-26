let num = [5, 8, 2, 9, 3]
num.sort()
num.push(1)

console.log(`${num}`)
console.log(`O vetor tem ${num.length} posições`)
console.log(`O primeiro valor do vetor é ${num[0]}`)
console.log(`chave do valor 9 = ${num.indexOf(9)}`) //se resultado for -1 posição não encontrada

/*for (let cont = 0; cont < num.length; cont++)
    console.log(`${num[cont]} -> `) */

for (let cont in num){
    console.log(num[cont])
}