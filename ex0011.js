var agora = new Date
var hora = agora.getHours()
console.log(`Agora são exatamente ${hora} horas. `)

if(hora < 12 && hora > 5){
    console.log('Bom Dia!')
} else if(hora <= 18){
    console.log('Boa Tarde!')
}else if(hora < 24){
    console.log('Boa Noite!')
}else{
    console.log('Boa Madrugada')
}