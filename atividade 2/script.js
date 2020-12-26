function verificar(){
    var data = new Date()
    var ano = data.getFullYear()
    var fano = document.getElementById('txtano')
    var res = document.querySelector('div#res')
    
    if(fano.value.length == 0 || Number (fano.value) > ano){
        window.alert('[ERRO] Verifique os dados e tente novamente!')
    
    } else {
        var fsex = document.getElementsByName('radsex')
        var idade = ano - Number(fano.value)
        var genero = ''
        var img = document.createElement('img')
        img.setAttribute('id', 'foto') //<img id = 'foto'> 
        
        if (fsex[0].checked){
            genero = 'Homem'
            //img.setAttribute('src', 'fotohomem.png')
        }
        else if (fsex[1].checked){
            genero = 'Mulher'
        }

        res.style.textAlign  = 'center'
        res.innerHTML = `Detectamos ${genero} com ${idade} anos`

        
    }

}