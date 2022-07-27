#define _VALES_H

#include "utils.h"
#include "limits.h"

int vetor_vales[20];
int tam_vetor_vales = 0;
void armazena_vales(int, int);
int escolhe_vale(int);
int center;

/*------------------------encontra vales--------------------------------------*/

int encontraVales(float *vfh, int steering){
    int flag = 0;
    int inicio = 0;
    int fim = 0;
    int *obst;
    int vale = 180;

    obst = (int *) malloc(sizeof(int) * QTD_LEITURAS_LASER);
    /*cria um vetor com 1 se há obstáculos e 0 se está livre,
     respeitando o limiar*/
    for(unsigned int i = 0; i < QTD_LEITURAS_LASER; i++){
		obst[i] = 0;
        if(vfh[i] < LIMIAR){
           obst[i] = 1;
        }
    }
    /* encontra inicio e fim dos vales */
    for(int i = 0; i < QTD_LEITURAS_LASER; i++){
        if((obst[i] == 0)&&(flag == 0)){ /* entrou no vale */
            inicio = i;
            flag = 1;
        }
        /* testa se saiu do vale */
        if((obst[i] == 1)&&(flag == 1)){
            fim = i - 1;
            flag = 0;
        }
        armazena_vales(inicio, fim);
    }
        
    /* Chama a função que escolhe o melhor vale, de acordo com o destino GPS */
    if(tam_vetor_vales > 0){
        vale = escolhe_vale(steering);
        //printf("\n Vale: %d", vale);
    }
    
    free(obst);
    return vale;
}


/* -------------------------------------------------------------------------- */
void armazena_vales(int inicio, int fim){
    /* testa se o carro passa no vale e atribui o centro ao vetor de vales */
	if(fim - inicio >= TAM_MIN_VALE){
       	center = (fim + inicio) / 2;
        vetor_vales[tam_vetor_vales] = center;
        tam_vetor_vales++;
    }
}

/* -------------------------------------------------------------------------- */
int escolhe_vale(int estercamento){
    int menor = vetor_vales[0];
    int diferenca;
    int resultado = 180;
    for(int i = 1; i < tam_vetor_vales; i++){
        diferenca = (abs(vetor_vales[i]) - abs(estercamento));
	    if(diferenca < menor){
            menor = diferenca;
            resultado = vetor_vales[i];
	    }
    }
    tam_vetor_vales = 0;
    return resultado;
}

