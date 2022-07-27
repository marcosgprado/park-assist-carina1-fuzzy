#define _VFH_H 
 
#include "gps_angle.h"
#include "vales.h"
#include "utils.h"

double GPSangle;
int steering;
int vale_livre = 180;
double resultado;


/*------------------------------------------------------------------------------------------------
       				 Função para padronizar a saida de estercamento          						          
------------------------------------------------------------------------------------------------*/

double steering_fix(double steering){
    if(steering >= 30.0){
        return 30.0;
    }else if(steering <= -30.0){
        return -30.0;
    }else{
        return steering;
    }
}


/*------------------------------------------------------------------------------------------------
       				 Função para retornar o esterçamento a ser realizado          						          
------------------------------------------------------------------------------------------------*/ 
 
double vfh_steering(float *vfh, int x, int y, int yaw){

        GPSangle = calculaAngulo(x,y);

        steering = steering_fix(GPSangle - ((yaw * 180)/PI));
        
		/* Retorna qual o angulo a ser utilizado para alcançar o vale escolhido */
		vale_livre = encontraVales(vfh, steering);
		
       	/* Regras de esterçamento para desviar de obstáculos e alcançar o destino GPS */
  
  		if(vale_livre >= 180){
  	        if(steering >= 0){	
		        resultado = -20;
    		}else{
    			resultado = 20;
    		}
    	}else{
    		if(steering >= 0){	
		        resultado = 20;
    		}else{
    			resultado = 20;
    		}
    	}
    		
        vale_livre = 180;
        
        return resultado;
}





