#define _GPS_ANGLE_H

#include "utils.h"

/*--------------------------------calculaAngulo-------------------------------*/

/* Pega a posição x e y do robô e subtrai do destino GPS. Calcula Arco Tangente
do resultado e tranforma de radianos em graus. De acordo com o quadrante, a
tranformação é feita. */

double calculaAngulo(double x, double y){
    double rad, angle;
    double dx = GOAL_X - x;
    double dy = GOAL_Y - y;

    rad = atan((dy)/(dx));
    if(dx < 0){
        if(dy > 0){
            angle = ((rad * 180) / PI) + 180;
        }else{
            angle = ((rad * 180) / PI) - 180;
        }
    }else{
            angle = ((rad * 180) / PI);
    }
    return angle;
}
/* Retorna o ângulo de destino do robô de acordo com o ponto GPS */

