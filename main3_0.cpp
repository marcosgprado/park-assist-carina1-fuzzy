// main.cpp


#include "utils.h"
#include "fis.h"
#include <cstdio>
#include <cstdarg>
#include <libplayerc++/playerc++.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "vfh.h"

using namespace PlayerCc;
using namespace std;

#define dist_min01 3.5
#define dist_min02 2.0
#define dist_min03 2.5
#define x_goal 32
#define y_goal 0  
#define deslocamento 40

int erroX=0, erroY=0, erroPeriodo = 0, estado = 0; 
float senDir0,senDir1, senDir2, senEsq0,senEsq1,senEsq2, senFrt; //armazena a média de conjuntos de feixes
float senDir00,senDir01, senDir02, senEsq00,senEsq01,senEsq02, senFrt00; // armazena o valor de um unico feixe


void erroGPS(float ampErro){
    if(erroPeriodo == 0){
        srand(time(NULL));
        
        while(erroPeriodo == 0)
            erroPeriodo = rand() % 15;
            
        // generates a randomic number between -maxError and maxError
        float ex = rand() % (int)(ampErro * 2 * 1000);
        ex /= 1000;
        erroX = ex - ampErro;

        // generates a randomic number between -maxError and maxError
        float ey = rand() % (int)(ampErro * 2 * 1000);
        ey /= 1000;
        erroY = ey - ampErro;
    }

    erroPeriodo --;
}

void estercar(float &estercamentoReal, float estercamentoDesejado){

         if((estercamentoReal <= 27) && (estercamentoReal >= -27)){

             if(estercamentoReal > estercamentoDesejado){

                estercamentoReal = estercamentoReal - 3.0;
 		usleep(150000);
             } 
  	     else{
     		estercamentoReal = estercamentoReal + 3.0;
		usleep(150000);
	     }
  	  }

         
}

///*
void readLaserFeixe(LaserProxy &lp){
         senDir0 = 0;
         senDir1 = 0;
         senDir2 = 0;
         senEsq0 = 0;
         senEsq1 = 0;
         senEsq2 = 0;
         senFrt = 0;
         for(int i=0;i<5;i++){
           
            senDir0 += lp.GetRange(i);
            senDir1 += lp.GetRange(55+i);
            senDir2 += lp.GetRange(115+i);
            senEsq0 += lp.GetRange(355+i);
            senEsq1 += lp.GetRange(295+i);
            senEsq2 += lp.GetRange(235+i);
            senFrt += lp.GetRange(178+i);
         
         }
         
            senDir0 = senDir0/5;
            senDir1 = senDir1/5;
            senDir2 = senDir2/5;
            senEsq0 = senEsq0/5;
            senEsq1 = senEsq1/5;
            senEsq2 = senEsq2/5;
            senFrt = senFrt/5;
}



void readLaser(LaserProxy &lp){
         senDir00 = 0;
         senDir01 = 0;
         senDir02 = 0;
         senEsq00 = 0;
         senEsq01 = 0;
         senEsq02 = 0;
         senFrt00 = 0;
                    
            senDir00 = lp.GetRange(0);
            senDir01 = lp.GetRange(60);
            senDir02 = lp.GetRange(120);
            senEsq00 = lp.GetRange(360);
            senEsq01 = lp.GetRange(300);
            senEsq02 = lp.GetRange(240);
            senFrt00 = lp.GetRange(180);
}

bool fim(float x, float y){

 if((x < (x_goal + .75)) && (x > (x_goal - .75)) && (y < (y_goal + .75)) && (y > (y_goal - .75)))
    return true;
 else
    return false;

/*
     if (senDir2 <= 0.4 && senEsq2 <= 0.4 && senFrt00 <= 0.3)
        return true;
     else
        return false;   
*/
}


bool limiteMax(void){

     if (senDir0 < 1.2 || senDir1 < 1. || senDir2 < 0.5 || senEsq0 < 1.2 || senEsq1 < 1. || senEsq2 < 0.5 || senFrt00 < 0.4)
        return true;
     else
        return false;   

}



bool limiteMed(void){

     if (senDir00 > 1. && senDir01 > 1. && senDir02 > 0.7 && senEsq00 > 1. && senEsq01 > 1. && senEsq02 > 0.7 && senFrt00 > 0.5)
        return true;
     else
        return false;   

}

bool livre(void){

     if ((senDir00 > dist_min01) && (senDir01 > dist_min01) && (senDir02 > dist_min01) && (senEsq00 > dist_min01) && (senEsq01 > dist_min01) && (senEsq02 > dist_min01) && (senFrt00 > dist_min01))
        return true;
     else
        return false;   

}


bool livre01(void){

     if ((senDir00 > dist_min03) && (senDir01 > dist_min03) && (senDir02 > dist_min03) && (senEsq00 > dist_min03) && (senEsq01 > dist_min03) && (senEsq02 > dist_min03) && (senFrt00 > dist_min03))
        return true;
     else
        return false;   

}


bool direitaLivre(void){

     if (senDir00 > dist_min02 && senDir01 > dist_min02 && senDir02 > dist_min02 && senFrt00 > dist_min02)
        return true;
     else
        return false;   

}


bool esquerdaLivre(void){

     if (senEsq00 > dist_min02 && senEsq01 > dist_min02 && senEsq02 > dist_min02 && senFrt00 > dist_min02)
        return true;
     else
        return false;   

}


bool frenteLivre01(void){

     if (senEsq02 > 1.4 && senDir02 > 1.4 && senFrt00 > 1)
        return true;
     else
        return false;   

}


bool frenteLivre02(void){

     if (senEsq02 > 3.4 && senDir02 > 3.4 && senFrt00 > 3)
        return true;
     else
        return false;   

}

//*/
int main(int argc, char **argv)
{

        Fis_Node *fis;
	int i, j;
	int debug = 0;

	double **dataMatrix, **fisMatrix, **outputMatrix;
	char *fis_file, *data_file;
	int data_row_n, data_col_n, fis_row_n, fis_col_n;
	float mediaDir, mediaEsq;
        float posicao_y,posicao_x, angulo, estercamento = 0, estercamentoRe = 0, estercamentoDesejado = 0;
   	float *vfh;
	int flag = 0;
   
   try
   {	
	
	

        PlayerClient robot(PlayerCc::PLAYER_HOSTNAME, PlayerCc::PLAYER_PORTNUM);

        SimulationProxy sp(&robot, 0);

        Position2dProxy pp(&robot,0);

        LaserProxy lp00(&robot, 0);

        lp00.RequestGeom();

        player_pose3d_t laserPose = lp00.GetPose();

        pp.SetMotorEnable(true);

        double x,y,z,roll,pitch,yaw, time;

        
        data_file = "fis_in";
        fis_file = "regrasFuzzy01.fis";
    
     
        for(;;){
     
        robot.Read();   


	/* obtain data matrix and FIS matrix */
	dataMatrix = fis->returnDataMatrix(data_file, &data_row_n, &data_col_n);  
	fisMatrix = fis->returnFismatrix(fis_file, &fis_row_n, &fis_col_n);

	/* build FIS data structure */
	fis = (Fis_Node *)fisCalloc(1, sizeof(Fis_Node));
	fis->fisBuildFisNode(fis, fisMatrix, fis_col_n, MF_POINT_N);

	/* error checking */
	if (data_col_n < fis->in_n) {
		printf("Given FIS is a %d-input %d-output system.\n",
			fis->in_n, fis->out_n);
		printf("Given data file does not have enough input entries.\n");
		fis->fisFreeMatrix((void **)dataMatrix, data_row_n);
		fis->fisFreeMatrix((void **)fisMatrix, fis_row_n);
		fis->fisFreeFisNode(fis);
		fisError("Exiting ...");
	}

	/* debugging */
	if (debug)
		fis->fisPrintData(fis);


        

	/* create output matrix */
	outputMatrix = (double **)fis->fisCreateMatrix(data_row_n, fis->out_n, sizeof(double));



	/* evaluate FIS on each input vector */
	
	dataMatrix[0][0] = posicao_y;
        dataMatrix[0][1] = angulo;
        fis->getFisOutput(dataMatrix[0], fis, outputMatrix[0]);



       // printf("\e[H\e[2J");

        sp.GetPose3d((char*)"Carina_model", x, y, z, roll, pitch, yaw, time);
           
 //      printf("Robot Pose3d: XYZ(m)[%.3f %.3f %.3f] angulo(rad):[%.3f] laser[dir:%.3f - esq:%.3f - frt:%.3f] estercamento: %.12f \n",x, y, z, yaw, lp00.GetRange(0),lp00.GetRange(360),lp00.GetRange(180),outputMatrix[0][0]);

 printf("angulo(rad):[%.3f] [ esq0:%.3f  esq1:%.3f esq2:%.3f frt:%.3f dir2:%.3f dir1:%.3f dir0:%.3f steering: %.3f flag: %d]  \n",yaw, senEsq00,senEsq01 , senEsq02, senFrt00, senDir02, senDir01, senDir00, estercamento, flag);



         estercamento = outputMatrix[0][0];
         // estercamentoDesejado =  outputMatrix[0][0];
         // estercar(estercamento, estercamentoDesejado); 
         erroGPS(3);
         posicao_y = y - y_goal; 
         //posicao =+ erroY;//acrescimo do erro do gps
         angulo = yaw;
         
         //captura dos daddos do laser
         //sete conjunto de feixes de laser distribuidos simetricamente.
                  
         readLaserFeixe(lp00);
         readLaser(lp00);
         
         mediaDir = (senDir1 + senDir2 + senDir0)/3;
         mediaEsq = (senEsq1 + senEsq2 + senEsq0)/3;
         //flag = 0;
        ///* 
         if(!fim(x,y)){
           
             // Se todos os sete feixes capturarem distancia maiores que 3.5m.
             if(livre()){ 
              flag = 0;
                //Velocidade constante e esterção calculada pelo sistema fuzzy.
                pp.SetSpeed(0.2, estercamento);
             }
             
             else{
           ///*   
flag = 5;

                   // Veículo em frente à vaga. 
                   if(frenteLivre01() && (y < (y_goal+4.)) && (y > (y_goal-4.)) && (x > (x_goal-7.)) && (x < (x_goal+7.))){   
                   	pp.SetSpeed(0.1, estercamento);
                   flag = 2;
                   
                   }
            
                   //verifica se veículo está muito perto de um obstáculo e  está perto  da vaga.
                   if(limiteMax() && (x > (x_goal-4))&& (x < (x_goal+4)) && (y < (y_goal+4)) && (y > (y_goal-4)) /*&& senDir00 > 1.5 && senEsq00 > 1.5 && senFrt00 < 0.7*/){
                        flag = 3;

                   	if(estercamento > -2 && estercamento < 2){
				if(mediaDir > mediaEsq){                       
                                	robot.Read(); 
	                            	readLaser(lp00);
					readLaserFeixe(lp00);
                                        mediaDir = (senDir1 + senDir2 + senDir0)/3;
         				mediaEsq = (senEsq1 + senEsq2 + senEsq0)/3;
	                            	pp.SetSpeed(-0.1, -20);
                                }
                                else{                       
                                	robot.Read(); 
	                            	readLaser(lp00);
					readLaserFeixe(lp00);
                                        mediaDir = (senDir1 + senDir2 + senDir0)/3;
         				mediaEsq = (senEsq1 + senEsq2 + senEsq0)/3;
	                            	pp.SetSpeed(-0.1, 20);
                                }
 
                        }
                        else{
	                       while(!limiteMed() || fim(x, y)){
	                            robot.Read(); 
	                            readLaser(lp00);
	                            pp.SetSpeed(-0.1, -estercamento);         
                               }
                	}
                        
                   }
                   
                   //verifica se veículo está muito perto de um obstáculo e está dentro da vaga.                  
                /*   else if(limiteMax() && x < (x_goal-1)){
                   
                       if(mediaDir > mediaEsq)
                          estercamentoRe = -30;
                       else
                          estercamentoRe = 30;  
                      
                       while(!limiteMed()){
                            flag = 4;
                            robot.Read(); 
                            readLaser(lp00);
                            pp.SetSpeed(-0.1, estercamentoRe);         
                        }
                
                   }  */ 
//*/
                   //Desvio de obstaculos por VFH
                   else if(!((y < (y_goal+4.)) && (y > (y_goal-4.)) && (x > (x_goal-7.)) && (x < (x_goal+7.)))) {
                         
				//pp.SetSpeed(0,0);
				//return(0);
                      
			/* aloca memória para o vetor VFH */
        		vfh = (float *) malloc(sizeof(float) * QTD_LEITURAS_LASER);
			/* Preenche o vetor vfh com a leitura do laser, gravando os valores * 10,
		         para aumentar a visualização do histograma */
			flag = 1;
		        for (int i = 0; i <= QTD_LEITURAS_LASER; i++){
		            vfh[i] = (lp00.GetRange(i)*10);
		        }
			estercamento = vfh_steering(vfh, x, y, yaw);
			free(vfh);
                        //estercar(estercamento, estercamento);

                        if (y < 0){
                         estercamento = -estercamento;

                        }
\
			pp.SetSpeed(0.1, estercamento); 

                   }
             
      /*
               

                   //verifica o lado que se encontra o obstaculos(quanto maior o valor da variavel, maior chance do lado, corespondente à variável, está livre)
                   if(mediaDir >= mediaEsq){
                      
                      if(direitaLivre()){
                      
                         
                          }                      
                         if((senEsq01 < 3.5 || senEsq00 < 3.5) && (posicao_y < -3.5 || posicao_y > 3.5) )       
                           pp.SetSpeed(0.1, 30);
                           
                           
                          //Se o veículo estiver perpendicular e próximo da vaga. 
                         if(posicao_y > y_goal && posicao_y <= (y_goal + 2) && senEsq00 > 4 && (posicao_x >= x_goal-7)) 
                                pp.SetSpeed(0.1, 30);
                                 
                          //Desviando de obstáulos.
                         if((senEsq02 < 3.5 || senFrt00 < 3 || senEsq00 < 3.5 || senEsq01 < 3.5) && (posicao_x < x_goal-7)){ 
                         
                        //pp.SetSpeed(0.1, 30);
                       
                          //   while(!frenteLivre02()){
                           //         robot.Read(); 
                           //         readLaser(lp00);
                           //         pp.SetSpeed(0.1, 30);
                       //      }  
                       
                                               
                         }
                   
                   }
                   else{
                   
                     if(esquerdaLivre()){
                      
                          
                                
                      
                         if((senDir01 < 3.5 || senDir00 < 3.5) && (posicao_y < -3.5 || posicao_y > 3.5) )
                           pp.SetSpeed(0.1, -30);
                           
                           
                           //Se o veículo estiver perpendicular e próximo da vaga.
                         if(posicao_y < y_goal && posicao_y >= (y_goal - 2) && senDir00 > 4 && (posicao_x >= x_goal-7)) 
                                pp.SetSpeed(0.1, -30); 
                                
                                
                                
                          if((senDir02 < 3.5 || senDir00 < 3 || senDir01 < 3.5) && (posicao_x < x_goal-7) ){ 
                          // pp.SetSpeed(0.1, -30);
                            
                           //  while(!frenteLivre02()){
                           //         robot.Read(); 
                            //        readLaser(lp00);
                             //       pp.SetSpeed(0.1, -30);
                            // }
                         
                          }
                 
                    }
                   
                   
                   }
             
           */      
                 
             } 
         
         }
         //}
         else{
         pp.SetSpeed(0,0);
         return 0;
         
         }
         
         
	/* clean up memory */
	fis->fisFreeFisNode(fis);
	fis->fisFreeMatrix((void **)dataMatrix, data_row_n);
	fis->fisFreeMatrix((void **)fisMatrix, fis_row_n);
	fis->fisFreeMatrix((void **)outputMatrix, data_row_n);
	
	}
    }
    catch (PlayerCc::PlayerError e)
    {
        std::cerr << "Erro encontrado: " << e << std::endl;
        return -1;
    }
    return 0;
}
