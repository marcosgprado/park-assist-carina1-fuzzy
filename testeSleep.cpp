#include <cstdio>
#include <iostream>
#include <stdlib.h>


void estercar(float &e){

                e = e + 2.0; 
 		usleep(100000);
                

	     }
  





int main(){
float a = 3.0;

for(int i =0; i<100; i++){

estercar(a);

printf("- %.3f -", a);

}
return(0);

}
