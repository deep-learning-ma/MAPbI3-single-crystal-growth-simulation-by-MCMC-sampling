#include<stdio.h>
#include<math.h>
#include <string.h>
#include<stdlib.h>
#include <time.h>
#define a 16807
#define m 2147483647
#define n 40000   //Number of KMC steps
#define alpha 1
#define beta 1
#define ft 1
#define bt 0.1
#define ct 0.1
#define N 30    //boundary of crystal frowth

int I(int x, int y, int h[N][N]){
    int i=0;
    if(h[x][y]==h[x-1][y]){
        i++;
    }
    if(h[x][y]==h[x+1][y]){
        i++;
    }
    if(h[x][y]==h[x][y-1]){
        i++;
    }
    if(h[x][y]==h[x][y+1]){
        i++;
    }
    return i;
}

double Kplus(int x, int y, int h[N][N]){
    return ft*exp(-alpha*(2-I(x,y, h))/4.0+beta/2.0);
}

double Kminus(int x, int y, int h[N][N]){
    return ft*exp(alpha*(2-I(x,y,h))/4.0+beta/2.0);
}

double SurfaceDiffusion(int x, int y, int h[N][N]){
    
    //direction=1,2,3,4 represent particles converge towards the center, separately
    
    double fij;
    int direction;
    direction=rand()%4+1;
    if(direction==1){
        fij=bt*Kplus(x, y, h)*Kminus(x-1,y, h);
    }
    if(direction==2){
        fij=bt*Kplus(x,y, h)*Kminus(x+1,y, h);
    }
    if(direction==3){
        fij=bt*Kplus(x,y, h)*Kminus(x,y-1, h);
    }
    if(direction==4){
        fij=bt*Kplus(x,y, h)*Kminus(x,y+1, h);
    }
    return fij;
}

double SurfaceSpread(int x, int y, int h[N][N]){
    double Sij;
    int direction;
    direction=rand()%4+1;
    if(direction==1){
        Sij=ct*Kminus(x, y, h)*Kplus(x-1,y, h);
    }
    if(direction==2){
        Sij=ct*Kminus(x,y, h)*Kplus(x+1,y, h);
    }
    if(direction==3){
        Sij=ct*Kminus(x,y, h)*Kplus(x,y-1, h);
    }
    if(direction==4){
        Sij=ct*Kminus(x,y, h)*Kplus(x,y+1, h);
    }
    return Sij;
}

double Kxy(int x, int y, int h[N][N]){
    return Kplus(x,y, h)+Kminus(x,y, h)+SurfaceDiffusion(x,y, h)+SurfaceSpread(x,y,h);
}

double Kx(int x, int h[N][N]) {
    double s=0;
    int y;
    for(y=0;y<N;y++){
            s+=Kxy(x,y, h);
        }
    return s;
}

double K(int h[N][N]) {
    double s=0;
    int x, y;
    for(x=0;x<N;x++){
        for(y=0;y<N;y++){
            s+=Kxy(x,y, h);
        }
    }
    return s;
}
int main(){
    // srand((unsigned)time( NULL ));
	freopen("test.dat","w",stdout);  //open a .txt document with write only access
	int q=0,r=0,i=0,h[N][N]={0}, x, y, R,P;
    long z[n+1]={1};
    double xi[n+1], Xi[n+1]={0};
	for(i=0;i<n;i++){             //RNG, with schrage method, produce the random number
        q=floor(m/a);
		r=m%a;
        z[i+1]=a*(z[i]%q)-r*floor(z[i]/q);
		if (z[i+1] < 0) z[i+1] += m;//make sure it's positive
	xi[i]=z[i]/(double)m;          //normalize
	//printf("%.15f\n",xi[i]);      //debugging
	}

	for(i=1;i<n;i++){
        // for(x=1;x<N-1;x++){
        //     for(y=1;y<N-1;y++){
        //         I(x,y,h);
        //         Kplus(x,y, h);
        //         Kminus(x,y, h);
        //     }
        // }
        // for(x=1;x<N-1;x++){
        //     for(y=1;y<N-1;y++){
        //         SurfaceDiffusion(x,y, h);
        //         Kxy(x,y, h);
        //         //printf("%.15f\n",Kxy(x,y,h));  //debugging
        //     }
        //     Kx(x, h);
        //     //printf("%.15f\n",Kx(x,h));  //debugging
        // }
        //K(h);
        //printf("%.15f\n",K(h));    //debugging
	    Xi[i]=xi[i]*K(h);
	    //printf("%lf\n",Xi[i]);       //debugging

        for(x=1;x<N-1;x++){         //randomly choose the suitable (x,y) as the growth sites of crystal
            if(Xi[i]>=Kx(x, h)){
            Xi[i]-=Kx(x, h);
        }
        else{
            for(y=1;y<N-1;y++){
                if(Xi[i]>Kxy(x,y, h)){
                    Xi[i]-=Kxy(x,y, h);
                }
                else{
                    break;
                }
            }
        }

	    if(Xi[i]<Kplus(x,y, h)){
            h[x][y]+=1;
            //printf("%d\n",h[x][y]);   //debugging
            //printf("%f\n",Kminus(x,y,h));
	    }
	    else if(Xi[i]< Kplus(x,y, h)+Kminus(x,y, h) && h[x][y]>0){
            h[x][y]-=1;
            //printf("%d\n",h[x][y]);
        }
        else if(Xi[i]<Kplus(x,y,h)+Kminus(x,y,h)+SurfaceSpread(x,y,h)){
            P=rand()%4+1;
            if(P==1 && h[x][y]>0){
                h[x][y]-=1;
                h[x-1][y]+=1;
            }
            if(P==2 && h[x][y]>0){
                h[x][y]-=1;
                h[x+1][y]+=1;
            }
            if(P==3 && h[x][y]>0){
                h[x][y]-=1;
                h[x][y-1]+=1;
            }
            if(P==4 && h[x][y]>0){
                h[x][y]-=1;
                h[x][y+1]+=1;
            }
            else{
                continue;
            }
        }
        else{
            R=rand()%4+1; //generate random number from (1,2,3,4)
            if(R==1 && h[x-1][y]>0){
                h[x][y]+=1;
                h[x-1][y]-=1;
            }
            if(R==2 && h[x+1][y]>0){
                h[x][y]+=1;
                h[x+1][y]-=1;
            }
            if(R==3 && h[x][y-1]>0){
                h[x][y]+=1;
                h[x][y-1]-=1;
            }
            if(R==4 && h[x][y+1]>0){
                h[x][y]+=1;
                h[x][y+1]-=1;
            }
            else{
                continue;
            }
        }
            continue;
        }
    }
    /*
    for(y=1;y<N-1;y++){         //print the shape of crystal while the growth is over
        for(x=1;x<N-1;x++){
            printf("%4d",h[x][y]);
        }
        printf("\n");
    }
    */
    for(y=1;y<N-1;y++){   //print the coordinates of every lattice sites
        for(x=1;x<N-1;x++){
            printf("%4d %4d %4d \n",x,y,h[x][y]);
        }
    }
}
















