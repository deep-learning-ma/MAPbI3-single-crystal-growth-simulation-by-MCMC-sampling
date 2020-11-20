#include<stdio.h>
#include<math.h>
#include <string.h>
#include<stdlib.h>
#define a 16807
#define m 2147483647
#define n 500   //KMC模拟的步数
#define alpha 1
#define beta 1
#define ft 1
#define bt 1
#define N 10    //晶体生长的边界

int I(int x, int y, int h[N][N]){
    int i=0;
    if(h[x][y]=h[x-1][y]){
        i++;
    }
    if(h[x][y]=h[x+1][y]){
        i++;
    }
    if(h[x][y]=h[x][y-1]){
        i++;
    }
    if(h[x][y]=h[x][y+1]){
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

double SurfaceDiffusion(int x, int y, int h[N][N]){//direction=1,2,3,4分别代表向四个方向的表面扩散
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
    else{
        fij=bt*Kplus(x,y, h)*Kminus(x,y+1, h);
    }
    return fij;
}

double Kxy(int x, int y, int h[N][N]){
    return Kplus(x,y, h)+Kminus(x,y, h)+SurfaceDiffusion(x,y, h);
}

double Kx(double x, int h[N][N]) {
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
	//freopen("test.dat","w",stdout);  //以只写访问方式打开一个txt文件
	int q=0,r=0,i=0,h[N][N]={0}, x, y, R;
    long z[n+1]={1};
    double xi[n+1], Xi[n+1]={0};
	for(i=0;i<n;i++){             //schrage方法生成MC模拟用的随机数
        q=floor(m/a);
		r=m%a;
        z[i+1]=a*(z[i]%q)-r*floor(z[i]/q);
		if (z[i+1] < 0) z[i+1] += m;//取为正值
	xi[i]=z[i]/(double)m;          //归一化
	//printf("%.15f\n",xi[i]);      //调试
	}

	for(i=1;i<n;i++){
            for(x=1;x<N-1;x++){
                for(y=1;y<N-1;y++){
                    I(x,y,h);
                    Kplus(x,y, h);
                    Kminus(x,y, h);
                }
            }
            for(x=1;x<N-1;x++){
                for(y=1;y<N-1;y++){
                    SurfaceDiffusion(x,y, h);
                    Kxy(x,y, h);
                    //printf("%.15f\n",Kxy(x,y,h));  //调试
                }
                Kx(x, h);
                //printf("%.15f\n",Kx(x,h));  //调试
            }
            K(h);
            //printf("%.15f\n",K(h));    //调试
	    Xi[i]=xi[i]*K(h);
	    //printf("%lf\n",Xi[i]);       //调试

        for(x=1;x<N-1;x++){         //随机选择合适的(x,y)作为晶体的生长点
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
        //printf("%d\n",h[x][y]);   //调试
	}
	else if(Xi[i]< Kplus(x,y, h)+Kminus(x,y, h) && h[x][y]>0){
            h[x][y]-=1;
            printf("%d\n",h[x][y]);
        }
        else{
            R=rand()%4+1; //产生1,2,3,4中的随机数
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
        for(y=1;y<N-1;y++){
            for(x=1;x<N-1;x++){
                printf("%d\t ",h[x][y]);
            }
            printf("\n");
        }
}
