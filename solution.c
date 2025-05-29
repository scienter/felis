#include "stdio.h"
#include "stdlib.h"
#include "math.h"

float whatBfield(float x,float y,float L,float B0)
{
  float result;
  
  if(x<=L) 
    result=B0;
  else
    result=0.0;
  return result;
}


void main(int argc, char *argv[])
{
  float v0,B,length,final,dt,BField;
  float x,y,t,vy,vx,oldVx,oldVy,midVx,midVy,v,alpha;
  if(argc<5)
  {
    printf("locus v0 B length finalTime dt\n");
    exit(0);
  }
  
  v0=atof(argv[1]);
  BField=atof(argv[2]);
  length=atof(argv[3]);
  final=atof(argv[4]);
  dt=atof(argv[5]);

  vy=0.0;
  vx=v0;
  t=0;
  while(t<=final)
  {
    B=whatBfield(x,y,length,BField);
    alpha=1.602e-19*B/1.672e-27*dt;
//printf("t=%g, vx=%g, vy=%g, alpha=%g\n",t,vx,vy,alpha);
    oldVy=vy;
    oldVx=vx;
    vy+=-alpha*vx;
    vx+=alpha*vy;
    midVx=0.5*(vx+oldVx);
    midVy=vy;
    x+=dt*midVx;
    y+=dt*midVy;
    v=sqrt(midVx*midVx+midVy*midVy);
    printf("%g %g %g\n",x,y,v);
    t+=dt;
  }
}
