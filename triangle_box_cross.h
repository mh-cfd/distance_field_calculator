#ifndef TRIANGLE_BOX_CROSS_H
#define TRIANGLE_BOX_CROSS_H


#define FINDMINMAX(x0,x1,x2,min,max) \
  min = max = x0;   \
  if(x1<min) min=x1;\
  if(x1>max) max=x1;\
  if(x2<min) min=x2;\
  if(x2>max) max=x2;

int triBoxOverlap(float boxcenter[3],float boxhalfsize[3],float triverts[3][3]);

#endif // TRIANGLE_BOX_CROSS_H
