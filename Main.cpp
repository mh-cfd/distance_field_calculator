
#include <stdio.h>
#include <stdlib.h>


#include  <math.h>
#include <time.h>
#include <string.h>

#include "display/GUI.h"
#include "model.h"

GUI* g_gui;
model* model_1;
int main(int argc, char** argv)
{
    char dir[1024],file_mod[1024];
    sprintf(dir,"%s",__FILE__);
    dir[strlen(dir)-9]=0; //hack to get source directory

    sprintf(file_mod,"%s/stl_files/dodeca_half_a.stl",dir);
   // sprintf(file_mod,"%s/stl_files/teapot.stl",dir);
    printf("hello world %s \n", file_mod);

    model_1= new model();

    model_1->load(file_mod);

  g_gui=new GUI();
  g_gui->init(argc,argv);
}
