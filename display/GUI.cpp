#include "GUI.h"
#include "model.h"
#include <sys/time.h>
#include <math.h>
#ifdef __linux__
#include  <GL/gl.h>
#include  <GL/glu.h>
#include  <GL/glut.h>
#elif _WIN32
#define GLUT_DISABLE_ATEXIT_HACK
#include <windows.h>
#include <my_include/gl.h>
#include <my_include/glu.h>
#include <my_include/glut.h>

#endif


int W_HEIGHT=900;
double w_x0=0.0;
double w_x1=1.0;

double w_y0=0.0;
double w_y1=1.0;


double w_z0=-10.0;
double w_z1=10.0;


int redr=0;
double color_scale=1.0;

double rx=0.0;
double ry=0.0;

double rx0=0.0;
double ry0=0.0;

double mx0=0.0;
double my0=0.0;

double mouse_x=0.0;
double mouse_y=0.0;

int rotate=0;
int draw_model=1;
extern model* model_1;
double scale =1.0;

int draw_dist=0;

double dist_cache[260][260];
double dist_fast_cache[260][260];
int N_cross=250;

double alpha_x=0.5;

double alpha_y=0.5;

double get_time(void) {
    struct timeval tv;
    struct timezone tz;
    gettimeofday(&tv, &tz);
    return ((double)(tv.tv_sec+tv.tv_usec*1.0e-6));
}

v3 get_color(double gval, double min, double max)
{
    const int nn=4;
    int i;
    double val;
    val=gval;
    if (val>max) val=max;
    if (val<min) val=min;

    v3 col_table[5];
    v3 res;

    col_table[0].m_x = 0.0; col_table[0].m_y = 0.0; col_table[0].m_z = 1.0;
    col_table[1].m_x = 0.0; col_table[1].m_y = 1.0; col_table[1].m_z = 1.0;
    col_table[2].m_x = 0.0; col_table[2].m_y = 1.0; col_table[2].m_z = 0.0;
    col_table[3].m_x = 1.0; col_table[3].m_y = 1.0; col_table[3].m_z = 0.0;
    col_table[4].m_x = 1.0; col_table[4].m_y = 0.0; col_table[4].m_z = 0.0;

    double alpha;
    if ((max-min) > 1e-35)
    {
        alpha=(val-min)/(max-min)*nn;
        i=(int)(alpha);
        alpha=alpha-i;
    }
    else
    {
        alpha=0.0;
        i=2;
    }
    res.m_x = col_table[i].m_x * (1 - alpha) + col_table[i+1].m_x * alpha;
    res.m_y = col_table[i].m_y * (1 - alpha) + col_table[i+1].m_y * alpha;
    res.m_z = col_table[i].m_z * (1 - alpha) + col_table[i+1].m_z * alpha;

    glColor3f(res.m_x,res.m_y,res.m_z);
    return res;

}
//draw call

void calc_distnce(double alpha,int N)
{
    for (int j=0;j<N+1;j++)
    {
        for (int i=0;i<N;i++)
        {
            v3 r;
            r.m_x=w_x0+i*(w_x1-w_x0)/N;
            r.m_y=w_y0*(1.0-alpha)+w_y1*alpha;//w_y0+j*(w_y1-w_y0)/N;
            r.m_z=w_z0+j*(w_z1-w_z0)/N;//w_z0*(1.0-alpha)+w_z1*alpha;
            dist_cache[i][j]=model_1->distP(r);
        }
    }
}

void calc_distnce_fast(double alpha,int N)
{
    for (int j=0;j<N+1;j++)
    {
        for (int i=0;i<N;i++)
        {
            v3 r;
            r.m_x=w_x0+i*(w_x1-w_x0)/N;
            r.m_y=w_y0+j*(w_y1-w_y0)/N;
            r.m_z=w_z0*(1.0-alpha)+w_z1*alpha;

            r.m_x=w_x0+i*(w_x1-w_x0)/N;
            r.m_y=w_y0*(1.0-alpha)+w_y1*alpha;//w_y0+j*(w_y1-w_y0)/N;
            r.m_z=w_z0+j*(w_z1-w_z0)/N;//w_z0*(1.0-alpha)+w_z1*alpha;

            dist_fast_cache[i][j]=model_1->distP_fast(r);
            //printf("%d  %d \n",i,j);
        }
    }
}



void draw_distance(double alpha, int axis, int N)
{
    for (int j=0;j<N;j++)
    {
        glBegin(GL_TRIANGLE_STRIP);
        for (int i=0;i<N;i++)
        {
            v3 r;

            //            r.m_x=w_x0+i*(w_x1-w_x0)/N;
            //          r.m_y=w_y0+j*(w_y1-w_y0)/N;
            //        r.m_z=w_z0*(1.0-alpha)+w_z1*alpha;

            r.m_x=w_x0+i*(w_x1-w_x0)/N;
            r.m_y=w_y0*(1.0-alpha)+w_y1*alpha;//w_y0+j*(w_y1-w_y0)/N;
            r.m_z=w_z0+j*(w_z1-w_z0)/N;//w_z0*(1.0-alpha)+w_z1*alpha;


            if (draw_dist==0){
                //get_color(dist_cache[i][j],0,scale*((w_x1-w_x0)+(w_y1-w_y0)+(w_z1-w_z0))/6.0);
                double d=dist_cache[i][j]*scale*((w_x1-w_x0)+(w_y1-w_y0)+(w_z1-w_z0))/6.0;
                glColor3f(-d,-d,d);
            }
            else
                get_color(dist_fast_cache[i][j],0,scale*((w_x1-w_x0)+(w_y1-w_y0)+(w_z1-w_z0))/6.0);
            glVertex3f(r.m_x,r.m_y,r.m_z);

            //            r.m_x=w_x0+(i)*(w_x1-w_x0)/N;
            //          r.m_y=w_y0+(j+1)*(w_y1-w_y0)/N;
            //        r.m_z=w_z0*(1.0-alpha)+w_z1*alpha;

            r.m_x=w_x0+i*(w_x1-w_x0)/N;
            r.m_y=w_y0*(1.0-alpha)+w_y1*alpha;//w_y0+j*(w_y1-w_y0)/N;
            r.m_z=w_z0+(j+1)*(w_z1-w_z0)/N;//w_z0*(1.0-alpha)+w_z1*alpha;


            if (draw_dist==0){
                //get_color(dist_cache[i][j+1],0,scale*((w_x1-w_x0)+(w_y1-w_y0)+(w_z1-w_z0))/6.0);
                double d=dist_cache[i][j]*scale*((w_x1-w_x0)+(w_y1-w_y0)+(w_z1-w_z0))/6.0;
                glColor3f(-d,-d,d);
            }
            else
                get_color(dist_fast_cache[i][j+1],0,scale*((w_x1-w_x0)+(w_y1-w_y0)+(w_z1-w_z0))/6.0);
            glVertex3f(r.m_x,r.m_y,r.m_z);
        }
        glEnd();
    }
}

void display(void)
{
    glEnable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    glRotatef(ry,1.0,0,0);
    glRotatef(rx,0.0,1.0,0);

    glColor3f(1,1,1);

    glBegin(GL_LINE_LOOP);
    glVertex3f(w_x0,w_y0,w_z0);
    glVertex3f(w_x0,w_y1,w_z0);
    glVertex3f(w_x1,w_y1,w_z0);
    glVertex3f(w_x1,w_y0,w_z0);
    glEnd();

    glBegin(GL_LINE_LOOP);
    glVertex3f(w_x0,w_y0,w_z1);
    glVertex3f(w_x0,w_y1,w_z1);
    glVertex3f(w_x1,w_y1,w_z1);
    glVertex3f(w_x1,w_y0,w_z1);
    glEnd();

    if (draw_model)
        model_1->draw();

    glPointSize(5.0);
    glBegin(GL_POINTS);
    glColor3f(0,1,0);
    double alpha =0.5;
    double x=w_x0*(1.0-alpha_x)+w_x1*alpha_x;
    double y=w_y0*(1.0-alpha)  +w_y1*alpha;//w_y0+j*(w_y1-w_y0)/N;
    double z=w_z0*(1.0-alpha_y)+w_z1*alpha_y;//w_z0*(1.0-alpha)+w_z1*alpha;
    glVertex3f(x,y,z);
    glEnd();

    v3 point;
    point.m_x=x;
    point.m_y=y;
    point.m_z=z;

    double min_dist =1e10;
    double sum_dot=0.0;//(0,0,0);
double delta=1e-5;
int min_ind=0;
for (int i=0; i<model_1->m_tris.size();i++)
{
    double dot_n;
    double l=model_1->m_tris[i].distP(point,dot_n);//*m_tris[i].getSign(point);
    if (fabs(min_dist)>=fabs(l))
    {
                   min_dist=l;
                   min_ind=i;
    }
}
glColor3f(1,0,0);
 model_1->m_tris[min_ind].draw();

/*
    for (int i=0; i<model_1->m_tris.size();i++)
    {
        double  dot_n;
        double l=m_tris[i].distP(point,dot_n);//*m_tris[i].getSign(point);
      //  if (fabs(l-min_dist)<1e-7)
        //    sum_dot+=dotProd(m_tris[i].normal,point-m_tris[i].m_p[0]);//*dot_n;

        if (fabs(min_dist)>=fabs(l)-delta)
        {
                       if  (fabs(l-min_dist)>=delta)
                sum_dot=dot_n;//dotProd(m_tris[i].normal,point-m_tris[i].m_p[0]);
            else
                sum_dot+=dot_n;//dotProd(m_tris[i].normal,point-m_tris[i].m_p[0]);

                       //sum_dot=dotProd(m_tris[i].normal,point-m_tris[i].m_p[0]);
                       min_dist=l;

        }

    }
  */
    draw_distance(0.5,2,N_cross);

    glutSwapBuffers();
    if (redr==1) glutPostRedisplay();
}

//mouse_move
void m_m(int x,int y) //mouse move
{
    if (rotate==1)
    {
        rx=rx0+0.5*(x-mx0);
        ry=ry0+0.5*(y-my0);
    }
    glutPostRedisplay();
}


//mouse down
void m_d(int button, int state,int x, int y)  //mouse down
{

    if (state==GLUT_UP)
    {
        rotate=0;
        rx0=rx;
        ry0=ry;
    }
    if (state==GLUT_DOWN)
    {
        rotate=1;
        mx0=x;
        my0=y;

    }
    double W_WIDTH=W_HEIGHT*(w_x1-w_x0)/(w_y1-w_y0);
    mouse_x=(1.0*x)/W_WIDTH;

    mouse_y=(W_HEIGHT-(1.0*y))/W_HEIGHT;

    glutPostRedisplay();
}



//keyboard handeling function
void kb(unsigned char key, int x, int y)
{
    if (key=='1')
    {
        draw_dist=0;
        printf("drawing distance from matrix \n");
    }

    if (key=='2')
    {
        draw_dist=1;
        printf("drawing distance from triangle centers \n");
    }

    if (key=='.')
    {
        scale*=1.1;
    }


    if (key==',')
    {
        scale/=1.1;
    }

    if (key==' ')
    {
        draw_model=!draw_model;
    }

    if (key=='w')
    {
        alpha_y=fmin(alpha_y+0.001,1.0);
    }
    if (key=='s')
    {
        alpha_y=fmax(alpha_y-0.001,0.0);
    }

    if (key=='a')
    {
        alpha_x=fmin(alpha_x+0.001,1.0);
    }
    if (key=='d')
    {
        alpha_x=fmax(alpha_x-0.001,0.0);
    }


    glutPostRedisplay();
}






void resize(int w, int h)
{
    glViewport(0, 0, w, h);
    glClearColor (0.0, 0.0, 0.0, 0.0);
    glColor3f(1.0, 1.0, 1.0);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    double h0,w0;
    w0=W_HEIGHT*(w_x1-w_x0)/(w_y1-w_y0);//w_x1*1.3-w_x0+(w_x1-w_x0)*0.3;
    h0=W_HEIGHT;//(w_y1-w_y0)*1.3;
    glOrtho(w_x0-(w_x1-w_x0)*0.3, w_x1+(w_x1-w_x0)*0.3, w_y0-(w_y1-w_y0)*0.3,w_y1 + (w_y1-w_y0)*0.3, w_z0-(w_z1-w_z0)*0.5, w_z1+(w_z1-w_z0)*0.5);
    glMatrixMode (GL_MODELVIEW);

}

GUI::GUI()
{

}

void GUI::init(int argc, char **argv)
{
    double dx,dy,dz;
    w_x0=model_1->w_x0;
    w_x1=model_1->w_x1;
    w_y0=model_1->w_y0;
    w_y1=model_1->w_y1;
    w_z0=model_1->w_z0;
    w_z1=model_1->w_z1;

    double t0=get_time();
    calc_distnce_fast(0.5,N_cross);
    double t1=get_time();
    printf("dist fast calc_time =%e \n",t1-t0);

    t0=get_time();
    calc_distnce(0.5,N_cross);
    t1=get_time();
    printf("dist full calc_time =%e \n",t1-t0);

    glutInit(&argc,argv);
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB |GLUT_DEPTH);
    glutInitWindowSize(W_HEIGHT*(w_x1-w_x0)/(w_y1-w_y0),W_HEIGHT);
    glutInitWindowPosition(0,0);
    glutCreateWindow("simple");
    glutDisplayFunc(display);
    glutReshapeFunc(resize);
    glutMotionFunc(m_m);
    glutMouseFunc(m_d);
    glutKeyboardFunc(kb);

    glClearColor (0.0, 0.0, 0.0, 0.0);
    glColor3f(1.0, 1.0, 1.0);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    glOrtho(w_x0-(w_x1-w_x0)*0.3, w_x1+(w_x1-w_x0)*0.3, w_y0-(w_y1-w_y0)*0.3,w_y1 + (w_y1-w_y0)*0.3, w_z0-(w_z1-w_z0)*1.5, w_z1+(w_z1-w_z0)*1.5);
    glMatrixMode (GL_MODELVIEW);

    glutMainLoop();
}


