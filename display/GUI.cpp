#include "GUI.h"
#include "model.h"
#ifdef __linux__
#include  <GL/gl.h>
#include  <GL/glu.h>
#include  <GL/glut.h>
#elif _WIN32
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

extern model* model_1;


int draw_dist=0;

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

void draw_distance(double alpha, int axis, int N)
{
    for (int j=0;j<N;j++)
    {
        glBegin(GL_TRIANGLE_STRIP);
        for (int i=0;i<N;i++)
        {
            v3 r;

            r.m_x=w_x0+i*(w_x1-w_x0)/N;
            r.m_y=w_y0+j*(w_y1-w_y0)/N;
            r.m_z=w_z0*(1.0-alpha)+w_z1*alpha;
            if (draw_dist==0)
                get_color(model_1->distP(r),0,((w_x1-w_x0)+(w_y1-w_y0)+(w_z1-w_z0))/6.0);
            else
                get_color(model_1->distP_naive(r),0,((w_x1-w_x0)+(w_y1-w_y0)+(w_z1-w_z0))/6.0);
            glVertex3f(r.m_x,r.m_y,r.m_z);

            r.m_x=w_x0+(i)*(w_x1-w_x0)/N;
            r.m_y=w_y0+(j+1)*(w_y1-w_y0)/N;
            r.m_z=w_z0*(1.0-alpha)+w_z1*alpha;
            if (draw_dist==0)
                get_color(model_1->distP(r),0,((w_x1-w_x0)+(w_y1-w_y0)+(w_z1-w_z0))/6.0);
            else
                get_color(model_1->distP_naive(r),0,((w_x1-w_x0)+(w_y1-w_y0)+(w_z1-w_z0))/6.0);
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


       model_1->draw();
       draw_distance(0.5,2,50);
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
    w_x0=model_1->m_xMin;
    w_x1=model_1->m_xMax;
    w_y0=model_1->m_yMin;
    w_y1=model_1->m_yMax;
    w_z0=model_1->m_zMin;
    w_z1=model_1->m_zMax;

    dx=w_x1-w_x0; dy=w_y1-w_y0; dz=w_z1-w_z0;

    w_x0-=dx; w_x1+=dx;
    w_y0-=dy; w_y1+=dy;
    w_z0-=dz; w_z1+=dz;

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


