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

//draw call
void display(void)
{
        glClear(GL_COLOR_BUFFER_BIT);
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
    if (key=='.')
    {

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
    w_x0=model_1->m_xMin;
    w_x1=model_1->m_xMax;
    w_y0=model_1->m_yMin;
    w_y1=model_1->m_yMax;
    w_z0=model_1->m_zMin;
    w_z1=model_1->m_zMax;


    glutInit(&argc,argv);
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
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


