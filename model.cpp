#include <model.h>
#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
#include <math.h>
#ifdef __linux__
#include  <GL/gl.h>
#include  <GL/glu.h>
#include  <GL/glut.h>
#elif _WIN32
#include <my_include/gl.h>
#include <my_include/glu.h>
#include <my_include/glut.h>
#endif

using namespace std;
v3::v3(char* facet)
{

    char f1[4] = {facet[0],
        facet[1],facet[2],facet[3]};

    char f2[4] = {facet[4],
        facet[5],facet[6],facet[7]};

    char f3[4] = {facet[8],
        facet[9],facet[10],facet[11]};

    float xx = *((float*) f1 );
    float yy = *((float*) f2 );
    float zz = *((float*) f3 );

    m_x = double(xx);
    m_y = double(yy);
    m_z = double(zz);

}

tri::tri(v3 p1, v3 p2, v3 p3)
{
    m_p[0]=p1;
    m_p[1]=p2;
    m_p[2]=p3;

    //cross product [p1p2 x p1p3]
    normal.m_x = (p2.m_y-p1.m_y)*(p3.m_z-p1.m_z) - (p2.m_z-p1.m_z)*(p3.m_y-p1.m_y);
    normal.m_y = (p2.m_z-p1.m_z)*(p3.m_x-p1.m_x) - (p2.m_x-p1.m_x)*(p3.m_z-p1.m_z);
    normal.m_z = (p2.m_x-p1.m_x)*(p3.m_y-p1.m_y) - (p2.m_x-p1.m_x)*(p3.m_z-p1.m_z);

    double l = sqrt(normal.m_x*normal.m_x + normal.m_y*normal.m_y + normal.m_z*normal.m_z);
    normal.m_x /= l;
    normal.m_y /= l;
    normal.m_z /= l;
}

void tri::draw()
{
 glColor3f(1.0,1.0,1.0);
 glBegin(GL_TRIANGLES);
 for (int i=0;i<3;i++)
    glVertex3f(m_p[i].m_x,m_p[i].m_y,m_p[i].m_z);
 glEnd();

}

void model::load(char *fname)
{

    m_tris.clear();
    ifstream myFile(fname, ios::in | ios::binary);

        char header_info[80] = "";
        char nTri[4];
        unsigned long nTriLong;

        //read 80 byte header
        if (myFile) {
            myFile.read (header_info, 80);
            cout <<"header: " << header_info << endl;
        }
        else{
            cout << "error" << endl;
        }

        //read 4-byte ulong
        if (myFile) {
            myFile.read (nTri, 4);
            nTriLong = *((unsigned long*)nTri) ;
            cout <<"n Tri: " << nTriLong << endl;
        }
        else{
            cout << "error" << endl;
        }

        //now read in all the triangles
        for(int i = 0; i < nTriLong; i++){

            char facet[50];

            if (myFile) {

            //read one 50-byte triangle
                myFile.read (facet, 50);

            //populate each point of the triangle
            //using v3::v3(char* bin);
                //facet + 12 skips the triangle's unit normal
                v3 p1(facet+12);
                v3 p2(facet+24);
                v3 p3(facet+36);

                //add a new triangle to the array
                m_tris.push_back( tri(p1,p2,p3) );

            }
        }

        getMinMax();

        printf("model loaded name=%s tri_count=%d \n xmin=%f xmax=%f \n ymin=%f ymax=%f \n zmin=%f zmax=%f \n",fname,m_tris.size(),
               m_xMin,m_xMax,m_yMin,m_yMax,m_zMin,m_zMax );
        return;

}

void model::draw()
{
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    for (int i=0; i<m_tris.size();i++)
    {
        m_tris[i].draw();
    }

}

void model::getMinMax()
{
    m_xMin = 1e10;    m_xMax = -1e10;
    m_yMin = 1e10;    m_yMax = -1e10;
    m_zMin = 1e10;    m_zMax = -1e10;

    for (int i=0; i<m_tris.size();i++)
    {
        for (int j=0;j<3;j++)
        {
            if (m_tris[i].m_p[j].m_x > m_xMax) m_xMax = m_tris[i].m_p[j].m_x;
            if (m_tris[i].m_p[j].m_x < m_xMin) m_xMin = m_tris[i].m_p[j].m_x;

            if (m_tris[i].m_p[j].m_y > m_yMax) m_yMax = m_tris[i].m_p[j].m_y;
            if (m_tris[i].m_p[j].m_y < m_yMin) m_yMin = m_tris[i].m_p[j].m_y;

            if (m_tris[i].m_p[j].m_z > m_zMax) m_zMax = m_tris[i].m_p[j].m_z;
            if (m_tris[i].m_p[j].m_z < m_zMin) m_zMin = m_tris[i].m_p[j].m_z;
        }

    }


}


