#include <assert.h>
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

void multiplyMatrVect(double **A, v3 *b) {
    v3 temp[3];

    for(int i=0;i<3;i++) {
        temp[i].m_x = A[0][0]*b[i].m_x + A[0][1]*b[i].m_y + A[0][2]*b[i].m_z + A[0][3];
        temp[i].m_y = A[1][0]*b[i].m_x + A[1][1]*b[i].m_y + A[1][2]*b[i].m_z + A[1][3];
        temp[i].m_z = A[2][0]*b[i].m_x + A[2][1]*b[i].m_y + A[2][2]*b[i].m_z + A[2][3];
    }

    for(int i=0;i<3;i++) {
        b[i] = temp[i];
    }
}

v3 multiplyMatrPoint(double **A, v3 b) {
    v3 res;

    res.m_x = A[0][0]*b.m_x + A[0][1]*b.m_y + A[0][2]*b.m_z + A[0][3];
    res.m_y = A[1][0]*b.m_x + A[1][1]*b.m_y + A[1][2]*b.m_z + A[1][3];
    res.m_z = A[2][0]*b.m_x + A[2][1]*b.m_y + A[2][2]*b.m_z + A[2][3];

    return res;
}

void multiplyMatrixes(double **A, double **B, int N) {
    double temp[N][N];

    for(int i=0;i<N;i++) {
        for(int j=0;j<N;j++) {
            temp[i][j] = 0;
            for(int m=0;m<N;m++) {
                    temp[i][j] += A[i][m]*B[m][j];
            }
        }
    }

    for(int i=0;i<N;i++) {
        for(int j=0;j<N;j++) {
            B[i][j] = temp[i][j];
        }
    }
}

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

double v3::len()
{
    return sqrt(m_x*m_x+m_y*m_y+ m_z*m_z);
}

double dotProd(v3 vec1, v3 vec2) {
    return vec1.m_x*vec2.m_x + vec1.m_y*vec2.m_y + vec1.m_z*vec2.m_z;
}

v3 operator+(v3 left, v3 right){
    v3 res;

    res.m_x = left.m_x + right.m_x;
    res.m_y = left.m_y + right.m_y;
    res.m_z = left.m_z + right.m_z;

    return res;
}

v3 operator-(v3 left, v3 right){
    v3 res;

    res.m_x = left.m_x - right.m_x;
    res.m_y = left.m_y - right.m_y;
    res.m_z = left.m_z - right.m_z;

    return res;
}

v3 operator/(v3 left, double right){
    v3 res;

    res.m_x = left.m_x / right;
    res.m_y = left.m_y / right;
    res.m_z = left.m_z / right;

    return res;
}


v3 operator*(v3 left, double right){
    v3 res;

    res.m_x = left.m_x * right;
    res.m_y = left.m_y * right;
    res.m_z = left.m_z * right;

    return res;
}


tri::tri(v3 p1, v3 p2, v3 p3)
{
    m_p[0]=p1;
    m_p[1]=p2;
    m_p[2]=p3;

    //cross product [p1p2 x p1p3]
    normal.m_x = (p2.m_y-p1.m_y)*(p3.m_z-p1.m_z) - (p2.m_z-p1.m_z)*(p3.m_y-p1.m_y);
    normal.m_y = (p2.m_z-p1.m_z)*(p3.m_x-p1.m_x) - (p2.m_x-p1.m_x)*(p3.m_z-p1.m_z);
    normal.m_z = (p2.m_x-p1.m_x)*(p3.m_y-p1.m_y) - (p2.m_y-p1.m_y)*(p3.m_x-p1.m_x);

    double l = sqrt(normal.m_x*normal.m_x + normal.m_y*normal.m_y + normal.m_z*normal.m_z);
    normal.m_x /= l;
    normal.m_y /= l;
    normal.m_z /= l;

    center = (p1+p2+p3)/3;

    shP = new v3[3]; //shifted points

    matrix = new double*[4];
    for(int i=0;i<4;i++) {
        matrix[i] = new double[4];
    }

    for(int i=0;i<4;i++) {
        for(int j=0;j<4;j++) {
            if(i!=j) matrix[i][j] = 0.;
            else matrix[i][j] = 1.;
        }
    }

    //translation for p1 is in the origin
    matrix[0][3] = -p1.m_x;
    matrix[1][3] = -p1.m_y;
    matrix[2][3] = -p1.m_z;

    for(int i=0;i<3;i++) shP[i] = m_p[i]-p1;

    //cout << "shP1=" << shP[0].m_x << "   " << shP[0].m_y << "   " << shP[0].m_z << endl;
    //cout << "shP2=" << shP[1].m_x << "   " << shP[1].m_y << "   " << shP[1].m_z << endl;

    double **rot_matrix;   //rotation matrix

    rot_matrix = new double*[4];
    for(int i=0;i<4;i++) {
        rot_matrix[i] = new double[4];
    }

    for(int i=0;i<4;i++) {
        for(int j=0;j<4;j++) {
            if(i!=j) rot_matrix[i][j] = 0.;
            else rot_matrix[i][j] = 1.;
        }
    }

    //angle between vector p1p2 and plane y=0
    double alpha = acos(shP[1].m_y/sqrt(shP[1].m_x*shP[1].m_x + shP[1].m_y*shP[1].m_y));
    if(shP[1].m_x > 0) alpha = -alpha;
    //cout << "alpha=" << alpha*180/M_PI << endl;
    //rotation for p1p2 is on the xz-plane
    rot_matrix[0][0] = cos(-alpha + M_PI/2.);
    rot_matrix[0][1] = -sin(-alpha + M_PI/2.);
    rot_matrix[1][0] = sin(-alpha + M_PI/2.);
    rot_matrix[1][1] = cos(-alpha + M_PI/2.);

    multiplyMatrVect(rot_matrix, shP);  //should be p2.y=0

    //cout << "shP2.y=" << shP[1].m_y << endl;

    multiplyMatrixes(rot_matrix, matrix, 4);

    for(int i=0;i<4;i++) {
        for(int j=0;j<4;j++) {
            if(i!=j) rot_matrix[i][j] = 0.;
            else rot_matrix[i][j] = 1.;
        }
    }

    //angle between vector p1p2 and plane x=0
    double beta = acos(shP[1].m_x/sqrt(shP[1].m_x*shP[1].m_x + shP[1].m_z*shP[1].m_z));
    if(shP[1].m_z < 0) beta = -beta;
    //cout << "beta=" << beta*180/M_PI << endl;
    //rotation for p1p2 and z-axis are co-directed
    rot_matrix[0][0] = cos(beta - M_PI/2.);
    rot_matrix[0][2] = sin(beta - M_PI/2.);
    rot_matrix[2][0] = -sin(beta - M_PI/2.);
    rot_matrix[2][2] = cos(beta - M_PI/2.);

    //cout << "shP2=" << shP[1].m_x << "   " << shP[1].m_y << "   " << shP[1].m_z << endl;

    multiplyMatrVect(rot_matrix, shP);  //should be p2.x=0

    //cout << "shP2.x=" << shP[1].m_x << endl;
    //cout << "shP2.y=" << shP[1].m_y << endl;

    multiplyMatrixes(rot_matrix, matrix, 4);

    for(int i=0;i<4;i++) {
        for(int j=0;j<4;j++) {
            if(i!=j) rot_matrix[i][j] = 0.;
            else rot_matrix[i][j] = 1.;
        }
    }

    //cout << "shP1=" << shP[0].m_x << "   " << shP[0].m_y << "   " << shP[0].m_z << endl;
    //cout << "shP2=" << shP[1].m_x << "   " << shP[1].m_y << "   " << shP[1].m_z << endl;
    //cout << "shP3=" << shP[2].m_x << "   " << shP[2].m_y << "   " << shP[2].m_z << endl;

    //angle between vector p1p3 and plane x=0
    double gamma = acos(shP[2].m_x/sqrt(shP[2].m_x*shP[2].m_x + shP[2].m_y*shP[2].m_y));
    if(shP[2].m_y < 0) gamma = -gamma;
    //rotation for p3p1 is on the yz-plane
    rot_matrix[0][0] = cos(-gamma + M_PI/2.);
    rot_matrix[0][1] = -sin(-gamma + M_PI/2.);
    rot_matrix[1][0] = sin(-gamma + M_PI/2.);
    rot_matrix[1][1] = cos(-gamma + M_PI/2.);

    multiplyMatrVect(rot_matrix, shP); //should be p3.x=0

    multiplyMatrixes(rot_matrix, matrix, 4);

    //cout << "shP3.x=" << shP[2].m_x << endl;

    //cout << "shP1=" << shP[0].m_x << "   " << shP[0].m_y << "   " << shP[0].m_z << endl;
    //cout << "shP2=" << shP[1].m_x << "   " << shP[1].m_y << "   " << shP[1].m_z << endl;
    //cout << "shP3=" << shP[2].m_x << "   " << shP[2].m_y << "   " << shP[2].m_z << endl;

    //for(int i=0;i<3;i++) shP[i] = m_p[i];

    //multiplyMatrVect(matrix, shP); //exam

    //cout << "shP1=" << shP[0].m_x << "   " << shP[0].m_y << "   " << shP[0].m_z << endl;
    //cout << "shP2=" << shP[1].m_x << "   " << shP[1].m_y << "   " << shP[1].m_z << endl;
    //cout << "shP3=" << shP[2].m_x << "   " << shP[2].m_y << "   " << shP[2].m_z << endl;

    //cout << endl;

    //delete[] shP;

    for(int i=0;i<4;i++) delete[] rot_matrix[i];
    delete[] rot_matrix;
}

double tri::distP(v3 point) {
    double dist = -1;
    v3 p = multiplyMatrPoint(this->matrix, point);

    double E12, E23, E31; //edges equation results: E<0 -> point on left side, E>0 -> point on right side
    //E(y,z) = (y-Y)*dZ - (z-Z)*dY

    E12 = p.m_y*shP[1].m_z;
    E23 = p.m_y*(shP[2].m_z-shP[1].m_z) - (p.m_z-shP[1].m_z)*shP[2].m_y;
    E31 = p.m_y*(-shP[2].m_z) - p.m_z*(-shP[2].m_y);

    //cout << "shP1=" << shP[0].m_x << "   " << shP[0].m_y << "   " << shP[0].m_z << endl;
    //cout << "shP2=" << shP[1].m_x << "   " << shP[1].m_y << "   " << shP[1].m_z << endl;
    //cout << "shP3=" << shP[2].m_x << "   " << shP[2].m_y << "   " << shP[2].m_z << endl;
    //cout << "p=" << p.m_x << "   " << p.m_y << "   " << p.m_z << endl;
    //cout << "Edges " << E12 << "   " << E23 << "   " << E31 << endl;

    if(fabs(E12) < 1.e-15 || fabs(E23) < 1.e-15 || fabs(E31) < 1.e-15) {    //on the edge
        dist = fabs(p.m_x);
        return dist;
    }

    if(E12>0 && E23>0 && E31>0) {   //inside
        dist = fabs(p.m_x);
        return dist;
    }

    double mult = 0;
    if(E12 < 0) { //left side of p1p2
        mult = dotProd(shP[1], p);
        if(mult < 0) {  //p1 is closest
            return p.len();
        }
        else {
            mult /= fabs(shP[1].m_z);
            if(mult <= fabs(shP[1].m_z)) {     //p1p2 is closest
                return sqrt(p.m_x*p.m_x + p.m_y*p.m_y);
            }
            else {  //p2 is closest
                return (p-shP[1]).len();
            }
        }
    }

    if(E31 < 0) {  //right side p1p2, left side p3p1
        mult = dotProd(shP[2], p);
        if(mult < 0) {  //p1 is closest
            return p.len();
        }
        else {
            mult /= fabs(shP[2].len());
            if(mult <= fabs(shP[2].len())) {     //p1p3 is closest
                double c2 = p.m_y*p.m_y + p.m_z*p.m_z - mult*mult;
                return sqrt(p.m_x*p.m_x + c2);
            }
            else {  //p3 is closest
                return (p-shP[2]).len();
            }
        }
    }

    if(E23 < 0) {  //left side p2p3, right side of other
        mult = dotProd(shP[2]-shP[1], p-shP[1]);
        if(mult < 0) {  //p2 is closest
            return (p-shP[1]).len();
        }
        else {
            mult /= (shP[2]-shP[1]).len();
            if(mult <= (shP[2]-shP[1]).len()) {     //p2p3 is closest
                double c2 = pow((p-shP[1]).m_y, 2.) + pow((p-shP[1]).m_z, 2) - mult*mult;
                return sqrt(p.m_x*p.m_x + c2);
            }
            else {  //p3 is closest
                return (p-shP[2]).len();
            }
        }
    }

    assert(dist < 0);

    return dist;
    //return (point - this->center).len();
}

double tri::distP_naive(v3 point)
{
    v3 delta=(m_p[0]+m_p[1]+m_p[2])*(1.0/3.0) - point;
    return delta.len();
}

void tri::draw()
{
     glColor3f(1.0,1.0,1.0);
     glBegin(GL_TRIANGLES);
            glNormal3f(normal.m_x,normal.m_y,normal.m_z);
        for (int i=0;i<3;i++)
            glVertex3f(m_p[i].m_x,m_p[i].m_y,m_p[i].m_z);
     glEnd();

     glBegin(GL_LINES);
        glVertex3f(center.m_x, center.m_y, center.m_z);
        glVertex3f(center.m_x+normal.m_x, center.m_y+normal.m_y, center.m_z+normal.m_z);
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

    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
}

double model::distP(v3 point)
{
    double min_dist =1e10;
    for (int i=0; i<m_tris.size();i++)
    {
        double l=m_tris[i].distP(point);
        if (fabs(min_dist)>fabs(l)) min_dist=l;
    }
    return min_dist;
}

double model::distP_naive(v3 point)
{
    double min_dist =1e10;
    for (int i=0; i<m_tris.size();i++)
    {
        double l=m_tris[i].distP_naive(point);
        if (min_dist>l) min_dist=l;
    }
    return min_dist;
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


