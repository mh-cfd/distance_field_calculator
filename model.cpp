#include <assert.h>
#include <model.h>
#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
#include <math.h>

#include "triangle_box_cross.h"

#ifdef __linux__
#include  <GL/gl.h>
#include  <GL/glu.h>
#include  <GL/glut.h>
#elif _WIN32
#include <windows.h>
#define GLUT_DISABLE_ATEXIT_HACK
#include <my_include/gl.h>
#include <my_include/glu.h>
#include <my_include/glut.h>
#endif

using namespace std;


std::vector<int> model::m_grid[NI_max][NJ_max][NK_max]; //triangles in cells
i3 model::m_nearest_cell[NI_max][NJ_max][NK_max];// nearest non-empty cell to the current cell

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

v3 multiplyMatrPoint3(double **A, v3 b) {
    v3 res;

    res.m_x = A[0][0]*b.m_x + A[0][1]*b.m_y + A[0][2]*b.m_z;
    res.m_y = A[1][0]*b.m_x + A[1][1]*b.m_y + A[1][2]*b.m_z;
    res.m_z = A[2][0]*b.m_x + A[2][1]*b.m_y + A[2][2]*b.m_z;

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

    //find incident angles at vertices;
    v3 a=m_p[1]-m_p[0];
    v3 b=m_p[2]-m_p[0];
    m_vert_angle[0]=acos( dotProd(a,b) / (a.len()*b.len()) );

     a=m_p[0]-m_p[1];
     b=m_p[2]-m_p[1];
    m_vert_angle[1]=acos( dotProd(a,b) / (a.len()*b.len()) );

     a=m_p[1]-m_p[2];
     b=m_p[0]-m_p[2];
    m_vert_angle[2]=acos( dotProd(a,b) / (a.len()*b.len()) );

    if (fabs(m_vert_angle[0] + m_vert_angle[1] +m_vert_angle[2]- 3.1415926535)>1e-5 )
        printf("WARNING triangle angles sum is %f \n",m_vert_angle[0] + m_vert_angle[1] +m_vert_angle[2]);

}

double tri::distP(v3 &point, double &dot_out) {

    double dist = -1;
    v3 p = multiplyMatrPoint(this->matrix, point);
    dot_out=p.m_x;//1.0;
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

   // if(fabs(E12) < 1.e-15 || fabs(E23) < 1.e-15 || fabs(E31) < 1.e-15) {    //on the edge
   //     dist = fabs(p.m_x);
  //      return dist;
  //  }

    if(E12>-1.e-15 && E23>-1.e-15 && E31>-1e-15) {   //inside
        dist = fabs(p.m_x);
        return dist;
    }

    double mult = 0;
    if(E12 < 0) { //left side of p1p2
        mult = dotProd(shP[1], p);
        if(mult < 0) {  //p1 is closest
            dot_out=dot_out*m_vert_angle[0];
            return p.len();
        }
        else {
            mult /= fabs(shP[1].m_z);
            if(mult <= fabs(shP[1].m_z)) {     //p1p2 is closest
                return sqrt(p.m_x*p.m_x + p.m_y*p.m_y);
            }
            else {  //p2 is closest
                dot_out=dot_out*m_vert_angle[1];
                return (p-shP[1]).len();
            }
        }
    }

    if(E31 < 0) {  //right side p1p2, left side p3p1
        mult = dotProd(shP[2], p);
        if(mult < 0) {  //p1 is closest
            dot_out=dot_out*m_vert_angle[0];
            return p.len();
        }
        else {
            mult /= fabs(shP[2].len());
            if(mult <= fabs(shP[2].len())) {     //p1p3 is closest
                double c2 = p.m_y*p.m_y + p.m_z*p.m_z - mult*mult;
                return sqrt(p.m_x*p.m_x + c2);
            }
            else {  //p3 is closest
                dot_out=dot_out*m_vert_angle[2];
                return (p-shP[2]).len();
            }
        }
    }

    if(E23 < 0) {  //left side p2p3, right side of other
        mult = dotProd(shP[2]-shP[1], p-shP[1]);
        if(mult < 0) {  //p2 is closest
            dot_out=dot_out*m_vert_angle[1];
            return (p-shP[1]).len();
        }
        else {
            mult /= (shP[2]-shP[1]).len();
            if(mult <= (shP[2]-shP[1]).len()) {     //p2p3 is closest
                double c2 = pow((p-shP[1]).m_y, 2.) + pow((p-shP[1]).m_z, 2) - mult*mult;
                return sqrt(p.m_x*p.m_x + c2);
            }
            else {  //p3 is closest
                dot_out=dot_out*m_vert_angle[2];
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

double tri::getSign(v3 point)
{
    return 2.0*(dotProd(point - m_p[0],normal)>0)-1.0;
}

void tri::draw()
{
    //glColor3f(1.0,1.0,1.0);
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
            v3 n(facet);
            v3 p1(facet+12);
            v3 p2(facet+24);
            v3 p3(facet+36);

            //add a new triangle to the array
            m_tris.push_back( tri(p1,p2,p3) );

           // m_tris[m_tris.size()-1].normal=n;
        }
    }

    getMinMax();


    printf("model loaded name=%s tri_count=%d \n xmin=%f xmax=%f \n ymin=%f ymax=%f \n zmin=%f zmax=%f \n",fname,m_tris.size(),
           m_xMin,m_xMax,m_yMin,m_yMax,m_zMin,m_zMax );


    double dx,dy,dz;
    w_x0=m_xMin;
    w_x1=m_xMax;
    w_y0=m_yMin;
    w_y1=m_yMax;
    w_z0=m_zMin;
    w_z1=m_zMax;

    dx=0.2*(w_x1-w_x0); dy=0.2*(w_y1-w_y0); dz=0.2*(w_z1-w_z0);

    w_x0-=dx; w_x1+=dx;
    w_y0-=dy; w_y1+=dy;
    w_z0-=dz; w_z1+=dz;

    m_ni=29;
    m_nj=29;
    m_nk=29;

    m_dx=(w_x1-w_x0)/m_ni;
    m_dy=(w_y1-w_y0)/m_nj;
    m_dz=(w_z1-w_z0)/m_nk;
    printf("start distributing \n");

    distributeTriangles();
    printf("end distributing \n");

    return;

}

void model::draw()
{
    glColor3f(0.5,0.5,0.5);
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
    double sum_dot=0.0;//(0,0,0);
double delta=1e-5;
    for (int i=0; i<m_tris.size();i++)
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

    return -min_dist*(2.0*(sum_dot>0.0)-1.0);
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

double model::distP_fast(v3 point)
{
    //printf("aa \n");


    int ix0,iy0,iz0;
    ix0=(int)((point.m_x - w_x0)/m_dx);
    iy0=(int)((point.m_y - w_y0)/m_dy);
    iz0=(int)((point.m_z - w_z0)/m_dz);

    if (ix0<0) ix0=0;
    if (iy0<0) iy0=0;
    if (iz0<0) iz0=0;

    if (ix0>=m_ni-1) ix0=m_ni-1;
    if (iy0>=m_nj-1) iy0=m_nj-1;
    if (iz0>=m_nk-1) iz0=m_nk-1;



    if (m_grid[ix0][iy0][iz0].size()==0)
    {
        i3 nearest=getNearest(ix0,iy0,iz0);
        ix0=nearest.m_i[0];
        iy0=nearest.m_i[1];
        iz0=nearest.m_i[2];
    }
    double min_dist =1e10;

    // printf("bb %d %d %d %d \n",ix,iy,iz,m_grid[ix][iy][iz].size());
    int im,ip,jm,jp,km,kp;
    im=ix0-1; ip=ix0+1;
    jm=iy0-1; jp=iy0+1;
    km=iz0-1; kp=iz0+1;

    if (im<0) im=0;
    if (jm<0) jm=0;
    if (km<0) km=0;

    if (ip>=m_ni-1) ip=m_ni-1;
    if (jp>=m_nj-1) jp=m_nj-1;
    if (kp>=m_nk-1) kp=m_nk-1;

    for (int ix=im;ix<=ip;ix++)
        for (int iy=jm;iy<=jp;iy++)
            for (int iz=km;iz<=kp;iz++)
            {
                for (int i=0; i<m_grid[ix][iy][iz].size();i++)
                {
                    int ind=m_grid[ix][iy][iz][i];
                    double norm;
                    double l=m_tris[ind].distP(point,norm);//*m_tris[ind].getSign(point);
                    if (fabs(min_dist)>fabs(l)) min_dist=l;
                }
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

i3 model::getNearest(int i0, int j0,int k0)
{
    double l2 = 1e20;
    i3 i3m;
    i3m.m_i[0]=i0;
    i3m.m_i[1]=j0;
    i3m.m_i[2]=k0;


    for (int i=0;i<m_ni;i++)
    {
        for (int j=0;j<m_nj;j++)
        {
            for (int k=0;k<m_nk;k++)
            {
                if (m_grid[i][j][k].size()>0)
                {
                    double Dx=(i-i0)*m_dx;
                    double Dy=(j-j0)*m_dy;
                    double Dz=(k-k0)*m_dz;

                    double li=(Dx*Dx + Dy*Dy + Dz*Dz);
                    if (l2>li)
                    {
                        l2=li;
                        i3m.m_i[0]=i;
                        i3m.m_i[1]=j;
                        i3m.m_i[2]=k;
                    }
                }
            }
        }
    }
    return i3m;
}

int model::is_inside(int i, int j, int k)
{

}

int  model::cell_tri_overlap(int i, int j, int k,tri &trian)
{
    float boxcent[3], boxhalf[3],triverts[3][3];
    boxcent[0]=w_x0+m_dx*(i+0.5);
    boxcent[1]=w_y0+m_dy*(j+0.5);
    boxcent[2]=w_z0+m_dz*(k+0.5);

    boxhalf[0]=m_dx*0.5;
    boxhalf[1]=m_dy*0.5;
    boxhalf[2]=m_dz*0.5;
    for (int i=0;i<3;i++)
    {
        triverts[i][0]=trian.m_p[i].m_x;
        triverts[i][1]=trian.m_p[i].m_y;
        triverts[i][2]=trian.m_p[i].m_z;
    }
    return triBoxOverlap(boxcent,boxhalf,triverts);
}

void model::distributeTriangles()
{
    for (int i=0;i<m_ni;i++)
    {
        for (int j=0;j<m_nj;j++)
        {
            for (int k=0;k<m_nk;k++)
            {
                m_grid[i][j][k].clear();
            }
        }
    }


    for (int i=0;i<m_tris.size();i++)
    {
                int ix[3],iy[3],iz[3],ix_min,ix_max,iy_min,iy_max,iz_min,iz_max;
        for (int nn=0; nn<3; nn++)
        {
            ix[nn]=(int)((m_tris[i].m_p[nn].m_x - w_x0)/m_dx);
            iy[nn]=(int)((m_tris[i].m_p[nn].m_y - w_y0)/m_dy);
            iz[nn]=(int)((m_tris[i].m_p[nn].m_z - w_z0)/m_dz);
        }
        FINDMINMAX(ix[0],ix[1],ix[2],ix_min,ix_max);
        FINDMINMAX(iy[0],iy[1],iy[2],iy_min,iy_max);
        FINDMINMAX(iz[0],iz[1],iz[2],iz_min,iz_max);

        for(int ii=ix_min;ii<=ix_max;ii++)
            for(int jj=iy_min;jj<=iy_max;jj++)
                for(int kk=iz_min;kk<=iz_max;kk++)
                {
                    if (cell_tri_overlap(ii,jj,kk,m_tris[i]))
                    {
                         m_grid[ii][jj][kk].push_back(i);
                    }
                }
    }

    /* //from centers
       for (int i=0;i<m_tris.size();i++)
      {
          int ix,iy,iz;
          v3 t0=(m_tris[i].m_p[0]+m_tris[i].m_p[1]+m_tris[i].m_p[2])*(1.0/3.0);
              ix=(int)((t0.m_x - w_x0)/dx);
              iy=(int)((t0.m_y - w_y0)/dy);
              iz=(int)((t0.m_z - w_z0)/dz);

                  m_grid[ix][iy][iz].push_back(i);


      }*/

    //now lets find nearest
    for (int i=0;i<m_ni;i++)
    {
        for (int j=0;j<m_nj;j++)
        {
            for (int k=0;k<m_nk;k++)
            {
                if (m_grid[i][j][k].size()==0)
                {
                    m_nearest_cell[i][j][k]=getNearest(i,j,k);
                }/*else
                {
                    printf("i=%d j=%d k=%d n=%d \n",i,j,k,m_grid[i][j][k].size());
                }*/
            }
        }
    }
}


