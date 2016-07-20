/**
* @file    Vertex.cpp
* @author  Marcos Lage      <mlage@mat.puc-rio.br>
* @author  Thomas Lewiner   <thomas.lewiner@polytechnique.org>
* @author  Hélio  Lopes     <lopes@mat.puc-rio.br>
* @author  Math Dept, PUC-Rio
* @author  Lab Matmídia
* @date    14/02/2006
*
* @brief   CHF: A Scalable topological data structure for tetrahedral meshes
* @brief  (Vertex Structure)
*/
//--------------------------------------------------//
#include <cfloat>
#include <vector>
#include <iostream>
#include <algorithm>

#include "Vertex.hpp"

using namespace std;
//--------------------------------------------------//
void Vertex::normalize_4( float *v )
//--------------------------------------------------//
/** Normalizes a vector in R^4*/
{
  float n= norm (v);

  if(n > FLT_EPSILON )
  {
    v[0]/=n;
    v[1]/=n;
    v[2]/=n;
    v[3]/=n;
  }
}
//--------------------------------------------------//
void Vertex::normalize( float *v )
//--------------------------------------------------//
/** Normalizes a vector in R^3*/
{
  float n= norm(v);

  if(n > FLT_EPSILON )
  {
    v[0]/=n;
    v[1]/=n;
    v[2]/=n;
  }
}
//--------------------------------------------------//
void Vertex::normal (const Vertex &v0, const Vertex &v1, const Vertex &v2, float *norm)
//--------------------------------------------------//
/** Computes the normal of a vector in R^3*/
{
  float a0[3], a1[3];

  norm[0]= 0;
  norm[1]= 0;
  norm[2]= 0;

  a0[0]= v1.x()-v0.x();
  a0[1]= v1.y()-v0.y();
  a0[2]= v1.z()-v0.z();

  a1[0]= v2.x()-v0.x();
  a1[1]= v2.y()-v0.y();
  a1[2]= v2.z()-v0.z();

  norm[0]= a0[1]*a1[2]-a1[1]*a0[2];
  norm[1]= a0[2]*a1[0]-a1[2]*a0[0];
  norm[2]= a0[0]*a1[1]-a1[0]*a0[1];

  normalize(norm);
}
//--------------------------------------------------//
const float Vertex::tetra_aspect_ratio( const Vertex &v0, const Vertex &v1, const Vertex &v2, const Vertex &v3 )
//--------------------------------------------------//
/** Computes a tetrahedron aspect ratio*/
{
  Vertex center;
  vector<float> a;
  float m, radio;

  a.resize(6, -1);
  
  a[0] = sqrt( (v1.x() - v0.x())*(v1.x() - v0.x()) + (v1.y() - v0.y())*(v1.y() - v0.y()) + (v1.z() - v0.z())*(v1.z() - v0.z()) );
  a[1] = sqrt( (v2.x() - v0.x())*(v2.x() - v0.x()) + (v2.y() - v0.y())*(v2.y() - v0.y()) + (v2.z() - v0.z())*(v2.z() - v0.z()) );
  a[2] = sqrt( (v3.x() - v0.x())*(v3.x() - v0.x()) + (v3.y() - v0.y())*(v3.y() - v0.y()) + (v3.z() - v0.z())*(v3.z() - v0.z()) );
 
  a[3] = sqrt( (v1.x() - v2.x())*(v1.x() - v2.x()) + (v1.y() - v2.y())*(v1.y() - v2.y()) + (v1.z() - v2.z())*(v1.z() - v2.z()) );
  a[4] = sqrt( (v2.x() - v3.x())*(v2.x() - v3.x()) + (v2.y() - v3.y())*(v2.y() - v3.y()) + (v2.z() - v3.z())*(v2.z() - v3.z()) );
  a[5] = sqrt( (v3.x() - v1.x())*(v3.x() - v1.x()) + (v3.y() - v1.y())*(v3.y() - v1.y()) + (v3.z() - v1.z())*(v3.z() - v1.z()) );
  
  m = *max_element(a.begin(), a.end());

  center.set_x( (v0.x() + v1.x() + v2.x() + v3.x())/4 );
  center.set_y( (v0.y() + v1.y() + v2.y() + v3.y())/4 );
  center.set_z( (v0.z() + v1.z() + v2.z() + v3.z())/4 );

  radio = sqrt( (center.x() - v1.x())*(center.x() - v1.x()) + (center.y() - v1.y())*(center.y() - v1.y()) + (center.z() - v1.z())*(center.z() - v1.z()) );

  return (float)( m/(2*sqrt(6.0)*radio) );
}
//--------------------------------------------------//
const float Vertex::tetra_volume(const Vertex &v0, const Vertex &v1, const Vertex &v2, const Vertex &v3)
//--------------------------------------------------//
/** Computes a tetrahedron volume*/
{
  float a[4], b[4], c[4], prod[4];
  float volume;

  a[0]=v1.x()-v0.x();
  b[0]=v2.x()-v0.x(); 
  c[0]=v3.x()-v0.x();

  a[1]=v1.y()-v0.y();
  b[1]=v2.y()-v0.y();
  c[1]=v3.y()-v0.y();

  a[2]=v1.z()-v0.z();
  b[2]=v2.z()-v0.z();
  c[2]=v3.z()-v0.z();

  a[3]=v1.f()-v0.f();
  b[3]=v2.f()-v0.f();
  c[3]=v3.f()-v0.f();

  prod[0]=  (a[1]*b[2]*c[3] + a[2]*b[3]*c[1] + a[3]*b[1]*c[2] - a[3]*b[2]*c[1] - a[2]*b[1]*c[3] - a[1]*b[3]*c[2]);
  prod[1]= -(a[0]*b[2]*c[3] + a[2]*b[3]*c[0] + a[3]*b[0]*c[2] - a[3]*b[2]*c[0] - a[2]*b[0]*c[3] - a[0]*b[3]*c[2]);
  prod[2]=  (a[0]*b[1]*c[3] + a[1]*b[3]*c[0] + a[3]*b[0]*c[1] - a[3]*b[1]*c[0] - a[1]*b[0]*c[3] - a[0]*b[3]*c[1]);
  prod[3]= -(a[0]*b[1]*c[2] + a[1]*b[2]*c[0] + a[2]*b[0]*c[1] - a[2]*b[1]*c[0] - a[1]*b[0]*c[2] - a[0]*b[2]*c[1]);
   
 volume= (float)((0.166667)*norm_4(prod));

 return volume;
}
//--------------------------------------------------//
const float Vertex::signed_tetra_volume(const Vertex &v0, const Vertex &v1, const Vertex &v2, const Vertex &v3)
//--------------------------------------------------//
/** Computes a tetrahedron signed volume*/
{
  float a[4], b[4], c[4];
  float volume;

  a[0]=v1.x()-v0.x();
  b[0]=v2.x()-v0.x();
  c[0]=v3.x()-v0.x();

  a[1]=v1.y()-v0.y();
  b[1]=v2.y()-v0.y();
  c[1]=v3.y()-v0.y();

  a[2]=v1.z()-v0.z();
  b[2]=v2.z()-v0.z();
  c[2]=v3.z()-v0.z();

  volume= a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0]) + a[2]*(b[0]*c[1]-b[1]*c[0]);

  return volume;
}
//--------------------------------------------------//
void Vertex::tetra_base(const Vertex &v0, const Vertex &v1, const Vertex &v2, const Vertex &v3, float *base )
//--------------------------------------------------//
/** Computes a basis for the tangent plane of a tetrahedron*/
{
  float a_0[4],a_1[4], a_2[4], temp[4];

  a_0[0]=v1.x()-v0.x();
  a_0[1]=v1.y()-v0.y();
  a_0[2]=v1.z()-v0.z();
  a_0[3]=v1.f() -v0.f();

  a_1[0]=v2.x()-v0.x();
  a_1[1]=v2.y()-v0.y();
  a_1[2]=v2.z()-v0.z();
  a_1[3]=v2.f() -v0.f();

  a_2[0]=v3.x()-v0.x();
  a_2[1]=v3.y()-v0.y();
  a_2[2]=v3.z()-v0.z();
  a_2[3]=v3.f() -v0.f();


  for(int i=0; i<4; i++)
    temp[i]= a_1[i]- (	(a_1[0]*a_0[0] + a_1[1]*a_0[1] + a_1[2]*a_0[2] + a_1[3]*a_0[3])/norm2_4(a_0)	)*a_0[i];

  for(int i=0; i<4; i++)
    a_1[i]=temp[i];

  for(int i=0; i<4; i++)
    temp[i]= a_2[i] - ( ( (a_2[0]*a_0[0] + a_2[1]*a_0[1] + a_2[2]*a_0[2] + a_2[3]*a_0[3] )/norm2_4(a_0) )*a_0[i]
    +( (a_2[0]*a_1[0] + a_2[1]*a_1[1] + a_2[2]*a_1[2] + a_2[3]*a_1[3] )/norm2_4(a_1) )*a_1[i] ) ;

    for(int i=0; i<4; i++)
      a_2[i]=temp[i];

    normalize_4( a_0 );
    normalize_4( a_1 );
    normalize_4( a_2 );

    base[0] =a_0[0];
    base[1] =a_0[1];
    base[2] =a_0[2];
    base[3] =a_0[3];

    base[4] =a_1[0];
    base[5] =a_1[1];
    base[6] =a_1[2];
    base[7] =a_1[3];

    base[8]  =a_2[0];
    base[9]  =a_2[1];
    base[10]=a_2[2];
    base[11]=a_2[3];
}
//--------------------------------------------------//
const float Vertex::trig_area(const Vertex &v0, const Vertex &v1, const Vertex &v2)
//--------------------------------------------------//
/** Computes the area of a triangle*/
{
  float a[3], b[3], prod[3];
  float area;

  a[0]=v1.x()-v0.x();
  b[0]=v2.x()-v0.x();

  a[1]=v1.y()-v0.y();
  b[1]=v2.y()-v0.y();

  a[2]=v1.z()-v0.z();
  b[2]=v2.z()-v0.z();

  prod[0]=    a[0]*b[1]- a[1]*b[0];
  prod[1]= -( a[0]*b[2]- a[2]*b[0] );
  prod[2]=    a[1]*b[2]- a[2]*b[1];

  area = norm( prod[0], prod[1], prod[2] )/2;

  return area;
}
//--------------------------------------------------//
void Vertex::trig_base(const Vertex &v0, const Vertex &v1, const Vertex &v2, float *base )
//--------------------------------------------------//
/** Computes a basis for the tangent plane of a triangle*/
{
  float a_0[4],a_1[4], temp[4];

  a_0[0]=v1.x()-v0.x();
  a_0[1]=v1.y()-v0.y();
  a_0[2]=v1.z()-v0.z();
  a_0[3]=v1.f()-v0.f();

  a_1[0]=v2.x()-v0.x();
  a_1[1]=v2.y()-v0.y();
  a_1[2]=v2.z()-v0.z();
  a_1[3]=v2.f()-v0.f();

  for(int i=0; i<4; i++)
    temp[i]= a_1[i]- (	(a_1[0]*a_0[0] + a_1[1]*a_0[1] + a_1[2]*a_0[2])/norm2_4(a_0)	)*a_0[i];

  for(int i=0; i<4; i++)
    a_1[i]=temp[i];

  normalize_4( a_0 );
  normalize_4( a_1 );

    base[0] =a_0[0];
    base[1] =a_0[1];
    base[2] =a_0[2];
    base[3] =a_0[3];

    base[4] =a_1[0];
    base[5] =a_1[1];
    base[6] =a_1[2];
    base[7] =a_1[3];
}
//--------------------------------------------------//
