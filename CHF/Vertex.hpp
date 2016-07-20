/**
* @file    Vertex.hpp
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
#ifndef _VERTEX_HPP_
#define _VERTEX_HPP_

#include <cmath>
#include <iostream>

/** \brief Vertex integer indice */
typedef int Vid;

/** \brief standart namespace definiton*/
using namespace std;
//--------------------------------------------------//
/** Vertex class for CHF data structure 
  * \brief Vertex class*/
class Vertex
//--------------------------------------------------//
{
//-- Vertex protected data.--//        
protected :
  /** \brief coordenate x of the vertex*/
  float  _x;    
  /** \brief coordenate y of the vertex*/
  float  _y;  
  /** \brief coordenate z of the vertex*/
  float  _z;
  /** \brief coordenate x of the normal vector of the vertex*/
  float _nx;
   /**\brief coordenate y x of the normal vector of the vertex*/
  float _ny;
  /** \brief coordenate z x of the normal vector of the vertex*/
  float _nz;
  /** \brief scalar field of the vertex*/
  float _f;

//-- Vertex constructors.--//        
public:
  /** \brief default constructor */
	Vertex():_x(FLT_MAX), _y(FLT_MAX), _z(FLT_MAX), _nx(FLT_MAX), _ny(FLT_MAX), _nz(FLT_MAX), _f(FLT_MAX) {}
	
  /**
    * \brief copy constructor
    * \param v - const Vertex object*/
    Vertex( const Vertex& v ): _x(v.x()), _y(v.y()), _z(v.z()), _nx(v.nx()), _ny(v.ny()), _nz(v.nz()), _f(v.f()) {} 
	
	/** \brief destructor */
	~Vertex(){}
//--------------------------------------------------//
// Accessors
//--------------------------------------------------//
public:
  /** \brief access to the coordinate x of the vertex*/
  inline const float     x() const{ return _x ; }

  /** \brief access to the coordinate y of the vertex*/
  inline const float     y() const{ return _y ; }

  /** \brief access to the coordinate z of the vertex*/
  inline const float     z() const{ return _z ; }
	
  /** \brief access to the coordenate x of the normal vector of the vertex*/
  inline const float    nx() const{ return _nx; } 

  /** \brief access to the coordenate y of the normal vector of the vertex*/
  inline const float    ny() const{ return _ny; }

  /** \brief access to the coordenate z of the normal vector of the vertex*/
  inline const float    nz() const{ return _nz; }

  /** \brief access to the scalar field of the vertex*/
  inline const float     f() const{ return  _f; }


  /** \brief sets the coordenate x of the vertex
    * \param x - const float */
  inline const void  set_x (const float x ) { _x=x;   }
 
  /** \brief sets the coordenate y of the vertex
    * \param y - const float */
  inline const void  set_y (const float y ) { _y=y;   }
	
  /** \brief sets the coordenate z of the vertex
    * \param z - const float */
  inline const void  set_z (const float z ) { _z=z;   }
	
  /** \brief sets the coordenate x of the normal vector of the vertex
    * \param nx - const float */
  inline const void  set_nx(const float nx) { _nx=nx; }

  /** \brief sets the coordenate y of the normal vector of the vertex
    * \param ny - const float */
  inline const void  set_ny(const float ny) { _ny=ny; }

  /** \brief sets the coordenate z of the normal vector of the vertex
    * \param nz - const float */
  inline const void  set_nz(const float nz) { _nz=nz; }

  /** \brief sets the field the vertex
    * \param f - const float */
  inline const void  set_f(const float f) { _f=f; }

/** Static Methods*/
public:
  /** \brief computes the squared norm of a vector in R^4
    * \param v[4] - const float */
  inline static float norm2_4( const float v[4] ) { return (v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3]) ; }

  /** z\brief computes the squared norm of a vector in R^4
    * \param a - const float 
    * \param b - const float 
    * \param c - const float 
    * \param d - const float */
  inline static float norm2_4( const float a, const float b, const float c, const float d ) { return a*a + b*b + c*c + d*d ; }

  /** \brief computes the norm of a vector in R^4
    * \param v[4] - const float */
  inline static float  norm_4( const float v[4] ) { return (float) sqrt( norm2_4(v) ) ; }

  /** \brief computes the norm of a vector in R^4
    * \param a - const float 
    * \param b - const float 
    * \param c - const float 
    * \param d - const float */
  inline static float  norm_4( const float a, const float b, const float c, const float d ) { return (float) sqrt( norm2_4(a,b,c,d) ) ; }

  /** \brief normalize a vector in R^4
    * \param v - const float* */
  static void normalize_4( float *v );

  /** \brief computes the squared norm of a vector in R^3
    * \param v[3] - const float */
  inline static float norm2( const float v[3] ) { return v[0]*v[0] + v[1]*v[1] + v[2]*v[2] ; }
	
  /** \brief computes the squared norm of a vector in R^3
    * \param a - const float 
    * \param b - const float 
    * \param c - const float */
  inline static float norm2( const float a, const float b, const float c ) { return a*a + b*b + c*c ; }

  /** \brief computes the norm of a vector in R^3
    * \param v[3] - const float */
  inline static float norm( const float v[3] ) { return (float) sqrt( norm2(v) ) ; }

  /** \brief computes the norm of a vector in R^3
    * \param a - const float 
    * \param b - const float 
    * \param c - const float */ 
  inline static float norm( const float a, const float b, const float c ) { return (float) sqrt( norm2(a,b,c) ) ; }

  /** \brief normalize a vector in R^3
    * \param - v const float* */
  static void normalize( float *v );

  /** c\brief computes the normal of a triangle
    * \param v0 - const Vertex &
    * \param v1 - const Vertex &
    * \param v2 - const Vertex &
    * \param n  - float* */
  static void normal( const Vertex &v0, const Vertex &v1, const Vertex &v2, float *n);

  /** \brief computes the aspect ratio of a tetrahedron
    * \param v0 - const Vertex&
    * \param v1 - const Vertex&
    * \param v2 - const Vertex&
    * \param v3 - const Vertex& */
  static const float tetra_aspect_ratio( const Vertex &v0, const Vertex &v1, const Vertex &v2, const Vertex &v3 );
	
   /** \brief computes the volume of a tetrahedron
     * \param v0 - const Vertex&
     * \param v1 - const Vertex&
     * \param v2 - const Vertex&
     * \param v3 - const Vertex& */
  static const float tetra_volume( const Vertex &v0, const Vertex &v1, const Vertex &v2, const Vertex &v3 );

  /** \brief computes the signed volume of a tetrahedron
    * \param v0 - const Vertex&
    * \param v1 - const Vertex&
    * \param v2 - const Vertex&
    * \param v3 - const Vertex& */
  static const float signed_tetra_volume( const Vertex &v0, const Vertex &v1, const Vertex &v2, const Vertex &v3 );
  
  /** \brief computes the basis the tangent space of a tetrahedron
    * \param v0 - const Vertex&
    * \param v1 - const Vertex&
    * \param v2 - const Vertex&
    * \param base - const float**/
  static void tetra_base( const Vertex &v0, const Vertex &v1, const Vertex &v2, const Vertex &v3,  float *base);
  
  /** \brief computes the area of a tetrahedron
    * \param v0 - const Vertex&
    * \param v1 - const Vertex&
    * \param v2 - const Vertex&*/
  static const float  trig_area( const Vertex &v0, const Vertex &v1, const Vertex &v2 );
    
  /** Computes a basis for the tangent space of a triangle
    * \brief computes a basis for the tangent space of a triangle
    * \param v0 - const Vertex&
    * \param v1 - const Vertex&
    * \param v2 - const Vertex&
    * \param base - const float* */
  static void  trig_base( const Vertex &v0, const Vertex &v1, const Vertex &v2, float *base);

//-- Vertex bool operators --//
public: 
  /** \brief Equal operator for Vertex objects 
    * \param Vertex& v */
   bool operator == (const Vertex& v){ return (x()==v.x() && y()==v.y() && z()==v.z() && nx()==v.nx() && ny()==v.ny() && nz()==v.nz() && f()==v.f());}
  
  /** \brief Equal operator for Vertex objects 
    * \param Vertex& v */
   bool operator != (const Vertex& v){ return (x()!=v.x() || y()!=v.y() || z()!=v.z() || nx()!=v.nx() || ny()!=v.ny() || nz()!=v.nz() || f()!=v.f());}
};
#endif
//--------------------------------------------------//
