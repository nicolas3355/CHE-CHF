/**
* @file    Vertex.hpp
* @author  Marcos Lage      <mlage@mat.puc-rio.br>
* @author  Thomas Lewiner   <thomas.lewiner@polytechnique.org>
* @author  HŽlio  Lopes     <lopes@mat.puc-rio.br>
* @author  Math Dept, PUC-Rio
* @author  Lab Matmídia
* @date    14/02/2006
*
* @brief  (Vertex Class).
*/
//--------------------------------------------------//

#ifndef _VERTEX_HPP_
#define _VERTEX_HPP_

#include <math.h>
#include <float.h>
#include <iostream>

/** Vertex id type; */
typedef int Vid;

using namespace std;
//--------------------------------------------------//
/** Vertex class for CHE data structure 
  * \brief Vertex class*/
class Vertex
//--------------------------------------------------//
{
//-- Vertex protected data --//
protected:
  /** Vertex coordenates
    * \brief coordenate x of the vertex*/
	double  _x;    
  /** Vertex coordenates
    * \brief coordenate y of the vertex*/
  double  _y;  
  /** Vertex coordenates
    * \brief coordenate z of the vertex*/
  double  _z;
  /** Vertex normal coordenates
    * \brief coordenate x of the vertex normal*/
  double _nx;
   /** Vertex normal coordenates
    * \brief coordenate y of the vertex normal*/
  double _ny;
  /** Vertex normal coordenates
    * \brief coordenate z of the vertex normal*/
	double _nz;
  /** Vertex scalar field
    * \brief scalar field of the vertex*/
	double _field;

//-- Vertex constructors.--// 
public:
  /** \brief Default constructor.*/
	Vertex():_x(0), _y(0), _z(0), _nx(0), _ny(0), _nz(0), _field(0) {}

  /** \brief First constructor
    * \param double x  Vertex _x
    * \param double y  Vertex _y
    * \param double z  Vertex _z */
	Vertex( double  x, double  y, double  z ):_x(x), _y(y), _z(z), _nx(0), _ny(0), _nz(0), _field(0) {}

  /** \brief Second constructor
    * \param double x  Vertex _x
    * \param double y  Vertex _y
    * \param double z  Vertex _z 
    * \param double nx Vertex _nx 
    * \param double ny Vertex _ny 
    * \param double nz Vertex _nz 
    * \param double f  Vertex _field */  
	Vertex( double  x, double  y, double  z, double nx, double ny, double nz, double f ):
              _x(x), _y(y), _z(z), _nx(nx), _ny(ny), _nz(nz), _field(f) {}

  /** \brief Copy constructor
    * \param Vertex& v.*/
  Vertex(const Vertex& v):
             _x(v.x()), _y(v.y()), _z(v.z()), _nx(v.nx()), _ny(v.ny()), _nz(v.nz()), _field(v.field()) {}

  /** \brief Destructor.*/
  ~Vertex(){}

//-- Vertex methods.--//
public:
	/** \brief Access to the coordinate x of the vertex*/
  inline const double     x() const{ return _x ; }

	/** \brief Access to the coordinate y of the vertex*/
  inline const double     y() const{ return _y ; }

	/** \brief Access to the coordinate z of the vertex*/
  inline const double     z() const{ return _z ; }
	
	/** \brief Access to the coordenate nx of the normal*/
  inline const double    nx() const{ return _nx; } 

	/** \brief Access to the coordenate ny of the normal*/
  inline const double    ny() const{ return _ny; }

	/** \brief Access to the coordenate nz of the normal*/
  inline const double    nz() const{ return _nz; }

	/** \brief Access to the scalar of the vertex*/
	inline const double field() const{ return _field; }

public:
	/** \brief Sets the coordenate x of the vertex
	  * \param const double x */
	inline const void  set_x (const double x ) { _x=x  ; }
	
	/** \brief Sets the coordenate y of the vertex
    * \param const double y */
  inline const void  set_y (const double y ) { _y=y  ; }
	
	/** \brief Sets the coordenate z of the vertex
	  * \param const double z */
  inline const void  set_z (const double z ) { _z=z  ; }
	
	/** \brief Sets the coordenate nx of the normal
	  * \param const double nx */
  inline const void  set_nx(const double nx) { _nx=nx; }

	/** \brief Sets the coordenate ny of the normal
    * \param const double ny */
  inline const void  set_ny(const double ny) { _ny=ny; }

	/** \brief Sets the coordenate nz of the normal
	  * \param const double nz */
  inline const void  set_nz(const double nz) { _nz=nz; }

	/** \brief Sets the field the vertex normal
	  * \param const double f */
  inline const void  set_field(const double f) { _field=f; }

//-- Geometric Operations --//
public:
	/** \brief Computes the normal of a triangle
    * \param const Vertex &v0
    * \param const Vertex &v1
    * \param const Vertex &v2
    * \param float* n*/
  static void normal ( const Vertex &v0, const Vertex &v1, const Vertex &v2, double *n)
  {
      double a0[3], a1[3];
      n[0]= 0;
      n[1]= 0;
      n[2]= 0;

      a0[0]= v1.x()-v0.x();
      a0[1]= v1.y()-v0.y();
      a0[2]= v1.z()-v0.z();

      a1[0]= v2.x()-v0.x();
      a1[1]= v2.y()-v0.y();
      a1[2]= v2.z()-v0.z();

      n[0]= a0[1]*a1[2]-a1[1]*a0[2];
      n[1]= a0[2]*a1[0]-a1[2]*a0[0];
      n[2]= a0[0]*a1[1]-a1[0]*a0[1];
  }
  /** Computes the squared norm of a vector in R^3
    * \brief computes the squared norm of a vector in R^3
    * \param const float a
    * \param const float b
    * \param const float c*/
  inline static double norm2 ( const double a, const double b, const double c ) { return a*a + b*b + c*c ; }

  /** \brief computes the norm of a vector in R^3
    * \param const float a
    * \param const float b
    * \param const float c*/
  inline static double norm  ( const double a, const double b, const double c ) { return sqrt( norm2(a,b,c) ) ; }

  /** \brief Normalizes a vector in R^3
    * \param float* v*/
  static void normalize( double *v )
  {
    double n= norm (v[0], v[1], v[2]);
    if(n > FLT_EPSILON )
    {
      v[0]/=n;
      v[1]/=n;
      v[2]/=n;
    }
  }

//-- Vertex i/o operators --//
public:
  /** \brief Read operator for Vertex objects 
    * \param istream& s
    * \param Vertex& v */
  friend istream& operator >> (istream& s,  Vertex& v)
  {
    double x, y, z, nx, ny, nz, f;

    s >> x >> y >> z >> nx >> ny >> nz >> f ;

    v.set_x(x); v.set_y(y); v.set_z(z);
    v.set_nx(nx); v.set_ny(ny); v.set_nz(nz);
    v.set_field(f);
	  
	 return s; 
  }
  /** \brief Write operator for Vertex objects 
    * \param ostream& s
    * \param Vertex& v */
  friend ostream& operator << (ostream& s, const Vertex& v)
  {
    s << "  coord: ( " << v.x()     << " , " << v.y()  << " , " << v.z() << " ) " << endl ;
    s << " normal: ( " << v.nx()    << " , " << v.ny() << " , " << v.nz() << " ) " << endl ;
    s << "  field:   " << v.field() << endl;

	return s;
  }

//-- Vertex bool operators --//
public: 
  /** \brief Equal operator for Vertex objects 
    * \param Vertex& v */
   bool operator == (const Vertex& v){ return (x()==v.x() && y()==v.y() && z()==v.z() && nx()==v.nx() && ny()==v.ny() && nz()==v.nz() && field()==v.field());}
  
  /** \brief Equal operator for Vertex objects 
    * \param Vertex& v */
   bool operator != (const Vertex& v){ return (x()!=v.x() || y()!=v.y() || z()!=v.z() || nx()!=v.nx() || ny()!=v.ny() || nz()!=v.nz() || field()!=v.field());}
};
#endif
//----------------------------------------------------------//
