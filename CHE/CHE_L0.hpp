/**
* @file    CHE_L0.hpp
* @author  Marcos Lage         <mlage@mat.puc-rio.br>
* @author  Thomas Lewiner      <thomas.lewiner@polytechnique.org>
* @author  Helio  Lopes        <lopes@mat.puc-rio.br>
* @author  Math Dept, PUC-Rio
* @author  Lab Matmidia
* @date    14/02/2006
*
* @brief  (Compact Half-Edge Structure - Level 0)
*/

#ifndef _CHE_L0_HPP_

#define _CHE_L0_HPP_



#include <vector>

#include <cfloat>

#include <iostream>

#include "Vertex.hpp"



/** \brief Invalid integer*/
#define INV -10000
/** \brief Pet_triangle id*/

typedef int  TRid; 

/** \brief Pet_half-edge id*/

typedef int  HEid;



/** \brief Invalid vertice*/
static Vertex V_INV( FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX );



/** \brief standart namespace definition*/

using namespace std;



/** \brief CHE_L0 Structure

  *

  * Creates a Soup of triangles, storing the

  * coordinate of the vertices of a triangulated

  * mesh (vector _G) and the indices os each 

  * triangle (vector _V).*/

class CHE_L0

{

protected:

  /** \brief Number of vertices in the mesh

    *

    * Protected data that stores the number

    * of vertices of the mesh*/

  Vid   _nvert;



  /** \brief Number of triangles in the mesh     

    *

    * Protected data that stores the number

    * of triangles of the mesh*/

  TRid  _ntrig;

  

  /** \brief Geometry table 

    * 

    * Stores the coordinates of the vertices 

    * of the mesh*/

  vector<Vertex> _G;

  

  /** \brief Vertices table 

    * 

    * Stores the indice of the vertices 

    * of the triangles of the mesh*/

  vector<Vid>    _V;





public:

  /** \brief Default constructor.

    * 

    * Constructs a CHE_L0 object with

    * _nvert= 0 and _ntrig= 0.*/

  CHE_L0(): _nvert(0), _ntrig(0) {}

  

  /** \brief First constructor.

    *

    * Constructs a CHE_L0 object with

    * _nvert = nvert and _ntrig = ntrig .

    *

    * \param nvert - Vid  CHE_L0 _nvert.

    * \param ntrig - TRid CHE_L0 _ntrig.*/

  CHE_L0(Vid nvert, TRid ntrig): _nvert(nvert), _ntrig(ntrig){ _V.resize( 3*ntrig, -1 ); _G.resize( nvert ); }

  

  /** \brief Copy constructor

    *

    * Constructs a CHE_L0 object equal to c.

    *

    * \param c - CHE_L0&.*/

  CHE_L0(const CHE_L0& c): _nvert( c.nvert() ), _ntrig( c.ntrig() ) { _V= c._V; _G=c._G; }



  /** \brief Destructor.

    *

    * Frees the memory used by _V and _G*/

  virtual ~CHE_L0(){ _V.clear(); _G.clear(); }



public:

	/** \brief Access to the number of vertices of the model*/

	inline const  Vid nvert() const{ return _nvert; }

	

	/** \brief Access to the number of triangles of the model*/

	inline const TRid ntrig() const{ return _ntrig; }



 	/** \brief Tests if a vertex is valid

    * \param const Vid v*/

	inline const bool  v_valid( const  Vid v )  const { return (v >= 0 && v<nvert() && (Vertex)_G[v]!= V_INV); }

 	

  /** \brief Tests if a half-edge is valid

    * \param const HEid h*/

	inline const bool he_valid( const  HEid h ) const { return (h >= 0 && h<3*ntrig() && _V[h] != INV); }



  /** \brief Tests if a triangle is valid

    * \param const HEid h*/

	inline const bool tr_valid( const  TRid t ) const { return (he_valid(3*t) && he_valid(3*t+1) && he_valid(3*t+2)); }



  /** \breaf Access to the geometry of a vertex in the model 

	  * \param const Vid v*/

	inline const Vertex G( const Vid  v ) const { if( !v_valid(v)  ) return V_INV;  return _G[v] ; }

  

  /** \breaf Access to the geometry of a vertex in the model 

	  * \param const Vid v*/

	inline Vertex &G( const Vid  v ) { if( !v_valid(v)  ) return V_INV;  return _G[v] ; }



  /** \brief Access to the vertex index of a half-edge

    * \param const Vid v*/

	inline const    Vid V( const HEid h ) const { if( !he_valid(h) ) return   INV;  return _V[h] ; }



public:

	/** \brief Sets to the number of vertices in the model  

    * \param nvert - Vid */

	inline const void set_nvert( Vid nvert ){ _nvert=nvert; }

	

	/** \brief Sets to the number of triangles in the model  

    * \param ntrig - TRid*/

  inline const void set_ntrig( TRid ntrig ){ _ntrig=ntrig; }



 	/** \brief Makes a vertex invalid

    * \param const Vid v */

	inline const void v_invalid ( const   Vid v ){ if( v_valid(v) ) _G[v]= V_INV; return; }

 	

  /** \brief Makes a half-edge invalid

    * \param const HEid h  */

	inline const void he_invalid( const  HEid h ){ if( he_valid(h) ) _V[h]= INV; return;}



  /** \brief Makes a triangle invalid

    * \param t - const TRid */

  inline const void tr_invalid( const  TRid t ){ if( tr_valid(t) ){he_invalid(3*t); he_invalid(3*t+1); he_invalid(3*t+2);} return;}



  /** \brief Sets the geometry of a vertex in the model 

	  * \param const Vid v  

    * \param Vertex p*/

	inline const void set_G( const Vid  v, Vertex p ) { if(v>=0 && v<nvert()) _G[v]=p ; }



 	/** \brief Sets the vertex of a half-edge

    * \param const HEid v  

    * \param Vid v*/

	inline const void set_V( const  HEid h, Vid v ) { if( h>=0 && h<3*ntrig()) _V[h]=v ; }



public:

	/** \brief Access to the triangle of a given half-edge 

	  * \param const HFid h */

	inline const TRid trig ( const HEid h ) const { if( !he_valid(h) ) return INV ; else return h/3 ; }

	

	/** \brief Access to the half-edge "next" of a given half-edge 

	  * \param const HFid h */

	inline const HEid next ( const HEid h ) const { if( !he_valid(h) ) return INV ; else return 3*(h/3) + (h+1) % 3 ; }

	

	/** \brief Access to the half-edge "previous" of a given half-edge 

	  * \param const HFid h */

	inline const HEid prev ( const HEid h ) const { if( !he_valid(h) ) return INV ; else return 3*(h/3) + (h+2) % 3 ; }



public:

	/** \brief Computes the vertices in the star of a given vertex 

    * \param v - const Vid  */

	virtual vector<Vid>  R_00( const Vid v  );

	/** \brief Computes the triangles in the star of a given vertex 

	  * \param v - const Vid  */

	virtual vector<TRid> R_02( const Vid v  );

	/** \brief Computes the vertices in the star of a given edge 

	  * \param h - const HEid */

	virtual vector<Vid>  R_10( const HEid h );

	/** \brief Computes the triangles in the star of a given edge 

	  * \param h - const HRid */

	virtual vector<TRid> R_12( const HEid h );

	/** \brief Computes the triangles adjacents of a given triangle 

	  * \param r - const TRid */

	virtual vector<TRid> R_22( const TRid t );



public:

	/** \brief Computes the normal of each triangle*/

	void compute_normals () ;

  /** \brief Gets the model bouding_box

    * \param min - float*.

		* \param max - float*.*/

	void bounding_box  ( float *min, float *max );

  /** \brief Legalizes the model

    * \param min - float*.

		* \param max - float*.*/

	void legalize_model( float *min, float *max );



	/** \brief Checks the mesh*/

	void check () ;



  /** \brief Draws the smooth surface with opengl*/

	void draw_smooth() ;

	/** \brief Draws the surface in wireframe with opengl*/

	virtual void draw_wire() ;

	/** \brief Draws the vertices of the surface with opengl*/

	void draw_verts () ;



public:

	/** \brief Reads a 3D model from the .ply format

	  * \param file - const char* */

  void  read_ply( const char* file );

	/** \brief Writes a 3D model in the .ply format

	  * \param file - const char* 

	  * \param bin= false - bool*/

  void write_ply( const char* file, bool bin=false );

};

#endif
