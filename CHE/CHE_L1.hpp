/**
* @file    CHE_L1.hpp
* @author  Marcos Lage         <mlage@mat.puc-rio.br>
* @author  Thomas Lewiner      <thomas.lewiner@polytechnique.org>
* @author  Helio  Lopes        <lopes@mat.puc-rio.br>
* @author  Math Dept, PUC-Rio
* @author  Lab Matmidia
* @date    14/02/2006
*
* @brief  (Compact Half-Edge Structure - Level 1)
*/

#ifndef _CHE_L1_HPP_
#define _CHE_L1_HPP_

#include <map>
#include <vector>
#include <iostream>
#include "CHE_L0.hpp"

/** \brief Connected Compound id type */
typedef int  Cid;

/** \brief standart namespace definition*/
using namespace std;

/** CHE_L1 class
  *
  * Creates a structure for triangulated meshes,
  * storing the coordinates of its vertices (vector _G), 
  * the indices of the vertices of each triangle (vector _V),
  * the relations between the triangles of the mesh (vector _O)
  * and the connected compound of each vertex.(vector _C)
  * 
  * The class inherits the informations of CHE_L0.*/
class CHE_L1:public CHE_L0
{
protected:
  /** \brief Number of connected compounds 
    *
    * Protected data that stores the number of 
    * connected componds of the mesh*/
  Cid  _ncomp;

  /** \brief Opposite table 
    * 
    * For each half-edge we store the other half-edge 
    * related to the same edge, with opposite orientation*/ 
  vector<HEid> _O;

  /** \brief Connected compound Table
    * 
    * For each vertex we store its connected compound */
  vector<Cid>  _C;

public:
  /** \brief Default constructor.
    *
    * Constructs a CHE_L1 object with 
    * _nvert = 0, _ntrig = 0 and _ncomp = 0.*/
  CHE_L1(): CHE_L0(), _ncomp(0) {}

  /** \brief First constructor.
    * 
    * Constructs a CHE_L1 object with
    * _nvert = nvert, _ntrig = ntrig and _ncomp = 0.
    * 
    * \param nvert -  Vid  CHE_L0 _nvert.
    * \param ntrig - TRid  CHE_L0 _ntrig.*/
  CHE_L1(Vid nvert, TRid ntrig): CHE_L0(nvert, ntrig), _ncomp(0){ _O.resize( 3*ntrig, -1 ); _C.resize( nvert, -1 ); }
  
  /** \brief Copy constructor
    *
    * Constructs a CHE_L1 object equal to c.
    *
    * \param c - const CHE_L1&.*/
  CHE_L1(const CHE_L1& c): CHE_L0( c ), _ncomp( c.ncomp() ) { _O = c._O; _C = c._C; }

  /** \brief Destructor.
    *
    * Frees the memory used by _V, _G, _O, _C*/
  virtual ~CHE_L1(){ _V.clear(); _G.clear(); _O.clear(); _C.clear(); }

public:
	/** \brief Access to the number of connected compounds of the model*/
	inline const  Cid ncomp() const{ return _ncomp; }

 	/** \brief Access to the "opposite" of a half-edge
    * \param  h - const HEid*/
	inline const HEid O( const HEid h ) const {  if( !he_valid(h) ) return INV ;  return _O[h] ; }
 	
  /** \brief Access to the boundary compound of a vertex
    * \param  v - const Vid*/
	inline const  Cid C( const  Vid v ) const {  if( !v_valid(v)  ) return INV ;  return _C[v] ; }

 	/** \brief Tests if a vertex is valid
    * \param  v - const Vid*/
  inline const bool  v_valid( const   Vid v ) const { return (CHE_L0::v_valid( v )  && _C[v]!= INV); }
 	 
  /** \brief Tests if a half-edge is valid
    * \param  h - const HEid*/
	inline const bool he_valid( const  HEid h ) const { return (CHE_L0::he_valid( h ) && _O[h]!= INV); }
  
  /** \brief Tests if a triangle is valid
    * \param  t - const TRid*/
	inline const bool tr_valid( const  TRid t ) const { return (he_valid(3*t) && he_valid(3*t+1) && he_valid(3*t+2)); }

public:
	/** \brief Sets the number of connected compounds of the model
    * \param ncomp - Cid*/
  inline const void set_nbound( Cid ncomp ){ _ncomp=ncomp; }

  /** \breaf Sets the opposite of a half_edge of the model 
	  * \param h - const HEid  
    * \param o - HEid */
	inline const void set_O( const HEid h, HEid o ) { if(h>=0 && h<3*ntrig()) _O[h]=o ; }

 	/** \brief Sets the bound of a vertex
      * \param v - const Vid  
      * \param b - Cid */
	inline const void set_C( const Vid v, Cid b ) {  if( v>=0 && v<nvert()) _C[v]=b ; }

  /** \brief Makes a vertex invalid
    * \param v - const Vid */
  inline const void  v_invalid( const   Vid v ){ if(  v_valid(v) ){CHE_L0::v_invalid(v);  _C[v]= INV;} return; }
 	
  /** \brief Makes a half-edge invalid
    * \param h - const HEid */
  inline const void he_invalid( const  HEid h ){ if( he_valid(h) ){CHE_L0::he_invalid(h); _O[h]= INV;} return;}
  
  /** \brief Makes a triangle invalid
    * \param t - const TRid */
  inline const void tr_invalid( const  TRid t ){ if( tr_valid(t) ){he_invalid(3*t); he_invalid(3*t+1); he_invalid(3*t+2);} return;}

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
	  * \param h - const HEid */
  virtual vector<TRid> R_12( const HEid h );
  /** \brief Computes the triangles in the star of a given triangle 
	  * \param r - const TRid*/
  virtual vector<TRid> R_22( const TRid t );

public:
	/** \brief Computes the opposite of each half-edge*/
  void compute_opposites();
 	/** \brief Computes the connected compound of each vertex*/
  void compute_connected();
	/** \brief Checks the mesh*/
 void check () ;

private:
  /** \brief Orients the mesh*/
  void orient();

  /** \brief Changes the orientation of a triangle
	  * \param t - const TRid */
  void orient_change(const TRid t);

  /** \brief Checks the orientation between two half-edges
	  * \param h - const HEid 
	  * \param o - const HEid */ 
  const bool orient_check(const HEid h, const HEid o);
	
  /** \brief Gets the connected compound
    * \param i - Cid*/
  Cid get_component( Cid i ) ;

public:
  /** \brief Reads a 3D model in the .ply format
    * \param file - const char* */
  void  read_ply( const char* file );
};
#endif
//-----------------------------------------------//
