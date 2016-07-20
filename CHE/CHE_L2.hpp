/**
* @file    CHE_L2.hpp
* @author  Marcos Lage         <mlage@mat.puc-rio.br>
* @author  Thomas Lewiner      <thomas.lewiner@polytechnique.org>
* @author  Helio  Lopes        <lopes@mat.puc-rio.br>
* @author  Math Dept, PUC-Rio
* @author  Lab Matmidia
* @date    14/02/2006
*
* @brief  (Compact Half-Edge Structure - Level 2)
*/

#ifndef _CHE_L2_HPP_
#define _CHE_L2_HPP_

#include <map>
#include <vector>
#include <iostream>
#include "Edge.hpp"
#include "CHE_L1.hpp"

/** \brief Edge iterator */
typedef  map<HEid, HEid>::iterator Eit;
/** \brief Edge constant iterator */
typedef  map<HEid, HEid>::const_iterator Ecit;

/** \brief standart namespace definition*/
using namespace std;
/** CHE_L2 class
  *
  * Creates a structure for triangulated meshes,
  * storing the coordinates of its vertices (vector _G), 
  * the indices of the vertices of each triangle (vector _V),
  * the relations between the triangles of the mesh (vector _O),
  * the connected compound of each vertex.(vector _C),
	* a half-edge for each vertex (vector _VH)
	* and the edges of the model (map _EH).
  * 
  * The class inherits the informations of CHE_L1.*/
class CHE_L2:public CHE_L1
{
protected:
  /** \brief Vertex Half-Edge Table: For each vertex 
     * we store a half-edge associated*/
  vector<HEid>    _VH;

  /** \brief Edge Table: Map with the edges of the model */
  map<HEid,HEid>  _EH;

public:
  /** \brief Default constructor.*/
  CHE_L2(): CHE_L1() {}

  /** \brief First constructor.
    * \param nvert  - Vid  CHE_L0 _nvert.
    * \param ntrig  - TRid CHE_L0 _ntrig.*/
  CHE_L2(Vid nvert, TRid ntrig): CHE_L1(nvert, ntrig){ _VH.resize( nvert, -1 ); }

  /** \brief Copy constructor
    * \param c - const CHE_L2&.*/
  CHE_L2(const CHE_L2& c): CHE_L1( c ) { _EH= c._EH; _VH=c._VH; }

   /** \brief Destructor.*/
  virtual ~CHE_L2(){ _V.clear(); _G.clear(); _O.clear(); _C.clear(); _VH.clear(); _EH.clear(); }

public:
   /** \brief Access to the edge table of the model*/
   inline const map<HEid, HEid> &EH() const { return _EH ; }
   /** \brief Access to an edge of the model  
     * \param  e - const HEid*/
   inline Ecit EH( const HEid e ) const{ return _EH.find(e); }
   /** \brief Access a half-edge of a vertex  
     * \param const Vid v    */
   inline HEid VH( const Vid v ) const{if( !v_valid(v) ) return INV; return _VH[v]; }
   /** \brief Tests if a vertex is valid
     * \param  v - const Vid */
   inline const bool v_valid( const  Vid v ) const {return (CHE_L1::v_valid( v ) && _VH[v]!= INV); }
   /** \brief Tests if a vertex is on bound
     * \param v - const Vid  */
   inline const bool v_bound( const  Vid v ) const {return ( O(VH(v))==-1); }

public:
   /** \brief Sets an edge of the model  
     * \param e  - const HEid 
     * \param ed - const HEid */
   inline void set_EH( const HEid e, const HEid ed ) { _EH[e] = ed; }
   /** \brief Sets an half-edge for a vertex  
     * \param v - const Vid
     * \param h - const HEid */
   inline void set_VH( const Vid v, const HEid h ){ if( v_valid(v) ) _VH[v] = h; return; }
   /** \brief Sets a vertex invalid
     * \param v - const Vid  */
   inline const void v_invalid( const  Vid v ){ if( v_valid(v) ){CHE_L0::v_invalid(v); _VH[v] = INV;} return; }

public:
  /** \brief Computes the vertices in the star of a given vertex 
	  * \param v - const Vid  */
  virtual vector<Vid>  R_00( const Vid v  );
  /** \brief Computes the triangles in the star of a given vertex 
	  * \param v - const Vid  */
  virtual vector<TRid> R_02( const Vid v  );

public:
  /** \brief Computes the Edge table*/
  void compute_EH();
  /** \brief Computes the Half-edge table*/
  void compute_VH();
  /** \brief Checks the mesh*/
  void check();
	/** \brief Draws the surface in wireframe with opengl*/
	virtual void draw_wire() ;

public:
	/** \brief Reads a 3D model in the .ply format
	  * \param file - const char* */
  void  read_ply( const char* file );
};
#endif
//-------------------------------------//
