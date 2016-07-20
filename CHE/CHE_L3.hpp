/**
* @file    CHE_L3.hpp
* @author  Marcos Lage         <mlage@mat.puc-rio.br>
* @author  Thomas Lewiner      <thomas.lewiner@polytechnique.org>
* @author  Helio  Lopes        <lopes@mat.puc-rio.br>
* @author  Math Dept, PUC-Rio
* @author  Lab Matmidia
* @date    14/02/2006
*
* @brief  (Compact Half-Edge Structure - Level 3)
*/

#ifndef  _CHE_L3_HPP_
#define  _CHE_L3_HPP_

#include <map>
#include <vector>
#include <iostream>
#include "CHE_L2.hpp"

 /** \brief standart namespace definition*/
using namespace std;

/** CHE_L3 class
  *
  * Creates a structure for triangulated meshes,
  * storing the coordinates of its vertices (vector _G), 
  * the indices of the vertices of each triangle (vector _V),
  * the relations between the triangles of the mesh (vector _O),
  * the connected compound of each vertex.(vector _C),
  * a half-edge for each vertex (vector _VH)
  * the edges of the model (map _EH),
  * and the boundary curves (vector _CH)
  * 
  * The class inherits the informations of CHE_L2.*/
class CHE_L3:public CHE_L2
{
protected:
  /** \brief Boundary Curves Table: 
    * For each boundary curve we store one 
    * representative half-edge*/
  vector<HEid>  _CH;

  /** \brief Number Boundary curves
    *
    * Protected data that stores the number
    * of boundary curves of the mesh*/
  Cid   _ncurves;
public:
  /** \brief Default constructor.*/
  CHE_L3(): CHE_L2() {}

  /** \brief First constructor.
    * \param nvert  - Vid  CHE_L0 _nvert.
    * \param ntrig  - TRid CHE_L0 _ntrig.*/
  CHE_L3(Vid nvert, TRid ntrig): CHE_L2(nvert, ntrig){ }

  /** \brief Copy constructor
    * \param c - const CHE_L3&.*/
  CHE_L3(const CHE_L3& c): CHE_L2( c ) { _CH= c._CH; }

   /** \brief Destructor.*/
  virtual ~CHE_L3(){ _V.clear(); _G.clear(); _O.clear(); _C.clear(); _VH.clear(); _EH.clear(); _CH.clear(); }

public:
	/** \brief Access to the number of boundary curves of the model*/
	inline const  Cid ncurves() const{ return _ncurves; }

  /** \brief Access the representative of a boundary compound  
    * \param const Vid v*/
   inline HEid CH( const Vid v ) const{if( !v_valid(v) ) return INV; return _VH[v]; }
 	
   /** \brief Tests if a boundary is valid
     * \param const Vid v*/
	inline const bool b_valid( const  Cid b ) const {return (b >= 0 && b<ncurves()); }

public:
	/** \brief Sets to the number of boundary curves of the model
       * \param nc - Cid*/
  inline const void set_ncurves( Cid nc ){ _ncurves = nc; }

  /** \brief Sets an half-edge for a vertex  
     * \param b - const Cid
     * \param h - const HEid */
   inline void set_CH( const Cid b, const HEid h ){ if( b_valid(b) ) _CH[b]=h; return; }

   /** \brief Sets a boundary curve invalid
      * \param b - const Cid  */
   inline const void b_invalid( const Cid b ){ if( b_valid(b) ) _CH[b] = INV; return; }

public:
  /** \brief Computes the Boundary Curves table*/
  void compute_CH();
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
//-----------------------------------------------------------------------//
