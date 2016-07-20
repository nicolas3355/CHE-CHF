/**
* @file    CHF_L2.hpp
* @author  Marcos Lage      <mlage@mat.puc-rio.br>
* @author  Thomas Lewiner   <thomas.lewiner@polytechnique.org>
* @author  Hélio  Lopes     <lopes@mat.puc-rio.br>
* @author  Math Dept, PUC-Rio
* @author  Lab Matmídia
* @date    14/02/2006
*
* @brief   CHF: A Scalable topological data structure for tetrahedral meshes
* @brief  (Level 2)
*/
//--------------------------------------------------//
#ifndef _CHF_L2_HPP_
#define _CHF_L2_HPP_

#include <map>
#include <ctime>
#include <vector>
#include <iostream>
#include "CHF_L1.hpp"

using namespace std;

/** \brief Edge id type */
typedef pair<Vid,Vid> Eid;

/** \brief Edge iterator*/
typedef map<Eid, HFid>::iterator Eit;
/** \brief Edge const iterator*/
typedef map<Eid, HFid>::const_iterator Ecit;

/** \brief Face Map iterator*/
typedef map<HFid, HFid>::iterator Fit;
/** \brief Face Map const iterator*/
typedef map<HFid, HFid>::const_iterator Fcit;
//--------------------------------------------------//
/** CHF data-structure for tetrahedral meshes
  * \brief CHF: Level 2*/
class CHF_L2:public CHF_L1
//--------------------------------------------------//
{
//-- CHF_L2 protected data.--//        
protected:
  /** \brief Half-Face of vertices*/
  vector<HFid>       _VH ;
  /** \brief Map of edges*/
  map<  Eid, HFid >  _EH ;
  /** \brief Map of faces*/
  map< HFid, HFid >  _FH ;

public:
  /** \brief Default constructor.*/
  CHF_L2():CHF_L1() {}

  /** \brief First constructor.
    * \param nv   -  const Vid.
    * \param ntet -  const TEid. */
  CHF_L2(const Vid nv, const TEid ntet): CHF_L1(nv,ntet) { _VH.resize( nvert(), -1 ); }

  /** \brief Copy constructor
    * \param h  -  const CHF_L2 object.*/
  CHF_L2(const CHF_L2& h): CHF_L1(h) { _EH=h._EH; _VH=h._VH; _FH=h._FH; }

  /** \brief Destructor.*/
  ~CHF_L2() { _FH.clear(); _EH.clear(); _VH.clear(); }

public:
  /** \brief Access to the half-face of a vertex. 
    * \param v - const Vid  */  
  inline const HFid  VH( const Vid v ) const { if( !v_valid(v) )  return INV;  return _VH[v] ; }

  /** \brief Access to the half-face of an edge. 
    * \param e - const Eid  */  
  inline const HFid *EH( const Eid e ) const { Ecit i = _EH.find( e ) ; if( i == _EH.end() ) return NULL ; return &i->second ; }

  /** \brief Access to the edge map.*/
  inline const map<Eid, HFid> &EH() const { return _EH; }

  /** \brief Access to one face of the map. 
    * \param h - const HFid  */  
  inline const HFid *FH( const HFid h ) const { Fcit i = _FH.find( h ) ; if( i == _FH.end() ) return NULL ; return &i->second ; }

  /** \brief Access to the face map.*/
  inline const map<HFid, HFid> &FH() const { return _FH; }

  /** \brief Tests if a vertex is valid
    * \param v - const Vid */
  inline const bool v_valid( const  Vid v ) const { return ( CHF_L1::v_valid(v) && _VH[v] != INV); }

  /** \brief Tests if an edge is valid
    * \param e - const Eid */
  inline const bool e_valid( const  Eid e ) const { Ecit i = _EH.find( e ) ;   return i != _EH.end() ; }

  /** \brief Tests if an face is valid
    * \param h - const HFid */
  inline const bool f_valid( const HFid f ) const { Fcit i = _FH.find( f ) ;   return i != _FH.end() ; }

public:
  /** \breaf Sets the half-face of a vertex. 
    * \param v - const Vid  
    * \param h - const HFid*/
  inline const void  set_VH( const Vid v, const HFid h ) { if( v_valid(v) ) _VH[v]=h; return ;}
  /** \brief Sets a vertex as invalid
    * \param v - const Vid */
  inline const void  v_invalid( const  Vid v ) { CHF_L1::v_invalid(v); _VH[v]=INV; }

public:
  /** \brief Computes the vertices in the star of a vertex 
    * \param v - const Vid */
  virtual vector<Vid>  R_00( const Vid v  );
  /** \brief Computes the tetrahedrons in the star of a vertex 
    * \param v - const Vid */
  virtual vector<TEid> R_03( const Vid v  );
  /** \brief Computes the vertices in the star of an edge 
    * \param a - const Vid  
    * \param b - const Vid */
  virtual vector<Vid>  R_10( const Vid a, const Vid b  );
  /** \brief Computes the tetrahedrons in the star of an edge 
    * \param a - const Vid  
    * \param b - const Vid */
  virtual vector<TEid> R_13( const Vid a, const Vid b  );

public:
  /** \brief creates the VH table. */
  void create_VH ();
  /** \brief creates the EH map. */
  void create_EH ();
  /** \brief creates the FH map. */
  void create_FH ();

  /** \brief Checks mesh validation*/
  void check ();

public:
  /** \brief Draws the vertices of the mesh. 
    * \param t= 0 - const int*/
  virtual void draw_vert ( const int t=1, const int p=false  );
  /** \brief Draws the mesh in wireframe. 
    * \param t= 0 - const int*/
  virtual void draw_wire ( const int t=4 );

public:
  /** \brief Reads the model from a file
    * \param fn - const char* */
  void  read_ply  ( const char* fn )
  { 
    // Stores Start && L1 && L2 time;
    clock_t start_time = static_cast<clock_t>(0.0),
               L2_time = static_cast<clock_t>(0.0);
    
    start_time = clock();
    
    CHF_L1::read_ply(fn); 

    create_VH(); 
    create_EH(); 
    create_FH();

    L2_time = clock();

    //cout << "L2 load time:" << static_cast<double>(L2_time-start_time)/static_cast<double>(CLOCKS_PER_SEC) << endl;
  }
};
#endif
