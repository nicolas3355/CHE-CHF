/**
* @file    CHF_L3.hpp
* @author  Marcos Lage      <mlage@mat.puc-rio.br>
* @author  Thomas Lewiner   <thomas.lewiner@polytechnique.org>
* @author  Hélio  Lopes     <lopes@mat.puc-rio.br>
* @author  Math Dept, PUC-Rio
* @author  Lab Matmídia
* @date    14/02/2006
*
* @brief   CHF: A Scalable topological data structure for tetrahedral meshes
* @brief  (Level 3)
*/
//--------------------------------------------------//
#ifndef _CHF_L3_HPP_
#define _CHF_L3_HPP_

#include <ctime>
#include <vector>
#include "CHF_L2.hpp"

using namespace std;

/** \brief Boundary  id type*/
typedef int  Bid;
/** \brief Triangle  id type*/
typedef int TRid;
/** \brief Half-Edge  id type*/
typedef int HEid;
//--------------------------------------------------//
/** CHF data-structure for tetrahedral meshes
  * \brief CHF: Level 3*/
class CHF_L3:public CHF_L2
{
//--------------------------------------------------//
//-- CHF_L3 protected data.--//        
protected:
  /** \brief number of boundary surfaces*/
  Bid  _bnsurf;
  
  /** \brief number of boundary triangles*/
  TRid _bntrig;

  /** \brief Boundary Vertex Table */
  vector<Vid> _bV;
  
  /** \brief Boundary Opposite Table */
  vector<Vid> _bO;
	
  /** \brief Vertex Surface Table */
  vector<Bid> _bS;

public:
  /** \brief Default constructor.*/
  CHF_L3():CHF_L2() {}

  /** \brief First constructor.
    * \param nv   -  const Vid.
    * \param ntet -  const TEid. */
  CHF_L3(const Vid nv, const TEid ntet): CHF_L2(nv,ntet) { _bS.resize( nvert(), -1 ); }

  /** \brief Copy constructor
    * \param h  -  const CHF_L2 object.*/
  CHF_L3(const CHF_L3& h): CHF_L2(h) { _bV=h._bV; _bO=h._bO; _bS= h._bS; _bnsurf = h._bnsurf; _bntrig= h._bntrig; }

  /** \brief Destructor.*/
  ~CHF_L3() { _bV.clear(); _bO.clear(); _bS.clear(); _bnsurf=0; _bntrig=0; }

public:
  /** \brief Accesses the number of boundary surfaces  of the model  */
  inline const  Bid bnsurf() const{ return _bnsurf; } 
	
  /** \brief Accesses the number of boundary triangles of the model */
  inline const TRid bntrig() const{ return _bntrig; }

  /** \brief Accesses the vertex of a boundary half-edge 
    * \param h - const HEid */
  inline const  Vid bV ( const HEid h ) const { if( !he_valid(h) ) return INV;  return _bV[h] ; }
  /** \brief Accesses the opposite of a boundary half-edge 
    * \param h - const HFid */
  inline const HEid bO ( const HEid h ) const { if( !he_valid(h) ) return INV;  return _bO[h] ; }
  /** \brief Accesses the boundary surface os a vertex 
    * \param h - const HFid */
  inline const  Bid bS ( const  Vid v ) const { if(  !v_valid(v) ) return INV;  return _bS[v] ; }

  /** \brief Tests if a boundary half-edge is valid
    * \param h - const HEid */
  inline const bool he_valid ( const HEid h ) const { return ( h >= 0 && h < 3*bntrig() && !(_bO[h] == INV) && !(_bV[h] == INV) ); }
  /** \brief Tests if a boundary triangle is valid
    * \param t - const TRid */
  inline const bool tr_valid ( const TRid t ) const { return ( t >= 0 && t < bntrig() && he_valid(3*t) && he_valid(3*t+1) && he_valid(3*t+2) ); }

protected:
  /** \brief Access to the vertex of a boundary half-edge. 
    * \param h - const HFid  */  
  inline void  set_bV ( const HEid h, const  Vid v )  { if( he_valid(h) ) _bV[h]=v ; return; }
  /** \brief Access to the opposite of a boundary half-edge. 
    * \param h - const HFid  */  
  inline void  set_bO ( const HEid h, const HEid m )  { if( he_valid(h) ) _bO[h]=m ; return; }
  /** \brief Access to the boundary surface of a vertex. 
    * \param h - const HFid  */  
  inline void  set_bS ( const  Vid v, const  Bid b )  { if(  v_valid(v) ) _bS[v]=b ; return; }

  /** \brief Sets a boundary half-edge as invalid
    * \param v - const Vid */
  inline void  he_invalid ( const HEid h )  { if( h >= 0 && h< 3*bntrig() ) {_bO[h] = INV; _bV[h] = INV; } }
  /** \brief Sets a boundary triangle as invalid
    * \param v - const Vid */
  inline void  tr_invalid ( const TRid t )  { if( t >= 0 && t<   bntrig() ) { he_invalid(3*t); he_invalid(3*t+1); he_invalid(3*t+2); } }

public:
  /** \brief Accesses the triangle of a half-edge 
    * \param h - const HEid*/
  inline const TRid btrig   (const HEid h) const { if( !he_valid(h) ) return INV ; return h/3 ; }
  /** \brief Accesses the next of a boundary half-edge 
    * \param h - const HEid*/
  inline const HEid bnexthe (const HEid h) const { if( !he_valid(h) ) return INV ; return 3*(h/3) + (h+1) % 3 ; }
  /** \brief Accesses the previous of a boundary half-edge 
    * \param h - const HEid*/
  inline const HEid bprevhe (const HEid h) const { if( !he_valid(h) ) return INV ; return 3*(h/3) + (h+2) % 3 ; }

public:
  /** \brief Creates the bS table. */
  void create_bS () ;
  /** \brief Creates the bV table. */
  void create_bV () ;
  /** \brief Creates the bO table. */
  void create_bO () ;

  /** \brief Computes bound faces' normal.*/  
  virtual void compute_normals();

private :
  /** \brief Sets bS and returns the value*/
  Bid  get_bS( Bid i ) ;

public :
  /** \brief Checks the orientation beetwen two boundary triangles. 
    * \param c - HEid 
    * \param t - HEid  */
  const bool borient_check(const HEid c, const HEid t);
  /** \brief Checks mesh validation*/
  void check ();

public:
  /** \brief Draws the bound surface of the mesh. 
    * \param t= 0 - const int*/
  virtual void draw_smooth ( const int t=0 );

public:
  /** \brief Reads the model from a file
    * \param fn - const char* */
  void read_ply ( const char* fn )
  {
    // Stores Start && L2 && L3 time;
    clock_t start_time = static_cast<clock_t>(0.0),
                     L3_time = static_cast<clock_t>(0.0);

    start_time = clock();

    CHF_L2::read_ply(fn); 
 
    create_bV(); 
    create_bO(); 
    create_bS(); 
    CHF_L3::compute_normals();

    L3_time = clock();

    //cout << "L3 load time:" << static_cast<double>(L3_time-start_time)/static_cast<double>(CLOCKS_PER_SEC) << endl;
  }
};
#endif
