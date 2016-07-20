/**
* @file    CHF_L1.hpp
* @author  Marcos Lage      <mlage@mat.puc-rio.br>
* @author  Thomas Lewiner   <thomas.lewiner@polytechnique.org>
* @author  Hélio  Lopes     <lopes@mat.puc-rio.br>
* @author  Math Dept, PUC-Rio
* @author  Lab Matmídia
* @date    14/02/2006
*
* @brief   CHF: A Scalable topological data structure for tetrahedral meshes
* @brief  (Level 1)
*/
//--------------------------------------------------//
#ifndef _CHF_L1_HPP_
#define _CHF_L1_HPP_

#include <map>
#include <ctime>
#include <algorithm>
#include "CHF_L0.hpp"

/** \brief standart namespace definiton*/
using namespace std;
//--------------------------------------------------//
/** CHF data-structure for tetrahedral meshes
  * \brief CHF: Level 1*/
class CHF_L1:public CHF_L0 
//--------------------------------------------------//
{
//-- CHF_L1 protected data.--//        
protected:
  /** \brief Opposite container */
  vector<HFid> _O;

public:
  /** \brief Default constructor.*/
  CHF_L1():CHF_L0() {}

  /** \brief First constructor.
    * \param nv   -  const Vid.
    * \param ntet -  const TEid. */
  CHF_L1(const Vid nv, const TEid ntet): CHF_L0(nv,ntet) {_O.resize(ntetra()<<2); }

  /** \brief Copy constructor
    * \param h  -  const CHF_L1 object.*/
  CHF_L1(const CHF_L1& h): CHF_L0(h) { _O=h._O; }

  /** \brief Destructor.*/
  ~CHF_L1() { _O.clear(); }

public:
 /** \brief Access to the opposite of a half-face. 
   * \param h - const HFid  */  
  inline const HFid  O( const HFid h ) const { if( !hf_valid(h) )  return INV;  return _O[h] ; }
	
 /** \brief Tests if a tetrahedron is valid
   * \param t - const TEid */
  inline const bool  te_valid( const TEid t ) const { return ( t >= 0 && t < ntetra() && hf_valid(t<<2) && hf_valid(t<<2 | 1) && hf_valid(t<<2 | 2)  && hf_valid(t<<2 | 3) ) ;  }
	
  /** Tests if a half-face is valid
    * \param h - const HFid */
  inline const bool  hf_valid( const HFid h ) const { return (CHF_L0::hf_valid(h) && _O[h] != INV); }

public:
  /** \breaf Sets the opposite of a half-face 
    * \param h - const HFid  
    * \param o - HFid */
  inline const void set_O( const HFid h, HFid o ) { if( hf_valid(h) )_O[h]=o ; }

  /** \brief Sets a tetrahedron as invalid
    * \param t - const TEid */
  inline const void te_invalid( const  TEid t ){ if( te_valid(t) ){hf_invalid(t<<2); hf_invalid(t<<2 | 1); hf_invalid(t<<2 | 2); hf_invalid(t<<2 | 3);} return;}

  /** \brief Sets a half-face as invalid
    * \param h - const HEid */
  inline const void hf_invalid( const  HFid h ){ if( hf_valid(h) ){CHF_L0::hf_invalid(h); _O[h]= INV;} return;}
  
  /** \brief Accesses the mate of a half-edge he in a half-face h 
    * \param h  - const HFid
    * \param he - const HFid*/
  inline const pair<HFid,HFid> radialhe( const HFid hf, const HFid he ) const 
  { 
    if( !hf_valid(hf) || !hf_valid(he) || (hf>>2) != (he>>2) || he == hf ) { cout<< "CHF ERRO: radialhe" << endl; return make_pair(INV,INV); }
    
    HFid n = nexthe(hf, he).second , o = O(hf);
    for(HFid i= 4*tetra(o); i< 4*(tetra(o)+1); ++i )
    if( V(i) == V(n) ) { n = i; break; }
   return ( make_pair(o, n) );
  }

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
  /** \brief Computes the tetrahedrons in the star of a tetrahedron 
    * \param r - const TEid */
  virtual vector<TEid> R_33( const TEid t );

  /** \brief creates the O table. */
  void create_O();

  /** \brief Computes bound faces' normal.*/  
  virtual void compute_normals();

protected:
  /** \brief Changes the orientation of a tetrahedron
    * \param t - TEid */
  void change_orientation(TEid t) ;

public:
  /** \brief Checks the orientation beetwen two tetrahedrons. 
    * \param c - HFid 
    * \param t - HFid  */
  const bool orient_check(const HFid c, const HFid t) ;

  /** \brief Orients the mesh.*/  
  void  orient();

  /** \brief Checks mesh validation*/
  void   check(); 

public:
  /** \brief Draws the bound surface of the mesh. 
    * \param t= 0 - const int*/
	virtual void draw_smooth ( const int t=0 );
    
public:
  /** \brief Reads the model from a file
    * \param fn - const char* */
  void  read_ply ( const char* fn )
  { 
    // Stores Start && L1 time;
    clock_t start_time  = static_cast<clock_t>(0.0),
                      L1_time = static_cast<clock_t>(0.0);

    start_time = clock();

    CHF_L0::read_ply(fn); 
 
    create_O(); 
    CHF_L1::compute_normals();

    L1_time = clock();

    // cout << "L1 load time:" << static_cast<double>(L1_time-start_time)/static_cast<double>(CLOCKS_PER_SEC) << endl;
  }
};
#endif
