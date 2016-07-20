//------------------------------------------------
// ColorRamp creation
//------------------------------------------------
//
// Colors for Morse function visualization
// Version 0.1 - 23/10/2002
//
// Thomas Lewiner: thomas.lewiner@polytechnique.org
// INRIA - Geometrica Project
//
//________________________________________________



#ifndef _COLORRAMP_H_
#define _COLORRAMP_H_

//#pragma interface


//_____________________________________________________________________________
// ColorMaps available
enum ColorMap
//-----------------------------------------------------------------------------
{
  COLOR_NONE = 0 ,
  COLOR_GRAYSCALE,
  COLOR_RAINBOW,
  COLOR_AUTUMN, // varies smoothly from red, through orange, to yellow.
  COLOR_BONE, // is a grayscale colormap with a higher value for the blue component. This colormap is useful for adding an "electronic" look to grayscale images.
  COLOR_COLORCUBE, // contains as many regularly spaced colors in RGB colorspace as possible, while attempting to provide more steps of gray, pure red, pure green, and pure blue
  COLOR_COOL, // consists of colors that are shades of cyan and magenta. It varies smoothly from cyan to magenta.
  COLOR_COPPER, // varies smoothly from black to bright copper.
  COLOR_FLAG, // consists of the colors red, white, blue, and black. This colormap completely changes color with each index increment.
  COLOR_HOT, // varies smoothly from black, through shades of red, orange, and yellow, to white.
  COLOR_HSV, // varies the hue component of the hue-saturation-value color model. The colors begin with red, pass through yellow, green, cyan, blue, magenta, and return to red. The colormap is particularly appropriate for displaying periodic functions. hsv(m) is the same as hsv2rgb([h ones(m,2)]) where h is the linear ramp, h = (0:m-1)'/m.
  COLOR_JET, // ranges from blue to red, and passes through the colors cyan, yellow, and orange. It is a variation of the hsv colormap. The jet colormap is associated with an astrophysical fluid jet simulation from the National Center for Supercomputer Applications.
  COLOR_LINES, // produces a colormap of colors specified by the axes ColorOrder property and a shade of gray.
  COLOR_PINK, // contains pastel shades of pink. The pink colormap provides sepia tone colorization of grayscale photographs.
  COLOR_PRISM, // repeats the six colors red, orange, yellow, green, blue, and violet.
  COLOR_SPRING, // consists of colors that are shades of magenta and yellow.
  COLOR_SUMMER, // consists of colors that are shades of green and yellow.
  COLOR_WINTER, // consists of colors that are shades of blue and green.
  COLOR_NUMBER
};
//_____________________________________________________________________________




//_____________________________________________________________________________
// Application
class ColorRamp
//-----------------------------------------------------------------------------
{
private :
  float  color_random   [256*3] ;
  float  color_grayscale[256*3] ;
  float  color_rainbow  [256*3] ;
  float *colormaps[COLOR_NUMBER];

private :
  static void rainbow( float v, float *cols );

public:
  ColorRamp() ;
  ~ColorRamp() ;

  void reset_random() ;

  void get_color( unsigned char i, ColorMap cmap, float *cols, bool inv = false ) ;
  void set_GLcolor( double v, double vmax, ColorMap cmap, float alpha = 1.0, bool inv = false ) ;
};
//_____________________________________________________________________________


#endif // _COLORRAMP_H_
