// Source file for image class



// Include files 

#include "R2/R2.h"
#include "R2Pixel.h"
#include "R2Image.h"
#include "svd.h"
#include <vector>
#include <math.h>


////////////////////////////////////////////////////////////////////////
// Constructors/Destructors
////////////////////////////////////////////////////////////////////////


R2Image::
R2Image(void)
  : pixels(NULL),
    npixels(0),
    width(0), 
    height(0)
{
}



R2Image::
R2Image(const char *filename)
  : pixels(NULL),
    npixels(0),
    width(0), 
    height(0)
{
  // Read image
  Read(filename);
}



R2Image::
R2Image(int width, int height)
  : pixels(NULL),
    npixels(width * height),
    width(width), 
    height(height)
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);
}



R2Image::
R2Image(int width, int height, const R2Pixel *p)
  : pixels(NULL),
    npixels(width * height),
    width(width), 
    height(height)
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels 
  for (int i = 0; i < npixels; i++) 
    pixels[i] = p[i];
}



R2Image::
R2Image(const R2Image& image)
  : pixels(NULL),
    npixels(image.npixels),
    width(image.width), 
    height(image.height)
    
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels 
  for (int i = 0; i < npixels; i++) 
    pixels[i] = image.pixels[i];
}



R2Image::
~R2Image(void)
{
  // Free image pixels
  if (pixels) delete [] pixels;
}



R2Image& R2Image::
operator=(const R2Image& image)
{
  // Delete previous pixels
  if (pixels) { delete [] pixels; pixels = NULL; }

  // Reset width and height
  npixels = image.npixels;
  width = image.width;
  height = image.height;

  // Allocate new pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels 
  for (int i = 0; i < npixels; i++) 
    pixels[i] = image.pixels[i];

  // Return image
  return *this;
}


void R2Image::
svdTest(void)
{
	// fit a 2D conic to five points
	R2Point p1(1.2,3.5);
	R2Point p2(2.1,2.2);
	R2Point p3(0.2,1.6);
	R2Point p4(0.0,0.5);
	R2Point p5(-0.2,4.2);

	// build the 5x6 matrix of equations
	double** linEquations = dmatrix(1,5,1,6);

	linEquations[1][1] = p1[0]*p1[0];
	linEquations[1][2] = p1[0]*p1[1];
	linEquations[1][3] = p1[1]*p1[1];
	linEquations[1][4] = p1[0];
	linEquations[1][5] = p1[1];
	linEquations[1][6] = 1.0;

	linEquations[2][1] = p2[0]*p2[0];
	linEquations[2][2] = p2[0]*p2[1];
	linEquations[2][3] = p2[1]*p2[1];
	linEquations[2][4] = p2[0];
	linEquations[2][5] = p2[1];
	linEquations[2][6] = 1.0;

	linEquations[3][1] = p3[0]*p3[0];
	linEquations[3][2] = p3[0]*p3[1];
	linEquations[3][3] = p3[1]*p3[1];
	linEquations[3][4] = p3[0];
	linEquations[3][5] = p3[1];
	linEquations[3][6] = 1.0;
	
	linEquations[4][1] = p4[0]*p4[0];
	linEquations[4][2] = p4[0]*p4[1];
	linEquations[4][3] = p4[1]*p4[1];
	linEquations[4][4] = p4[0];
	linEquations[4][5] = p4[1];
	linEquations[4][6] = 1.0;

	linEquations[5][1] = p5[0]*p5[0];
	linEquations[5][2] = p5[0]*p5[1];
	linEquations[5][3] = p5[1]*p5[1];
	linEquations[5][4] = p5[0];
	linEquations[5][5] = p5[1];
	linEquations[5][6] = 1.0;

	printf("\n Fitting a conic to five points:\n");
	printf("Point #1: %f,%f\n",p1[0],p1[1]);
	printf("Point #2: %f,%f\n",p2[0],p2[1]);
	printf("Point #3: %f,%f\n",p3[0],p3[1]);
	printf("Point #4: %f,%f\n",p4[0],p4[1]);
	printf("Point #5: %f,%f\n",p5[0],p5[1]);

	// compute the SVD
	double singularValues[7]; // 1..6
	double** nullspaceMatrix = dmatrix(1,6,1,6);
	svdcmp(linEquations, 5, 6, singularValues, nullspaceMatrix);

	// get the result
	printf("\n Singular values: %f, %f, %f, %f, %f, %f\n",singularValues[1],singularValues[2],singularValues[3],singularValues[4],singularValues[5],singularValues[6]);

	// find the smallest singular value:
	int smallestIndex = 1;
	for(int i=2;i<7;i++) if(singularValues[i]<singularValues[smallestIndex]) smallestIndex=i;

	// solution is the nullspace of the matrix, which is the column in V corresponding to the smallest singular value (which should be 0)
	printf("Conic coefficients: %f, %f, %f, %f, %f, %f\n",nullspaceMatrix[1][smallestIndex],nullspaceMatrix[2][smallestIndex],nullspaceMatrix[3][smallestIndex],nullspaceMatrix[4][smallestIndex],nullspaceMatrix[5][smallestIndex],nullspaceMatrix[6][smallestIndex]);

	// make sure the solution is correct:
	printf("Equation #1 result: %f\n",	p1[0]*p1[0]*nullspaceMatrix[1][smallestIndex] + 
										p1[0]*p1[1]*nullspaceMatrix[2][smallestIndex] + 
										p1[1]*p1[1]*nullspaceMatrix[3][smallestIndex] + 
										p1[0]*nullspaceMatrix[4][smallestIndex] + 
										p1[1]*nullspaceMatrix[5][smallestIndex] + 
										nullspaceMatrix[6][smallestIndex]);

	printf("Equation #2 result: %f\n",	p2[0]*p2[0]*nullspaceMatrix[1][smallestIndex] + 
										p2[0]*p2[1]*nullspaceMatrix[2][smallestIndex] + 
										p2[1]*p2[1]*nullspaceMatrix[3][smallestIndex] + 
										p2[0]*nullspaceMatrix[4][smallestIndex] + 
										p2[1]*nullspaceMatrix[5][smallestIndex] + 
										nullspaceMatrix[6][smallestIndex]);

	printf("Equation #3 result: %f\n",	p3[0]*p3[0]*nullspaceMatrix[1][smallestIndex] + 
										p3[0]*p3[1]*nullspaceMatrix[2][smallestIndex] + 
										p3[1]*p3[1]*nullspaceMatrix[3][smallestIndex] + 
										p3[0]*nullspaceMatrix[4][smallestIndex] + 
										p3[1]*nullspaceMatrix[5][smallestIndex] + 
										nullspaceMatrix[6][smallestIndex]);

	printf("Equation #4 result: %f\n",	p4[0]*p4[0]*nullspaceMatrix[1][smallestIndex] + 
										p4[0]*p4[1]*nullspaceMatrix[2][smallestIndex] + 
										p4[1]*p4[1]*nullspaceMatrix[3][smallestIndex] + 
										p4[0]*nullspaceMatrix[4][smallestIndex] + 
										p4[1]*nullspaceMatrix[5][smallestIndex] + 
										nullspaceMatrix[6][smallestIndex]);

	printf("Equation #5 result: %f\n",	p5[0]*p5[0]*nullspaceMatrix[1][smallestIndex] + 
										p5[0]*p5[1]*nullspaceMatrix[2][smallestIndex] + 
										p5[1]*p5[1]*nullspaceMatrix[3][smallestIndex] + 
										p5[0]*nullspaceMatrix[4][smallestIndex] + 
										p5[1]*nullspaceMatrix[5][smallestIndex] + 
										nullspaceMatrix[6][smallestIndex]);

	R2Point test_point(0.34,-2.8);

	printf("A point off the conic: %f\n",	test_point[0]*test_point[0]*nullspaceMatrix[1][smallestIndex] + 
											test_point[0]*test_point[1]*nullspaceMatrix[2][smallestIndex] + 
											test_point[1]*test_point[1]*nullspaceMatrix[3][smallestIndex] + 
											test_point[0]*nullspaceMatrix[4][smallestIndex] + 
											test_point[1]*nullspaceMatrix[5][smallestIndex] + 
											nullspaceMatrix[6][smallestIndex]);

	return;	
}



////////////////////////////////////////////////////////////////////////
// Image processing functions
// YOU IMPLEMENT THE FUNCTIONS IN THIS SECTION
////////////////////////////////////////////////////////////////////////

// Per-pixel Operations ////////////////////////////////////////////////

void R2Image::
Brighten(double factor)
{
  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      Pixel(i,j) *= factor;
      Pixel(i,j).Clamp();
    }
  }
}

//De-normalize kernel, de-gray val, and don't clamp for harris filter

void R2Image::
SobelX(void)
{
    R2Image tempImage(width,height);
    double kernel [3][3];
//    kernel [0][0] = 1.0/8.0;
//    kernel [0][1] = 2.0/8.0;
//    kernel [0][2] = 1.0/8.0;
//    kernel [1][0] = 0.0;
//    kernel [1][1] = 0.0;
//    kernel [1][2] = 0.0;
//    kernel [2][0] = -1.0/8.0;
//    kernel [2][1] = -2.0/8.0;
//    kernel [2][2] = -1.0/8.0;
    kernel [0][0] = 1.0;
    kernel [0][1] = 2.0;
    kernel [0][2] = 1.0;
    kernel [1][0] = 0.0;
    kernel [1][1] = 0.0;
    kernel [1][2] = 0.0;
    kernel [2][0] = -1.0;
    kernel [2][1] = -2.0;
    kernel [2][2] = -1.0;
    
    
    for (int i = 1; i<width-1; i++)
    {
        for (int j = 1; j<height-1; j++)
        {
            R2Pixel val(0.0, 0.0, 0.0, 1.0);
            for (int lx = -1; lx<2; lx++)
            {
                for (int ly = -1; ly<2; ly++)
                {
                    val += Pixel(i+lx,j+ly) * kernel[lx+1][ly+1];
                }
            }
            tempImage.Pixel(i,j) = val;
            //tempImage.Pixel(i,j).Clamp();
        }
    }
    
    (*this) = tempImage;
}


void R2Image::
SobelY(void)
{
    R2Image tempImage(*this);
    double kernel [3][3];
//    kernel [0][0] = 1.0/8.0;
//    kernel [0][1] = 0.0;
//    kernel [0][2] = -1.0/8.0;
//    kernel [1][0] = 2.0/8.0;
//    kernel [1][1] = 0.0;
//    kernel [1][2] = -2.0/8.0;
//    kernel [2][0] = 1.0/8.0;
//    kernel [2][1] = 0.0;
//    kernel [2][2] = -1.0/8.0;
    kernel [0][0] = 1.0;
    kernel [0][1] = 0.0;
    kernel [0][2] = -1.0;
    kernel [1][0] = 2.0;
    kernel [1][1] = 0.0;
    kernel [1][2] = -2.0;
    kernel [2][0] = 1.0;
    kernel [2][1] = 0.0;
    kernel [2][2] = -1.0;
    
    for (int i = 1; i<width-1; i++)
    {
        for (int j = 1; j<height-1; j++)
        {
            R2Pixel val(0.0, 0.0, 0.0, 1.0);
            for (int lx = -1; lx<2; lx++)
            {
                for (int ly = -1; ly<2; ly++)
                {
                    val += Pixel(i+lx,j+ly) * kernel[lx+1][ly+1];
                }
            }
            tempImage.Pixel(i,j) = val;
            //tempImage.Pixel(i,j).Clamp();
        }
    }
    
    (*this) = tempImage;
    
}



void R2Image::
LoG(void)
{
  fprintf(stderr, "LoG() not implemented\n");
}



void R2Image::
ChangeSaturation(double factor)
{
   fprintf(stderr, "ChangeSaturation() not implemented\n");
}


// Linear filtering ////////////////////////////////////////////////


double Gaussian(int x, double sigma)
{
    int k = floor(sigma*3.0);
    return exp(-1*pow(x-k, 2)/(2.0*sigma*sigma));
}

void R2Image::
Blur(double sigma)
{
    //make kernel
    int n = floor((sigma * 6.0 + 1));
    double kernel [n];
    double total = 0;
   
    for (int i = 0; i<n; i++)
    {
        kernel [i] = Gaussian(i,sigma);
        total += kernel [i];
    }

    //normalize kernel
    for (int i = 0; i<n; i++)
    {
        kernel [i] = kernel[i]/total;
    }
    
    
    //change the image
    R2Image tempImage (*this);
    
    for (int j = n; j<height-n; j++)
    {
        for (int i = n; i<width-n; i++)
        {
            R2Pixel val(0.0, 0.0, 0.0, 1.0);
            for (int lx = (-n/2); lx<=(n/2); lx++)
            {
                val += Pixel(i+lx,j) * kernel[lx+(n/2)];
                
            }
            tempImage.Pixel(i,j) = val;
            //tempImage.Pixel(i,j).Clamp();
        }
    }
    
    
    
    for (int j = n; j<height-n; j++)
    {
        for (int i = n; i<width-n; i++)
        {
            R2Pixel val(0.0, 0.0, 0.0, 1.0);
            for (int ly = (-n/2); ly<=(n/2); ly++)
            {
                val += tempImage.Pixel(i,j+ly) * kernel[ly+(n/2)];
            }
            Pixel(i,j) = val;
           // Pixel(i,j).Clamp();
        }
    }
}

void R2Image::
SharpenAlt(double sigma)
{
    R2Image tempImage = (*this);
    tempImage.Blur(sigma);
    for (int j = 0; j<height; j++)
    {
        for (int i = 0; i<width; i++)
        {
            tempImage.Pixel(i,j) = Pixel(i,j)-tempImage.Pixel(i,j);
            //tempImage.Pixel(i,j).Clamp();
        }
    }
    
    for (int j = 0; j<height; j++)
    {
        for (int i = 0; i<width; i++)
        {
            Pixel(i,j) = Pixel(i,j)-tempImage.Pixel(i,j);//can multiply tempImage.Pixel to get diff effect
            //Pixel(i,j).Clamp(); commenting out makes it sharper
        }
    }
}

void R2Image::
Sharpen()
{
    R2Image tempImage(width,height);
    double kernel [3][3];
    kernel [0][0] = -1.0;
    kernel [0][1] = -1.0;
    kernel [0][2] = -1.0;
    kernel [1][0] = -1.0;
    kernel [1][1] = 9.0;
    kernel [1][2] = -1.0;
    kernel [2][0] = -1.0;
    kernel [2][1] = -1.0;
    kernel [2][2] = -1.0;
    
    for (int i = 1; i<width-1; i++)
    {
        for (int j = 1; j<height-1; j++)
        {
            R2Pixel val(0.0, 0.0, 0.0, 1.0);
            for (int lx = -1; lx<2; lx++)
            {
                for (int ly = -1; ly<2; ly++)
                {
                    val += Pixel(i+lx,j+ly) * kernel[lx+1][ly+1];
                }
            }
            tempImage.Pixel(i,j) = val;
            tempImage.Pixel(i,j).Clamp();
        }
    }
    
    (*this) = tempImage;
}

int xLocationsOriginalImage [150];
int yLocationsOriginalImage [150];

int distanceToBlackOut(int x, int d)
{
    //int i = d, take out declaration in loop
    for (int i = d; i>-1; i--)
    {
        if ((sqrt(x*x+i*i))<=d)
        {
            return i;
        }
    }
    //return i?????
}

void R2Image::
Harris(double sigma)
{
    R2Image sobelxSquared (width,height);
    R2Image sobelySquared (width,height);
    R2Image sobelProduct (width,height);
    R2Image sobelx (*this);
    R2Image sobely (*this);
    
    sobelx.SobelX();
    sobely.SobelY();
    
    printf("SOBELS DONE\n");
    
    for (int i = 0; i<width; i++)
    {
        for (int j = 0; j<height; j++)
        {
            sobelxSquared.Pixel(i,j) = sobelx.Pixel(i,j)*sobelx.Pixel(i,j);
            sobelySquared.Pixel(i,j) = sobely.Pixel(i,j)*sobely.Pixel(i,j);
            sobelProduct.Pixel(i,j) = sobelx.Pixel(i,j)*sobely.Pixel(i,j);
        }
    }
    
    sobelxSquared.Blur(sigma);
    sobelySquared.Blur(sigma);
    sobelProduct.Blur(sigma);
    
    printf("SOBELS MATRIX DONE\n");
    
    R2Image harrisPicture (width,height);
    
    
    for (int i = 0; i<width; i++)
    {
        for (int j = 0; j<height; j++)
        {
            R2Pixel val(0.0, 0.0, 0.0, 1.0);
            val = sobelxSquared.Pixel(i,j)*sobelySquared.Pixel(i,j)-sobelProduct.Pixel(i,j)*sobelProduct.Pixel(i,j)
            -0.04*((sobelxSquared.Pixel(i,j)+sobelySquared.Pixel(i,j))*(sobelxSquared.Pixel(i,j)+sobelySquared.Pixel(i,j)));
            harrisPicture.Pixel (i,j) = val;
            
//            R2Pixel val(0.0, 0.0, 0.0, 1.0);
//            R2Pixel pointfive(0.5, 0.5, 0.5, 1.0);
//            val = sobelxSquared.Pixel(i,j)*sobelySquared.Pixel(i,j)-sobelProduct.Pixel(i,j)*sobelProduct.Pixel(i,j)
//            -0.04*((sobelxSquared.Pixel(i,j)+sobelySquared.Pixel(i,j))*(sobelxSquared.Pixel(i,j)+sobelySquared.Pixel(i,j)))+pointfive;
//            harrisPicture.Pixel (i,j) = val;
//            harrisPicture.Pixel (i,j).Clamp();
        }
    }
    
    printf("HARRIS PICTURE MADE\n");
    
    R2Image original (width,height);
    original = (*this);
    (*this) = harrisPicture;
    
    int count = 0;
    R2Pixel blackPixel(0.0,0.0,0.0,1.0);
    int xlocations [150];
    int ylocations [150];
    
    int radius = 20;//radius of pixels blacked out around harris point
    
    while (count < 150)
    {
        printf("Harris feature %d\n", count);
        //finds the brightest feature left
        int xloc = 0;
        int yloc = 0;
        R2Pixel brightest(0.0,0.0,0.0,1.0);
        for (int i = 0; i<width; i++)
        {
            for (int j = 0; j<height; j++)
            {
                if (harrisPicture.Pixel(i,j).Luminance()>brightest.Luminance())
                {
                    brightest = Pixel(i,j);
                    xloc = i;
                    yloc = j;
                }
            }
        }
        //blacks out pixels of radius d
        for (int k = 0; k<radius+1; k++)
        {
            for (int r = 0; r<distanceToBlackOut(k, radius)+1; r++)
            {
                    
                harrisPicture.Pixel(xloc+k,yloc+r) = blackPixel;
                harrisPicture.Pixel(xloc-k,yloc+r) = blackPixel;
                harrisPicture.Pixel(xloc-k,yloc-r) = blackPixel;
                harrisPicture.Pixel(xloc+k,yloc-r) = blackPixel;
            }
        }
        
        //stores locations
        xlocations[count] = xloc;
        ylocations[count] = yloc;
        xLocationsOriginalImage[count] = xloc;
        yLocationsOriginalImage[count] = yloc;
        count += 1;
    }
    

    
//    printf("HARRIS LOCATIONS\n");
//    for (int i = 0; i<150; i++)
//        {
//            printf ("(%d, %d)  \n", xLocationsOriginalImage[i], yLocationsOriginalImage[i]);
//        }
//    
//    printf("HARRIS LOCATIONS\n");
//    for (int i = 0; i<150; i++)
//    {
//        printf ("(%d, %d)  \n", xlocations[i], ylocations[i]);
//    }
    
    R2Image boxes (width,height);
    boxes = original;
    (*this) = original;
    
    //draws box;
        R2Pixel redPixel (1.0, 0.0, 0.0, 1.0);
        for (int i = 0; i<150; i++)
        {
            int x = xlocations[i];
            int y = ylocations[i];
            for (int z = -10; z<10; z++)
            {
                boxes.Pixel(x+z,y+10)=redPixel;
                boxes.Pixel(x+z,y-10)=redPixel;
                boxes.Pixel(x+10,y+z)=redPixel;
                boxes.Pixel(x-10,y+z)=redPixel;
            }
            boxes.Pixel(x,y) = redPixel;
        }
    
    
    //USE THIS TO SEE HARRIS IMAGE
    //(*this) = harrisPicture;
    
    //USE THISE TO SEE POINTS FOUND
    //(*this) = boxes;
}



void R2Image::
line(int x0, int x1, int y0, int y1, float r, float g, float b)
{
    if(x0>x1)
    {
        int x=y1;
        y1=y0;
        y0=x;
        
        x=x1;
        x1=x0;
        x0=x;
    }
    int deltax = x1 - x0;
    int deltay = y1 - y0;
    float error = 0;
    float deltaerr = 0.0;
    if(deltax!=0) deltaerr =fabs(float(float(deltay) / deltax));    // Assume deltax != 0 (line is not vertical),
    // note that this division needs to be done in a way that preserves the fractional part
    int y = y0;
    for(int x=x0;x<=x1;x++)
    {
        Pixel(x,y).Reset(r,g,b,1.0);
        error = error + deltaerr;
        if(error>=0.5)
        {
            if(deltay>0) y = y + 1;
            else y = y - 1;
            
            error = error - 1.0;
        }
    }
    if(x0>3 && x0<width-3 && y0>3 && y0<height-3)
    {
        for(int x=x0-3;x<=x0+3;x++)
        {
            for(int y=y0-3;y<=y0+3;y++)
            {
                Pixel(x,y).Reset(r,g,b,1.0);
            }
        }
    }
}

double SSD(double a, double b)
{
    double diff = a-b;
    double squared = diff*diff;
    return squared;
}


void R2Image::
blendOtherImageTranslated(R2Image * otherImage)
{
    R2Image imageAoriginal = (*this);
    R2Image imageBoriginal = (*otherImage);
    R2Image harrisPictureA = (*this);
    harrisPictureA.Harris(2);
    
    float windowWidth = width/5;
    float windowHeight = height/5;
    int location;
    double mindifference;
    double sumsquareddifference;
    int templateSize = 13;
    
    double imageATemplate [templateSize*templateSize];
    
    //the final arrays to store locations of imageB and temporary indices to look at different pixels of image B
    int finalxlocationsImageB [150];
    int finalylocationsImageB [150];
    int tempxlocationImageB;
    int tempylocationImageB;
    //loop through each location found in harris filter and work with each one to find a match
    for (int i = 0; i<150; i++)
    {
        //create image template for pixel in imageA
        int xlocationimageA = xLocationsOriginalImage[i];
        int ylocationimageA = yLocationsOriginalImage[i];
        location = 0;
        for (int n = -templateSize/2; n <(templateSize/2)+1; n++)
        {
            for (int m = -6; m<7; m++)
            {
                imageATemplate[location] = imageAoriginal.Pixel(xlocationimageA+n,ylocationimageA+m).Luminance();
                location++;
            }
        }
        //image template for location is now stored. Next, look at SSD in imageB
        mindifference = 1000000000;
        
        //looping through imageB in a window 20% size of original image
        for (int n = -windowWidth/2; n<windowWidth/2; n++)
        {
            for (int m = -windowHeight/2; m<windowHeight/2; m++)
            {
                //store each location I am about to make an image template for in imageB
                tempxlocationImageB = xlocationimageA+n;
                tempylocationImageB = ylocationimageA+m;
                //create array for image template to be created in imageB
                double imageBTemplate [templateSize*templateSize];
                location = 0;
                //make image template for each pixel in 13x13 window
                for (int j = -templateSize/2; j <(templateSize/2)+1; j++)
                {
                    for (int k = -templateSize/2; k<(templateSize/2)+1; k++)
                    {
                        imageBTemplate[location] = imageBoriginal.Pixel(tempxlocationImageB+j,tempylocationImageB+k).Luminance();
                        location++;
                    }
                }
                sumsquareddifference = 0;
                //image template for imageB is now created
                //compute sum of squared differences between imageA and imageB template
                for (int s = 0; s<templateSize*templateSize; s++)
                {
                    sumsquareddifference+=SSD(imageATemplate[s],imageBTemplate[s]);
                }
                //if I find a smaller SSD, store the coordinates of the pixel yielding that smaller SSD in finalimageBlocations array
                if (sumsquareddifference<mindifference)
                {
                    mindifference = sumsquareddifference;
                    finalxlocationsImageB[i] = tempxlocationImageB;
                    finalylocationsImageB[i] = tempylocationImageB;
                }
            }
        }
    }

    // RANSAC - 150 points have been found //////////////////////////////////////////////////////////////////////////////
    
    //Variables to use
    float threshold = 5.00;
    int indices [150]; //array to keep track of which indices have already been checked
    int inliers = 0; //final number of inliers when algorithm terminates
    int index = 0; //the index that I am checking with each iteration
    int finalindex = 0; //the index that gives the best RANSAC match
    
    //set up indices to randomly check. If the index of the array is -1, that index has already been checked
    for (int i = 0; i<150; i++)
    {
        indices[i]=i;
    }
    
    
    //The RANSAC part. Choose 50
    for (int q = 0; q<50; q++)
    {
        //ensure I get a valid index that has not already been checked
        while (indices[index]<0)
        {
            index = rand() % 150;
        }
        
        //sets up part of distance formula. changes vector to origin and finds how much it is displaced by the translation
        int xpicked = finalxlocationsImageB[index]-xLocationsOriginalImage[index];
        int ypicked = finalylocationsImageB[index]-yLocationsOriginalImage[index];
        
        //temporary variable to keep track of number of inliers for each random vector
        int tempinliers = 0;
        
        //I have the vector. Now run through all other correspondences and compare them to picked vector
        for (int i = 0; i<150; i++)
        {
            if (i!=index)
            {
                //compute distance between every vector and picked vectore
                int xcompare = finalxlocationsImageB[i]-xLocationsOriginalImage[i];
                int ycompare = finalylocationsImageB[i]-yLocationsOriginalImage[i];
                
                int xsquared = (xcompare-xpicked)*(xcompare-xpicked);
                int ysquared = (ycompare-ypicked)*(ycompare-ypicked);
                
                float distance = sqrt(xsquared+ysquared);
                
                if (distance<=threshold)
                {
                    tempinliers++;
                }
            }
        }
        
        //check which vector gives best set of inliers
        if (tempinliers>inliers)
        {
            finalindex = index;
            inliers = tempinliers;
        }
        
        indices[index]=-1;
    }
    
    //I have the index that gives the best vector
    
    //vectors to store locations of good vectors and bad vectors
    std::vector<int> xLocationsARed;
    std::vector<int> yLocationsARed;
    std::vector<int> xLocationsBRed;
    std::vector<int> yLocationsBRed;
    
    std::vector<int> xLocationsABlue;
    std::vector<int> yLocationsABlue;
    std::vector<int> xLocationsBBlue;
    std::vector<int> yLocationsBBlue;
    
    //loop that sets up which vectors will be blue and which will be red
    for (int i = 0; i<150; i++)
    {
        int xpicked = finalxlocationsImageB[finalindex]-xLocationsOriginalImage[finalindex];
        int ypicked = finalylocationsImageB[finalindex]-yLocationsOriginalImage[finalindex];
        int xcompare = finalxlocationsImageB[i]-xLocationsOriginalImage[i];
        int ycompare = finalylocationsImageB[i]-yLocationsOriginalImage[i];
        
        int xsquared = (xcompare-xpicked)*(xcompare-xpicked);
        int ysquared = (ycompare-ypicked)*(ycompare-ypicked);
    
        float distance = sqrt(xsquared+ysquared);
        
        if (distance<=threshold)
        {
            xLocationsARed.push_back(xLocationsOriginalImage[i]);
            yLocationsARed.push_back(yLocationsOriginalImage[i]);
            xLocationsBRed.push_back(finalxlocationsImageB[i]);
            yLocationsBRed.push_back(finalylocationsImageB[i]);
        }
        
        if (distance>threshold)
        {
            xLocationsABlue.push_back(xLocationsOriginalImage[i]);
            yLocationsABlue.push_back(yLocationsOriginalImage[i]);
            xLocationsBBlue.push_back(finalxlocationsImageB[i]);
            yLocationsBBlue.push_back(finalylocationsImageB[i]);
        }
    }
    
    
    //draw the good vectors in red
    for (int i = 0; i<xLocationsARed.size(); i++)
    {
        if (xLocationsBRed[i]<0)
            xLocationsBRed[i]=0;
        if (yLocationsBRed[i]<0)
            yLocationsBRed[i]=0;
        line(xLocationsARed[i],xLocationsBRed[i],
                yLocationsARed[i],yLocationsBRed[i],
            1.0,0.0,0.0);
    }
    
    
    //draw the bad vectors in blue
    for (int i = 0; i<xLocationsABlue.size(); i++)
    {
        if (xLocationsBBlue[i]<0)
            xLocationsBBlue[i]=0;
        if (yLocationsBBlue[i]<0)
            yLocationsBBlue[i]=0;
        line(xLocationsABlue[i],xLocationsBBlue[i],
             yLocationsABlue[i],yLocationsBBlue[i],
             0.0,0.0,1.0);
    }
    
}

double** HMatrix(std::vector<R2Point> pointsA, std::vector<R2Point> pointsB)
{
    std::vector<R2Point> aPoints;
    std::vector<R2Point> bPoints;
    
    aPoints = pointsA;
    bPoints = pointsB;
    
    double** A = dmatrix(1,aPoints.size()*2, 1, 9);
    
    for (int i = 1; i<=aPoints.size(); i++)
    {
        A[i*2-1][1] = 0;
        A[i*2][1] = aPoints[i-1].X();
        
        A[i*2-1][2] = 0;
        A[i*2][2] = aPoints[i-1].Y();
        
        A[i*2-1][3] = 0;
        A[i*2][3] = 1;
        
        A[i*2-1][4] = -1.0*aPoints[i-1].X();
        A[i*2][4] = 0;
        
        A[i*2-1][5] = -1.0*aPoints[i-1].Y();
        A[i*2][5] = 0;
        
        A[i*2-1][6] = -1;
        A[i*2][6] = 0;
        
        A[i*2-1][7] = bPoints[i-1].Y()*aPoints[i-1].X();
        A[i*2][7] = -1*bPoints[i-1].X()*aPoints[i-1].X();
        
        A[i*2-1][8] = bPoints[i-1].Y()*aPoints[i-1].Y();
        A[i*2][8] = -1*bPoints[i-1].X()*aPoints[i-1].Y();
        
        A[i*2-1][9] = bPoints[i-1].Y();
        A[i*2][9] = -1*bPoints[i-1].X();
        
    }
    
    double w[10];
    double** nullSpaceMatrix = dmatrix(1,9,1,9);
    svdcmp(A, aPoints.size()*2, 9, w, nullSpaceMatrix);
    
    double smallestSingularValue = w[1];
    int smallestIndex = 1;
    
    for(int i=2;i<=9;i++)
    {
        if(w[i]<w[smallestIndex])
        {
            smallestIndex=i;
            smallestSingularValue = w[smallestIndex];
        }
    }
    
    
    double** H = dmatrix(1,3, 1, 3);
    
    H[1][1] = nullSpaceMatrix[1][smallestIndex];
    H[1][2] = nullSpaceMatrix[2][smallestIndex];
    H[1][3] = nullSpaceMatrix[3][smallestIndex];
    H[2][1] = nullSpaceMatrix[4][smallestIndex];
    H[2][2] = nullSpaceMatrix[5][smallestIndex];
    H[2][3] = nullSpaceMatrix[6][smallestIndex];
    H[3][1] = nullSpaceMatrix[7][smallestIndex];
    H[3][2] = nullSpaceMatrix[8][smallestIndex];
    H[3][3] = nullSpaceMatrix[9][smallestIndex];
    
    return H;
}

float distance(R2Point a, R2Point b)
{

    float xpart = (a.X()-b.X())*(a.X()-b.X());
    float ypart = (a.Y()-b.Y())*(a.Y()-b.Y());
    float s = sqrt(xpart+ypart);
    return s;
}

double** HMatrixNormal(std::vector<R2Point> pointsA, std::vector<R2Point> pointsB)
{
    if (pointsA.size()!=pointsB.size())
    {
        printf("NOT VALID");
    }
    
    std::vector<R2Point> aPointsNotNormal;
    std::vector<R2Point> bPointsNotNormal;
    
    aPointsNotNormal = pointsA;
    bPointsNotNormal = pointsB;
    
    //aPoints and bPoints are the vectors I am using.
    float txA = 0;
    float tyA = 0;
    float txB = 0;
    float tyB = 0;
    
    for (int i = 0; i<aPointsNotNormal.size(); i++)
    {
        txA = txA + aPointsNotNormal[i].X();
        tyA = tyA + aPointsNotNormal[i].Y();
        
        txB = txB + bPointsNotNormal[i].X();
        tyB = tyB + bPointsNotNormal[i].Y();
    }
    
    txA = txA/aPointsNotNormal.size();
    tyA = tyA/aPointsNotNormal.size();
    txB = txB/aPointsNotNormal.size();
    tyB = tyB/aPointsNotNormal.size();
    
    R2Point centroidA (txA, tyA);
    R2Point centroidB (txB, tyB);
    
    float dA = 0;
    float dB = 0;
    
    for (int i = 0; i<aPointsNotNormal.size(); i++)
    {
        dA = dA + distance(aPointsNotNormal[i],centroidA);
        
        dB = dB + distance(bPointsNotNormal[i],centroidB);
    }
    
    dA = dA/aPointsNotNormal.size();
    dB = dB/bPointsNotNormal.size();
    
    float sA = 1.41421356237/dA;
    float sB = 1.41421356237/dB;
    
    std::vector<R2Point> aPoints;
    std::vector<R2Point> bPoints;
    
    //have txA tyA and sA  and txB tyB and sB
    
    for (int i = 0; i<aPointsNotNormal.size(); i++)
    {
        float xA = aPointsNotNormal[i].X();
        float yA = aPointsNotNormal[i].Y();
        xA = (xA-txA)*sA;
        yA = (yA-tyA)*sA;
        R2Point pA (xA,yA);
        
        float xB = bPointsNotNormal[i].X();
        float yB = bPointsNotNormal[i].Y();
        
        xB = (xB-txB)*sB;
        yB = (yB-tyB)*sB;
        
        R2Point pB (xB,yB);
        
        aPoints.push_back(pA);
        bPoints.push_back(pB);
    }
    
    double** A = dmatrix(1,aPoints.size()*2, 1, 9);
    
    for (int i = 1; i<=aPoints.size(); i++)
    {
        A[i*2-1][1] = 0;
        A[i*2][1] = aPoints[i-1].X();
        
        A[i*2-1][2] = 0;
        A[i*2][2] = aPoints[i-1].Y();
        
        A[i*2-1][3] = 0;
        A[i*2][3] = 1;
        
        A[i*2-1][4] = -1.0*aPoints[i-1].X();
        A[i*2][4] = 0;
        
        A[i*2-1][5] = -1.0*aPoints[i-1].Y();
        A[i*2][5] = 0;
        
        A[i*2-1][6] = -1;
        A[i*2][6] = 0;
        
        A[i*2-1][7] = bPoints[i-1].Y()*aPoints[i-1].X();
        A[i*2][7] = -1*bPoints[i-1].X()*aPoints[i-1].X();
        
        A[i*2-1][8] = bPoints[i-1].Y()*aPoints[i-1].Y();
        A[i*2][8] = -1*bPoints[i-1].X()*aPoints[i-1].Y();
        
        A[i*2-1][9] = bPoints[i-1].Y();
        A[i*2][9] = -1*bPoints[i-1].X();
    }
    
    double w[10];
    double** nullSpaceMatrix = dmatrix(1,9,1,9);
    svdcmp(A, aPoints.size()*2, 9, w, nullSpaceMatrix);
    
    double smallestSingularValue = w[1];
    int smallestIndex = 1;
    
    for(int i=2;i<=9;i++)
    {
        if(w[i]<w[smallestIndex])
        {
            smallestIndex=i;
            smallestSingularValue = w[smallestIndex];
        }
    }
    
    double** Hnorm = dmatrix(1,3, 1, 3);
    
    Hnorm[1][1] = nullSpaceMatrix[1][smallestIndex];
    Hnorm[1][2] = nullSpaceMatrix[2][smallestIndex];
    Hnorm[1][3] = nullSpaceMatrix[3][smallestIndex];
    Hnorm[2][1] = nullSpaceMatrix[4][smallestIndex];
    Hnorm[2][2] = nullSpaceMatrix[5][smallestIndex];
    Hnorm[2][3] = nullSpaceMatrix[6][smallestIndex];
    Hnorm[3][1] = nullSpaceMatrix[7][smallestIndex];
    Hnorm[3][2] = nullSpaceMatrix[8][smallestIndex];
    Hnorm[3][3] = nullSpaceMatrix[9][smallestIndex];
    
    //have txA tyA and sA  and txB tyB and sB
    double** T = dmatrix(1,3,1,3);
    double** Tprime = dmatrix(1,3,1,3);
    double** TprimeInv = dmatrix(1,3,1,3);
    
    T[1][1] = sA;
    T[1][2] = 0;
    T[1][3] = sA*txA*-1;
    T[2][1] = 0;
    T[2][2] = sA;
    T[2][3] = sA*tyA*-1;
    T[3][1] = 0;
    T[3][2] = 0;
    T[3][3] = 1;
    
    Tprime[1][1] = sB;
    Tprime[1][2] = 0;
    Tprime[1][3] = sB*txB*-1;
    Tprime[2][1] = 0;
    Tprime[2][2] = sB;
    Tprime[2][3] = sB*tyB*-1;
    Tprime[3][1] = 0;
    Tprime[3][2] = 0;
    Tprime[3][3] = 1;
    
    double determinant = Tprime[1][1]*(Tprime[2][2]*Tprime[3][3]-Tprime[3][2]*Tprime[2][3])
                        -Tprime[1][2]*(Tprime[2][1]*Tprime[3][3]-Tprime[2][3]*Tprime[3][1])
                        +Tprime[1][3]*(Tprime[2][1]*Tprime[3][2]-Tprime[2][2]*Tprime[3][1]);
    double invdet = 1/determinant;
    
    TprimeInv[1][1] = (Tprime[2][2]*Tprime[3][3]-Tprime[3][2]*Tprime[2][3])*invdet;
    TprimeInv[1][2] = -(Tprime[1][2]*Tprime[3][3]-Tprime[1][3]*Tprime[3][2])*invdet;
    TprimeInv[1][3] = (Tprime[1][2]*Tprime[2][3]-Tprime[1][3]*Tprime[2][2])*invdet;
    
    TprimeInv[2][1] = -(Tprime[2][1]*Tprime[3][3]-Tprime[2][3]*Tprime[3][1])*invdet;
    TprimeInv[2][2] = (Tprime[1][1]*Tprime[3][3]-Tprime[1][3]*Tprime[3][1])*invdet;
    TprimeInv[2][3] = -(Tprime[1][1]*Tprime[2][3]-Tprime[2][1]*Tprime[1][3])*invdet;
    
    TprimeInv[3][1] = (Tprime[2][1]*Tprime[3][2]-Tprime[3][1]*Tprime[2][2])*invdet;
    TprimeInv[3][2] = -(Tprime[1][1]*Tprime[3][2]-Tprime[3][1]*Tprime[1][2])*invdet;
    TprimeInv[3][3] = (Tprime[1][1]*Tprime[2][2]-Tprime[2][1]*Tprime[1][2])*invdet;
    
    double** Hhalf = dmatrix(1,3, 1, 3);
    double** H = dmatrix(1,3, 1, 3);
    
    H[1][1] = 0;
    H[1][2] = 0;
    H[1][3] = 0;
    H[2][1] = 0;
    H[2][2] = 0;
    H[2][3] = 0;
    H[3][1] = 0;
    H[3][2] = 0;
    H[3][3] = 0;
    
    Hhalf[1][1] = 0;
    Hhalf[1][2] = 0;
    Hhalf[1][3] = 0;
    Hhalf[2][1] = 0;
    Hhalf[2][2] = 0;
    Hhalf[2][3] = 0;
    Hhalf[3][1] = 0;
    Hhalf[3][2] = 0;
    Hhalf[3][3] = 0;

    for(int i = 1; i < 4; i++)
    {
        for(int j = 1; j < 4; j++)
        {
            for(int k = 1; k < 4; k++)
            {
                Hhalf[i][j] +=  Hnorm[i][k] *  T[k][j];
            }
        }
    }
    
    for(int i = 1; i < 4; i++)
    {
        for(int j = 1; j < 4; j++)
        {
            for(int k = 1; k < 4; k++)
            {
                H[i][j] +=  TprimeInv[i][k] *  Hhalf[k][j];
            }
        }
    }
    
    return H;
}




//THINGS THAT ALTER THE ALGORITHM ARE DISTANCE OF BLACK OUT PIXELS IN HARRIS AND THRESHOLD OF DISTANCE FOR INLIERS

void R2Image::
blendOtherImageHomography(R2Image * otherImage)
{
//	// find at least 100 features on this image, and another 100 on the "otherImage". Based on these,
//	// compute the matching homography, and blend the transformed "otherImage" into this image with a 50% opacity.
    
    R2Image imageAoriginal = (*this);
    R2Image imageBoriginal = (*otherImage);
    R2Image harrisPictureA = (*this);
    harrisPictureA.Harris(2);
    
    printf("Harris filter done\n");
    
    float windowWidth = width/5;
    float windowHeight = height/5;
    
    int templateWidth = 13;//how big the image template is
    
    int location;
    double mindifference;
    double sumsquareddifference;
    double imageATemplate [templateWidth*templateWidth];
    
    //the final arrays to store locations of imageB and temporary indices to look at different pixels of image B
    int finalxlocationsImageB [150];
    int finalylocationsImageB [150];
    int tempxlocationImageB;
    int tempylocationImageB;
    //loop through each location found in harris filter and work with each one to find a match
    for (int i = 0; i<150; i++)
    {
        printf("Working on feature %d\n", i);
        //create image template for pixel in imageA
        int xlocationimageA = xLocationsOriginalImage[i];
        int ylocationimageA = yLocationsOriginalImage[i];
        location = 0;
        for (int n = -templateWidth/2; n <(templateWidth/2)+1; n++)
        {
            for (int m = -templateWidth/2; m<(templateWidth/2)+1; m++)
            {
                imageATemplate[location] = imageAoriginal.Pixel(xlocationimageA+n,ylocationimageA+m).Luminance();
                location++;
            }
        }
        //image template for location is now stored. Next, look at SSD in imageB
        mindifference = 1000000000;
        
        //looping through imageB in a window 20% size of original image
        for (int n = -windowWidth/2; n<windowWidth/2; n++)
        {
            for (int m = -windowHeight/2; m<windowHeight/2; m++)
            {
                //store each location I am about to make an image template for in imageB
                tempxlocationImageB = xlocationimageA+n;
                tempylocationImageB = ylocationimageA+m;
                //create array for image template to be created in imageB
                double imageBTemplate [templateWidth*templateWidth];
                location = 0;
                //make image template for each pixel in 13x13 window
                for (int j = -templateWidth/2; j <(templateWidth/2)+1; j++)
                {
                    for (int k = -templateWidth/2; k<(templateWidth/2)+1; k++)
                    {
                        imageBTemplate[location] = imageBoriginal.Pixel(tempxlocationImageB+j,tempylocationImageB+k).Luminance();
                        location++;
                    }
                }
                sumsquareddifference = 0;
                //image template for imageB is now created
                //compute sum of squared differences between imageA and imageB template
                for (int s = 0; s<templateWidth*templateWidth; s++)
                {
                    sumsquareddifference+=SSD(imageATemplate[s],imageBTemplate[s]);
                }
                //if I find a smaller SSD, store the coordinates of the pixel yielding that smaller SSD in finalimageBlocations array
                if (sumsquareddifference<mindifference)
                {
                    mindifference = sumsquareddifference;
                    finalxlocationsImageB[i] = tempxlocationImageB;
                    finalylocationsImageB[i] = tempylocationImageB;
                }
            }
        }
    }
    
    printf("DONE WITH FEATURE PART\n");
    
    // RANSAC - 150 points have been found ///////////////////////////////////////////////////////////////////////
    
    //Variables to use
    float threshold = 2.00;
    int inliers = 0; //final number of inliers when algorithm terminates
    std::vector <int> allInlierIndices;//array of all the indices when RANSAC terminates
    int fourIndicesForRansac [4]; //the array of indices that I am checking with each iteration
    std::vector <R2Point> allInlierPoints;//vector of all points that are inliers
    
    std::vector <R2Point> RANSACPointsA;
    std::vector <R2Point> RANSACPointsB;
    
    //The RANSAC part. Choose 300
    for (int q = 0; q<300; q++)
    {
        //find 4 indices and make sure none are the same
        fourIndicesForRansac [0] = rand() % 150;
        fourIndicesForRansac [1] = rand() % 150;
        fourIndicesForRansac [2] = rand() % 150;
        fourIndicesForRansac [3] = rand() % 150;
        
        while (fourIndicesForRansac[0] == fourIndicesForRansac[1] || fourIndicesForRansac[0] == fourIndicesForRansac[2] ||
               fourIndicesForRansac[0] == fourIndicesForRansac[3] ||fourIndicesForRansac[1] == fourIndicesForRansac[2] ||
               fourIndicesForRansac[1] == fourIndicesForRansac[3] ||fourIndicesForRansac[2] == fourIndicesForRansac[3])
        {
            fourIndicesForRansac [0] = rand() % 150;
            fourIndicesForRansac [1] = rand() % 150;
            fourIndicesForRansac [2] = rand() % 150;
            fourIndicesForRansac [3] = rand() % 150;
        }
        
        //make points out of those indices
        RANSACPointsA.clear();
        RANSACPointsB.clear();
        
        for (int i = 0; i<4; i++)
        {
            R2Point A(xLocationsOriginalImage[fourIndicesForRansac[i]],yLocationsOriginalImage[fourIndicesForRansac[i]]);
            R2Point B(finalxlocationsImageB[fourIndicesForRansac[i]],finalylocationsImageB[fourIndicesForRansac[i]]);
            RANSACPointsA.push_back(A);
            RANSACPointsB.push_back(B);
        }
        
        //Vectors have been set up now. I have the 4 locations stored in vectors.
        //now make an H matrix from those 4 points
        double** H = HMatrixNormal(RANSACPointsA,RANSACPointsB);
        
        //temporary variables to compare to final values
        int tempinliers = 0;
        std:: vector <int> tempIndicies;
        tempIndicies.clear();
        
        for (int i = 0; i<150; i++)
        {
            //matrix multiplication
            double H11 = H[1][1];
            double H12 = H[1][2];
            double H13 = H[1][3];
            double H21 = H[2][1];
            double H22 = H[2][2];
            double H23 = H[2][3];
            double H31 = H[3][1];
            double H32 = H[3][2];
            double H33 = H[3][3];
            
            float xA = xLocationsOriginalImage[i];
            float yA = yLocationsOriginalImage[i];
            float zA = 1;
            float xB = finalxlocationsImageB[i];
            float yB = finalylocationsImageB[i];
            
            float xPrime = xA*H11 + yA*H12 + zA*H13;
            float yPrime = xA*H21 + yA*H22 + zA*H23;
            float zPrime = xA*H31 + yA*H32 + zA*H33;
            
            //need to translate the points to z=1 plane. HOMOGENEUOS COORDINATES
            xPrime = xPrime/zPrime;
            yPrime = yPrime/zPrime;
            zPrime = zPrime/zPrime;
            
            float xsquared = (xPrime-xB)*(xPrime-xB);
            float ysquared = (yPrime-yB)*(yPrime-yB);
            
            float distance = sqrt(xsquared+ysquared);
            
            if (distance<=threshold)
            {
                tempinliers++;
                tempIndicies.push_back(i);
            }
        }
        if (tempinliers>inliers)
        {
            inliers = tempinliers;
            allInlierIndices = tempIndicies;
        }
    }
    
    std::vector<int> xLocationsARed;
    std::vector<int> yLocationsARed;
    std::vector<int> xLocationsBRed;
    std::vector<int> yLocationsBRed;
    
    std::vector<int> xLocationsABlue;
    std::vector<int> yLocationsABlue;
    std::vector<int> xLocationsBBlue;
    std::vector<int> yLocationsBBlue;
    
    std::vector<int> allIndicesAfterDrawing;
    
    for (int i = 0; i<150; i++)
    {
        if (i == allInlierIndices[0])
        {
            xLocationsARed.push_back(xLocationsOriginalImage[i]);
            yLocationsARed.push_back(yLocationsOriginalImage[i]);
            xLocationsBRed.push_back(finalxlocationsImageB[i]);
            yLocationsBRed.push_back(finalylocationsImageB[i]);
            allInlierIndices.erase(allInlierIndices.begin());
            allIndicesAfterDrawing.push_back(i);
        }
        else
        {
            xLocationsABlue.push_back(xLocationsOriginalImage[i]);
            yLocationsABlue.push_back(yLocationsOriginalImage[i]);
            xLocationsBBlue.push_back(finalxlocationsImageB[i]);
            yLocationsBBlue.push_back(finalylocationsImageB[i]);
        }
    }
    
    R2Image originalImage = (*this);
    
    //draw the good vectors in red
    for (int i = 0; i<xLocationsARed.size(); i++)
    {
        if (xLocationsBRed[i]<0)
            xLocationsBRed[i]=0;
        if (yLocationsBRed[i]<0)
            yLocationsBRed[i]=0;
        line(xLocationsARed[i],xLocationsBRed[i],
             yLocationsARed[i],yLocationsBRed[i],
             1.0,0.0,0.0);
    }
    
    
    //draw the bad vectors in blue
    for (int i = 0; i<xLocationsABlue.size(); i++)
    {
        if (xLocationsBBlue[i]<0)
            xLocationsBBlue[i]=0;
        if (yLocationsBBlue[i]<0)
            yLocationsBBlue[i]=0;
        line(xLocationsABlue[i],xLocationsBBlue[i],
             yLocationsABlue[i],yLocationsBBlue[i],
             0.0,0.0,1.0);
    }
    
    //make the final H matrix
    std::vector <R2Point> FinalPointsA;
    std::vector <R2Point> FinalPointsB;
    
    for (int i = 0; i<xLocationsARed.size(); i++)
    {
        R2Point A(xLocationsOriginalImage[allIndicesAfterDrawing[i]],yLocationsOriginalImage[allIndicesAfterDrawing[i]]);
        R2Point B(finalxlocationsImageB[allIndicesAfterDrawing[i]],finalylocationsImageB[allIndicesAfterDrawing[i]]);
        FinalPointsA.push_back(A);
        FinalPointsB.push_back(B);
    }

    double** H = HMatrixNormal(FinalPointsA,FinalPointsB);
    
    printf("\n\nthe final H matrix is:\n");
    
    for (int i = 1; i<4; i++)
    {
        for (int j = 1; j<4; j++)
        {
            printf("%f   ", H[i][j]);
        }
        printf("\n");
    }
    
    printf("\n");
    
    //RANSAC DONE. NOW BLEND THE TWO IMAGES TOGETHER////////////////////////////////////////////
    
    double H11 = H[1][1];
    double H12 = H[1][2];
    double H13 = H[1][3];
    double H21 = H[2][1];
    double H22 = H[2][2];
    double H23 = H[2][3];
    double H31 = H[3][1];
    double H32 = H[3][2];
    double H33 = H[3][3];
    
    R2Image Correspondences = (*this);
    
    R2Image finalImage(width, height);
    R2Pixel blackPixel(0.0,0.0,0.0,1.0);
    
    //Bi-linear interpolation
    for (int x = 0; x<width; x++)
    {
        for (int y = 0; y<height; y++)
        {
            float tempxLocB = H11 * x + H12 * y + H13;
            float tempyLocB = H21 * x + H22 * y + H23;
            float tempzLocB = H31 * x + H32 * y + H33;
            
            float xLocBfloat = tempxLocB/tempzLocB;
            float yLocBfloat = tempyLocB/tempzLocB;
            
            int xLocB = floor(xLocBfloat);
            int yLocB = floor(yLocBfloat);
            
            float d = yLocBfloat - yLocB;
            float dprime = xLocBfloat - xLocB;
            
            if (xLocB>=0 && xLocB<width && yLocB>=0 && yLocB<height)
            {
                finalImage.Pixel(x,y) = imageBoriginal.Pixel(xLocB, yLocB) * (1-d) * (1-dprime)
                                      + imageBoriginal.Pixel(xLocB, yLocB+1) * d * (1-dprime)
                                      + imageBoriginal.Pixel(xLocB+1, yLocB) * (1-d) * dprime
                                      + imageBoriginal.Pixel(xLocB+1, yLocB+1) * d * dprime;
            }
            else
            {
                finalImage.Pixel(x,y) = (originalImage.Pixel(x,y)+blackPixel)*0.5;
                //finalImage.Pixel(x,y) = originalImage.Pixel(x,y);
            }
            
        }
    }
    
    (*this) = Correspondences;
    (*this) = finalImage;
    
    printf("\ninliers is %d\n", inliers);
    
    //Algorithm correctly warps image B onto image A and can blend them together
    
    //How to edit parameters to control sigma and distance to black out in Harris filter, and threshold in RANSAC
}

////////////////////////////////////////////////////////////////////////
// I/O Functions
////////////////////////////////////////////////////////////////////////

int R2Image::
Read(const char *filename)
{
  // Initialize everything
  if (pixels) { delete [] pixels; pixels = NULL; }
  npixels = width = height = 0;

  // Parse input filename extension
  char *input_extension;
  if (!(input_extension = (char*)strrchr(filename, '.'))) {
    fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
    return 0;
  }
  
  // Read file of appropriate type
  if (!strncmp(input_extension, ".bmp", 4)) return ReadBMP(filename);
  else if (!strncmp(input_extension, ".ppm", 4)) return ReadPPM(filename);
  else if (!strncmp(input_extension, ".jpg", 4)) return ReadJPEG(filename);
  else if (!strncmp(input_extension, ".jpeg", 5)) return ReadJPEG(filename);
  
  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



int R2Image::
Write(const char *filename) const
{
  // Parse input filename extension
  char *input_extension;
  if (!(input_extension = (char*)strrchr(filename, '.'))) {
    fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
    return 0;
  }
  
  // Write file of appropriate type
  if (!strncmp(input_extension, ".bmp", 4)) return WriteBMP(filename);
  else if (!strncmp(input_extension, ".ppm", 4)) return WritePPM(filename, 1);
  else if (!strncmp(input_extension, ".jpg", 5)) return WriteJPEG(filename);
  else if (!strncmp(input_extension, ".jpeg", 5)) return WriteJPEG(filename);

  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



////////////////////////////////////////////////////////////////////////
// BMP I/O
////////////////////////////////////////////////////////////////////////

#if (RN_OS == RN_LINUX) && !WIN32

typedef struct tagBITMAPFILEHEADER {
  unsigned short int bfType;
  unsigned int bfSize;
  unsigned short int bfReserved1;
  unsigned short int bfReserved2;
  unsigned int bfOffBits;
} BITMAPFILEHEADER;

typedef struct tagBITMAPINFOHEADER {
  unsigned int biSize;
  int biWidth;
  int biHeight;
  unsigned short int biPlanes;
  unsigned short int biBitCount;
  unsigned int biCompression;
  unsigned int biSizeImage;
  int biXPelsPerMeter;
  int biYPelsPerMeter;
  unsigned int biClrUsed;
  unsigned int biClrImportant;
} BITMAPINFOHEADER;

typedef struct tagRGBTRIPLE {
  unsigned char rgbtBlue;
  unsigned char rgbtGreen;
  unsigned char rgbtRed;
} RGBTRIPLE;

typedef struct tagRGBQUAD {
  unsigned char rgbBlue;
  unsigned char rgbGreen;
  unsigned char rgbRed;
  unsigned char rgbReserved;
} RGBQUAD;

#endif

#define BI_RGB        0L
#define BI_RLE8       1L
#define BI_RLE4       2L
#define BI_BITFIELDS  3L

#define BMP_BF_TYPE 0x4D42 /* word BM */
#define BMP_BF_OFF_BITS 54 /* 14 for file header + 40 for info header (not sizeof(), but packed size) */
#define BMP_BI_SIZE 40 /* packed size of info header */


static unsigned short int WordReadLE(FILE *fp)
{
  // Read a unsigned short int from a file in little endian format 
  unsigned short int lsb, msb;
  lsb = getc(fp);
  msb = getc(fp);
  return (msb << 8) | lsb;
}



static void WordWriteLE(unsigned short int x, FILE *fp)
{
  // Write a unsigned short int to a file in little endian format
  unsigned char lsb = (unsigned char) (x & 0x00FF); putc(lsb, fp); 
  unsigned char msb = (unsigned char) (x >> 8); putc(msb, fp);
}



static unsigned int DWordReadLE(FILE *fp)
{
  // Read a unsigned int word from a file in little endian format 
  unsigned int b1 = getc(fp);
  unsigned int b2 = getc(fp);
  unsigned int b3 = getc(fp);
  unsigned int b4 = getc(fp);
  return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void DWordWriteLE(unsigned int x, FILE *fp)
{
  // Write a unsigned int to a file in little endian format 
  unsigned char b1 = (x & 0x000000FF); putc(b1, fp);
  unsigned char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
  unsigned char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
  unsigned char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



static int LongReadLE(FILE *fp)
{
  // Read a int word from a file in little endian format 
  int b1 = getc(fp);
  int b2 = getc(fp);
  int b3 = getc(fp);
  int b4 = getc(fp);
  return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void LongWriteLE(int x, FILE *fp)
{
  // Write a int to a file in little endian format 
  char b1 = (x & 0x000000FF); putc(b1, fp);
  char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
  char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
  char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



int R2Image::
ReadBMP(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  /* Read file header */
  BITMAPFILEHEADER bmfh;
  bmfh.bfType = WordReadLE(fp);
  bmfh.bfSize = DWordReadLE(fp);
  bmfh.bfReserved1 = WordReadLE(fp);
  bmfh.bfReserved2 = WordReadLE(fp);
  bmfh.bfOffBits = DWordReadLE(fp);
  
  /* Check file header */
  assert(bmfh.bfType == BMP_BF_TYPE);
  /* ignore bmfh.bfSize */
  /* ignore bmfh.bfReserved1 */
  /* ignore bmfh.bfReserved2 */
  assert(bmfh.bfOffBits == BMP_BF_OFF_BITS);
  
  /* Read info header */
  BITMAPINFOHEADER bmih;
  bmih.biSize = DWordReadLE(fp);
  bmih.biWidth = LongReadLE(fp);
  bmih.biHeight = LongReadLE(fp);
  bmih.biPlanes = WordReadLE(fp);
  bmih.biBitCount = WordReadLE(fp);
  bmih.biCompression = DWordReadLE(fp);
  bmih.biSizeImage = DWordReadLE(fp);
  bmih.biXPelsPerMeter = LongReadLE(fp);
  bmih.biYPelsPerMeter = LongReadLE(fp);
  bmih.biClrUsed = DWordReadLE(fp);
  bmih.biClrImportant = DWordReadLE(fp);
  
  // Check info header 
  assert(bmih.biSize == BMP_BI_SIZE);
  assert(bmih.biWidth > 0);
  assert(bmih.biHeight > 0);
  assert(bmih.biPlanes == 1);
  assert(bmih.biBitCount == 24);  /* RGB */
  assert(bmih.biCompression == BI_RGB);   /* RGB */
  int lineLength = bmih.biWidth * 3;  /* RGB */
  if ((lineLength % 4) != 0) lineLength = (lineLength / 4 + 1) * 4;
  assert(bmih.biSizeImage == (unsigned int) lineLength * (unsigned int) bmih.biHeight);

  // Assign width, height, and number of pixels
  width = bmih.biWidth;
  height = bmih.biHeight;
  npixels = width * height;

  // Allocate unsigned char buffer for reading pixels
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = bmih.biSizeImage;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Read buffer 
  fseek(fp, (long) bmfh.bfOffBits, SEEK_SET);
  if (fread(buffer, 1, bmih.biSizeImage, fp) != bmih.biSizeImage) {
    fprintf(stderr, "Error while reading BMP file %s", filename);
    return 0;
  }

  // Close file
  fclose(fp);

  // Allocate pixels for image
  pixels = new R2Pixel [ width * height ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Assign pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      double b = (double) *(p++) / 255;
      double g = (double) *(p++) / 255;
      double r = (double) *(p++) / 255;
      R2Pixel pixel(r, g, b, 1);
      SetPixel(i, j, pixel);
    }
  }

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return success
  return 1;
}



int R2Image::
WriteBMP(const char *filename) const
{
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Compute number of bytes in row
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;

  // Write file header 
  BITMAPFILEHEADER bmfh;
  bmfh.bfType = BMP_BF_TYPE;
  bmfh.bfSize = BMP_BF_OFF_BITS + rowsize * height;
  bmfh.bfReserved1 = 0;
  bmfh.bfReserved2 = 0;
  bmfh.bfOffBits = BMP_BF_OFF_BITS;
  WordWriteLE(bmfh.bfType, fp);
  DWordWriteLE(bmfh.bfSize, fp);
  WordWriteLE(bmfh.bfReserved1, fp);
  WordWriteLE(bmfh.bfReserved2, fp);
  DWordWriteLE(bmfh.bfOffBits, fp);

  // Write info header 
  BITMAPINFOHEADER bmih;
  bmih.biSize = BMP_BI_SIZE;
  bmih.biWidth = width;
  bmih.biHeight = height;
  bmih.biPlanes = 1;
  bmih.biBitCount = 24;       /* RGB */
  bmih.biCompression = BI_RGB;    /* RGB */
  bmih.biSizeImage = rowsize * (unsigned int) bmih.biHeight;  /* RGB */
  bmih.biXPelsPerMeter = 2925;
  bmih.biYPelsPerMeter = 2925;
  bmih.biClrUsed = 0;
  bmih.biClrImportant = 0;
  DWordWriteLE(bmih.biSize, fp);
  LongWriteLE(bmih.biWidth, fp);
  LongWriteLE(bmih.biHeight, fp);
  WordWriteLE(bmih.biPlanes, fp);
  WordWriteLE(bmih.biBitCount, fp);
  DWordWriteLE(bmih.biCompression, fp);
  DWordWriteLE(bmih.biSizeImage, fp);
  LongWriteLE(bmih.biXPelsPerMeter, fp);
  LongWriteLE(bmih.biYPelsPerMeter, fp);
  DWordWriteLE(bmih.biClrUsed, fp);
  DWordWriteLE(bmih.biClrImportant, fp);

  // Write image, swapping blue and red in each pixel
  int pad = rowsize - width * 3;
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      const R2Pixel& pixel = (*this)[i][j];
      double r = 255.0 * pixel.Red();
      double g = 255.0 * pixel.Green();
      double b = 255.0 * pixel.Blue();
      if (r >= 255) r = 255;
      if (g >= 255) g = 255;
      if (b >= 255) b = 255;
      fputc((unsigned char) b, fp);
      fputc((unsigned char) g, fp);
      fputc((unsigned char) r, fp);
    }

    // Pad row
    for (int i = 0; i < pad; i++) fputc(0, fp);
  }
  
  // Close file
  fclose(fp);

  // Return success
  return 1;  
}



////////////////////////////////////////////////////////////////////////
// PPM I/O
////////////////////////////////////////////////////////////////////////

int R2Image::
ReadPPM(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Read PPM file magic identifier
  char buffer[128];
  if (!fgets(buffer, 128, fp)) {
    fprintf(stderr, "Unable to read magic id in PPM file");
    fclose(fp);
    return 0;
  }

  // skip comments
  int c = getc(fp);
  while (c == '#') {
    while (c != '\n') c = getc(fp);
    c = getc(fp);
  }
  ungetc(c, fp);

  // Read width and height
  if (fscanf(fp, "%d%d", &width, &height) != 2) {
    fprintf(stderr, "Unable to read width and height in PPM file");
    fclose(fp);
    return 0;
  }
	
  // Read max value
  double max_value;
  if (fscanf(fp, "%lf", &max_value) != 1) {
    fprintf(stderr, "Unable to read max_value in PPM file");
    fclose(fp);
    return 0;
  }
	
  // Allocate image pixels
  pixels = new R2Pixel [ width * height ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for PPM file");
    fclose(fp);
    return 0;
  }

  // Check if raw or ascii file
  if (!strcmp(buffer, "P6\n")) {
    // Read up to one character of whitespace (\n) after max_value
    int c = getc(fp);
    if (!isspace(c)) putc(c, fp);

    // Read raw image data 
    // First ppm pixel is top-left, so read in opposite scan-line order
    for (int j = height-1; j >= 0; j--) {
      for (int i = 0; i < width; i++) {
        double r = (double) getc(fp) / max_value;
        double g = (double) getc(fp) / max_value;
        double b = (double) getc(fp) / max_value;
        R2Pixel pixel(r, g, b, 1);
        SetPixel(i, j, pixel);
      }
    }
  }
  else {
    // Read asci image data 
    // First ppm pixel is top-left, so read in opposite scan-line order
    for (int j = height-1; j >= 0; j--) {
      for (int i = 0; i < width; i++) {
	// Read pixel values
	int red, green, blue;
	if (fscanf(fp, "%d%d%d", &red, &green, &blue) != 3) {
	  fprintf(stderr, "Unable to read data at (%d,%d) in PPM file", i, j);
	  fclose(fp);
	  return 0;
	}

	// Assign pixel values
	double r = (double) red / max_value;
	double g = (double) green / max_value;
	double b = (double) blue / max_value;
        R2Pixel pixel(r, g, b, 1);
        SetPixel(i, j, pixel);
      }
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R2Image::
WritePPM(const char *filename, int ascii) const
{
  // Check type
  if (ascii) {
    // Open file
    FILE *fp = fopen(filename, "w");
    if (!fp) {
      fprintf(stderr, "Unable to open image file: %s", filename);
      return 0;
    }

    // Print PPM image file 
    // First ppm pixel is top-left, so write in opposite scan-line order
    fprintf(fp, "P3\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255\n");
    for (int j = height-1; j >= 0 ; j--) {
      for (int i = 0; i < width; i++) {
        const R2Pixel& p = (*this)[i][j];
        int r = (int) (255 * p.Red());
        int g = (int) (255 * p.Green());
        int b = (int) (255 * p.Blue());
        fprintf(fp, "%-3d %-3d %-3d  ", r, g, b);
        if (((i+1) % 4) == 0) fprintf(fp, "\n");
      }
      if ((width % 4) != 0) fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    // Close file
    fclose(fp);
  }
  else {
    // Open file
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
      fprintf(stderr, "Unable to open image file: %s", filename);
      return 0;
    }
    
    // Print PPM image file 
    // First ppm pixel is top-left, so write in opposite scan-line order
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255\n");
    for (int j = height-1; j >= 0 ; j--) {
      for (int i = 0; i < width; i++) {
        const R2Pixel& p = (*this)[i][j];
        int r = (int) (255 * p.Red());
        int g = (int) (255 * p.Green());
        int b = (int) (255 * p.Blue());
        fprintf(fp, "%c%c%c", r, g, b);
      }
    }
    
    // Close file
    fclose(fp);
  }

  // Return success
  return 1;  
}



////////////////////////////////////////////////////////////////////////
// JPEG I/O
////////////////////////////////////////////////////////////////////////


// #define USE_JPEG
#ifdef USE_JPEG
  extern "C" { 
#   define XMD_H // Otherwise, a conflict with INT32
#   undef FAR // Otherwise, a conflict with windows.h
#   include "jpeg/jpeglib.h"
  };
#endif



int R2Image::
ReadJPEG(const char *filename)
{
#ifdef USE_JPEG
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Initialize decompression info
  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_decompress(&cinfo);
  jpeg_stdio_src(&cinfo, fp);
  jpeg_read_header(&cinfo, TRUE);
  jpeg_start_decompress(&cinfo);

  // Remember image attributes
  width = cinfo.output_width;
  height = cinfo.output_height;
  npixels = width * height;
  int ncomponents = cinfo.output_components;

  // Allocate pixels for image
  pixels = new R2Pixel [ npixels ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Allocate unsigned char buffer for reading image
  int rowsize = ncomponents * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
    fclose(fp);
    return 0;
  }

  // Read scan lines 
  // First jpeg pixel is top-left, so read pixels in opposite scan-line order
  while (cinfo.output_scanline < cinfo.output_height) {
    int scanline = cinfo.output_height - cinfo.output_scanline - 1;
    unsigned char *row_pointer = &buffer[scanline * rowsize];
    jpeg_read_scanlines(&cinfo, &row_pointer, 1);
  }

  // Free everything
  jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);

  // Close file
  fclose(fp);

  // Assign pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      double r, g, b, a;
      if (ncomponents == 1) {
        r = g = b = (double) *(p++) / 255;
        a = 1;
      }
      else if (ncomponents == 1) {
        r = g = b = (double) *(p++) / 255;
        a = 1;
        p++;
      }
      else if (ncomponents == 3) {
        r = (double) *(p++) / 255;
        g = (double) *(p++) / 255;
        b = (double) *(p++) / 255;
        a = 1;
      }
      else if (ncomponents == 4) {
        r = (double) *(p++) / 255;
        g = (double) *(p++) / 255;
        b = (double) *(p++) / 255;
        a = (double) *(p++) / 255;
      }
      else {
        fprintf(stderr, "Unrecognized number of components in jpeg image: %d\n", ncomponents);
        return 0;
      }
      R2Pixel pixel(r, g, b, a);
      SetPixel(i, j, pixel);
    }
  }

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return success
  return 1;
#else
  fprintf(stderr, "JPEG not supported");
  return 0;
#endif
}


	

int R2Image::
WriteJPEG(const char *filename) const
{
#ifdef USE_JPEG
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Initialize compression info
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);
  jpeg_stdio_dest(&cinfo, fp);
  cinfo.image_width = width; 	/* image width and height, in pixels */
  cinfo.image_height = height;
  cinfo.input_components = 3;		/* # of color components per pixel */
  cinfo.in_color_space = JCS_RGB; 	/* colorspace of input image */
  cinfo.dct_method = JDCT_ISLOW;
  jpeg_set_defaults(&cinfo);
  cinfo.optimize_coding = TRUE;
  jpeg_set_quality(&cinfo, 95, TRUE);
  jpeg_start_compress(&cinfo, TRUE);
	
  // Allocate unsigned char buffer for reading image
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
    fclose(fp);
    return 0;
  }

  // Fill buffer with pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      const R2Pixel& pixel = (*this)[i][j];
      int r = (int) (255 * pixel.Red());
      int g = (int) (255 * pixel.Green());
      int b = (int) (255 * pixel.Blue());
      if (r > 255) r = 255;
      if (g > 255) g = 255;
      if (b > 255) b = 255;
      *(p++) = r;
      *(p++) = g;
      *(p++) = b;
    }
  }



  // Output scan lines
  // First jpeg pixel is top-left, so write in opposite scan-line order
  while (cinfo.next_scanline < cinfo.image_height) {
    int scanline = cinfo.image_height - cinfo.next_scanline - 1;
    unsigned char *row_pointer = &buffer[scanline * rowsize];
    jpeg_write_scanlines(&cinfo, &row_pointer, 1);
  }

  // Free everything
  jpeg_finish_compress(&cinfo);
  jpeg_destroy_compress(&cinfo);

  // Close file
  fclose(fp);

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return number of bytes written
  return 1;
#else
  fprintf(stderr, "JPEG not supported");
  return 0;
#endif
}






