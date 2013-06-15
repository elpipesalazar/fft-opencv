#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <valarray>
#include <complex>
#include "highgui.h"
#include "cv.h"

const double PI = acos(-1);
using namespace std;
 
typedef complex<double> Complex;
typedef valarray<Complex> CArray;
typedef float pix;

// FFT Cooley–Tukey 
void fft(CArray& x)
{
    const size_t N = x.size();
    if (N <= 1) return;
 
    // divide
    CArray par = x[slice(0, N/2, 2)];
    CArray  impar = x[slice(1, N/2, 2)];
 
    // conquer
    fft(par);
    fft(impar);
 
    // combine
    for (size_t k = 0; k < N/2; ++k)
    {
        Complex t = polar(1.0, -2 * PI * k / N) * impar[k];
        x[k    ] = par[k] + t;
        x[k+N/2] = par[k] - t;
    }
}
 
// FFT Inversa
void ifft(CArray& x)
{
    // conjugate the complex numbers
    transform(&x[0], &x[x.size()], &x[0], conj<double>);
 
    // forward fft
    fft( x );
 
    // conjugate the complex numbers again
    transform(&x[0], &x[x.size()], &x[0], conj<double>);
 
    // scale the numbers
    x /= x.size();
}

int main(int argc, char **argv){
    if(argc < 2){
        cout<<endl<<"Usage: "<<argv[0]<<" <filename.jpg>"<<endl;
        return 1;
    }
    
    cv::Mat filter = cv::getGaussianKernel(9,0.5,CV_32F);
    cv::Mat image;
    char *imageName = argv[1];
    image = cv::imread( imageName, 1 );

    cv::Mat gray_image;
    cv::cvtColor( image, gray_image, CV_BGR2GRAY );

    cv::imwrite( "./Gray_Image.jpg", gray_image);

    /*
    * This section shows the imgage 
    cv::namedWindow( imageName, CV_WINDOW_AUTOSIZE );
    cv::namedWindow( "Gray image", CV_WINDOW_AUTOSIZE );

    cv::imshow( imageName, image );
    cv::imshow( "Gray image", gray_image );

    cv::waitKey(0);*/
    return 0;
}
