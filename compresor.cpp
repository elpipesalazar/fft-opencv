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


void img2mat(CvMat *mat, IplImage *img, int ch) {
    int n = img->height, m = img->width, nch = img->nChannels;
    for(int i = 0; i < n; ++i){
        uchar *orow = (uchar*)(img->imageData + i*img->widthStep);
        pix *drow = (pix*)(mat->data.ptr + i*mat->step);
        for(int j = 0; j< m; ++j)
            drow[j] = orow[j*nch + ch];
    }
}

void mat2img(IplImage *img, CvMat* mat, int ch) {
    int n = mat->height, m = mat->width, nch = img->nChannels;
    for(int i = 0; i< n; ++i){
        uchar *drow = (uchar*)(img->imageData + i*img->widthStep);
        pix *orow = (pix*)(mat->data.ptr + i*mat->step);
        for(int j = 0; j< m; ++j){
            pix x = orow[j];
            if (x<0) x=0; if(x>255) x=255;
            drow[j*nch + ch] = x;
        }
    }
}



int main(int argc, char **argv){
    if(argc < 2){
        printf("\nUsage: %s <filename.jpg>\n",argv[0]);
        return 1;
    }
    
    cv::Mat filter = cv::getGaussianKernel(9,0.5,CV_32F);
    cv::Mat image;
    char *imageName = argv[1];
    image = cv::imread( imageName, 1 );

    cv::Mat gray_image;
    cv::cvtColor( image, gray_image, CV_BGR2GRAY );

    cv::imwrite( "./Gray_Image.jpg", gray_image);

    cv::namedWindow( imageName, CV_WINDOW_AUTOSIZE );
    cv::namedWindow( "Gray image", CV_WINDOW_AUTOSIZE );

    cv::imshow( imageName, image );
    cv::imshow( "Gray image", gray_image );

    cv::waitKey(0);

    /*
    // Display kernel values
	cv::Mat_<float>::const_iterator it= filter.begin<float>();  
	cv::Mat_<float>::const_iterator itend= filter.end<float>();  
	cout << "Valores del Kernel\n" << "[";
	for ( ; it!= itend; ++it) 
		cout << *it << " ";
	cout << "]\n" << endl;

    /*

    // // forward fft
    // fft(data);
 
    // cout << "fft" << endl;
    // for (int i = 0; i < img->width; ++i)
    // {
    //     cout << data[i] << endl;
    // }
 
    // // inverse fft
    // ifft(data);
 
    // cout << endl << "ifft" << endl;
    // for (int i = 0; i < 8; ++i)
    // {
    //     cout << data[i] << endl;
    // }
    */
    return 0;
}
