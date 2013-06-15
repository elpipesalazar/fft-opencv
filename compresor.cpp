#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <valarray>
#include <complex>
#include "highgui.h"
#include "cv.h"

const double PI = 3.141592653589793238460;
 
typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

using namespace std;


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
    int pos = 0;
    IplImage *img = cvLoadImage(argv[1]);
    CvMat *matImg = cvCreateMat(img->height,img->width,CV_32FC3 );
    cvConvert( img, matImg );

    cv::Mat filter = cv::getGaussianKernel(9,0.5,CV_32F);
    
	Complex test[1000];

    for(int i=0;i<1000;i++)
    {
        for(int j=0;j<1000;j++)
        {
            CvScalar scal = cvGet2D( matImg,j,i);
            test[pos] =  std::complex<double>(0.0) ;
            test[pos + 1] =  std::complex<double>(scal.val[1]);
            test[pos + 2] =  std::complex<double>(scal.val[2]);
            pos = pos + 3 ;
        }
    }

    //CArray data(test, img->width*3);





    // Display kernel values
	cv::Mat_<float>::const_iterator it= filter.begin<float>();  
	cv::Mat_<float>::const_iterator itend= filter.end<float>();  
	cout << "Valores del Kernel\n" << "[";
	for ( ; it!= itend; ++it) {
		cout << *it << " ";
	}
	cout << "]\n" << endl;


	for(int i=0;i<10;i++)
	{
	    for(int j=0;j<10;j++)
	    {
	        CvScalar scal = cvGet2D( matImg,j,i);
	        printf( "(%.f,%.f,%.f)  ",scal.val[0], scal.val[1], scal.val[2] );
	    }
	    printf("\n");
	}

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
    return 0;
}