import java.io.*;
import static java.lang.Math.pow;

class AGF {
    // Read base from file.
    public static double[][] base;
    // The centers are represented as a two dimensional array.
    // The number of rows represent the number of centers, while
    // the number of columns represent the number of coordinates.
    // Las letras representan las coordenadas en el espacio 12-dimensional de los centros
    public  double F01(double A1, double A2, double A3, double A4,  double A5,  double A6,
                       double A7, double A8, double A9, double A10, double A11, double A12,
                       double B1, double B2, double B3, double B4,  double B5,  double B6,
                       double B7, double B8, double B9, double B10, double B11, double B12,
                       double C1, double C2, double C3, double C4,  double C5,  double C6,
                       double C7, double C8, double C9, double C10, double C11, double C12,
                       double D1, double D2, double D3, double D4,  double D5,  double D6,
                       double D7, double D8, double D9, double D10, double D11, double D12,
                       double E1, double E2, double E3, double E4,  double E5,  double E6,
                       double E7, double E8, double E9, double E10, double E11, double E12)
    {
        // Generate the centers matrix
        double[][] centers = {
            {A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12},
            {B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12},
            {C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11, C12},
            {D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12},
            {E1, E2, E3, E4, E5, E6, E7, E8, E9, E10, E11, E12},
        };
        int        nrow          = base.length;       // Number of rows of base.
        double[]   dist_ind      = new double[nrow];  // Array of distance between individuals and centers.
        int        distCluster   = 0;                 // Keeps record of final distance to centers.
        double     min           = 1e12;              // Minimum distance found between individuals and centers.
        int[]      clust_class   = new int[nrow];     // Array of  classes of individuals.
        int        sum           = 0;                 // Keeps record of distance to centers.
        double[]   dist_mean     = new double[nrow];  // Array of mean distance to centers.
        double[]   dist_var      = new double[nrow];  // Array of variance of distances to centers.
        int        in_cluster    = 0;                 // Keeps record of individuals in clusters.
        double     total_var     = 0;
        for(int i = 0; i < (base.length - 1); i++){                   // This index runs over all the individuals.
            for(int k = 0; k < (centers.length - 1); k++){            // This index runs over all the centers.
                for(int j = 0; j < (base[0].length - 1); j++){        // This index runs over all the observations of the base.
                    sum = sum + pow((base[i][j]) - centers[k][j], 2); // Distance from observation i to center k.
                }
                // Find the minimum distance among all the centers.
                if(sum < min){
                    min            = sum;
                    distCluster    = min;
                    clust_class[i] = k;
                }
                sum = 0;
            }
            dist_ind[i] = distCluster;
        }
        // Fill in dist_mean and dist_var
        for(int k = 0; k < (centers.length - 1); k++){ // This index runs over all the centers.
            dist_mean[k] = 0;
            dist_var[k]  = 0;
        }
        // Calculate mean per cluster.
        for(int k = 0; k < (centers.length - 1); k++){   // This index runs over all the centers.
            for(int i = 0; i < (base.length - 1); i++){  // This index runs over all the individuals.
                if(clust_class[i] == k){
                    dist_mean[k] =  dist_mean[k] + dist_ind[i]; // Sums over all the members of the cluster.
                    in_cluster   = in_cluster + 1;              // Counts every individual in the cluster.
                }
            }
            in_cluster   = 0;
            dist_mean[k] = dist_mean[k] / in_cluster;
        }
        // Calculate var per cluster.
        for(int k = 0; k < (centers.length - 1); k++){   // This index runs over all the centers.
            for(int i = 0; i < (base.length - 1); i++){  // This index runs over all the individuals.
                if(clust_class[i] == k){
                    dist_var[k] =  dist_var[k] + Math.pow(dist_mean[k] - dist_ind[i], 2); // Sums the square of differences
                    in_cluster   = in_cluster + 1;                                        // Counts every individual in the cluster.
                }
            }
            in_cluster   = 0;
            dist_var[k]  = dist_var[k] / in_cluster;
            total_var    = total_var + dist_var[k];
        }
        return total_var;
    }  //endF01


    
  public  double F02(double X){
  double Y;
  /*
   * (2)
   *
   *	Resolver la siguiente ecuacion cubica:
   *
   *	Y=X^3+2.5X^2-2X+1
   *
   *	MINIMIZAR CON LAS RESTRICCIONES
   *
   *  X=-3.2180565 --> Y=0.0000000
   *  N=200
   *  E=4
   *  D=40
   *  Pc=1.000
   *  Pm=.005
   *  G=1000
   */
   Y=X*X*X+2.5*X*X-2*X+1;
   if (Y>0&&Y<=1) return Y;
   if (Y<1)     return 1000d;
   if (Y<10)    return 10000d;
   if (Y<100)   return 100000d;
   if (Y<1000)  return 1000000d;
   if (Y<10000) return 100000000d;
   return 10000000000d;
  }//endF02

  public  double F03(double X,double Y){
/*
 * (3)
 *
 *	Maximizar
 *
 *	N=200
 *	E=4
 *	D=40
 *	Pc=1
 *	Pm=.005
 *	G=5000
 *
 *	X=5
 *	Y=3
 *	F(X,Y)=0
 *
 */
 	return -(X-5)*(X-5)-(Y-3)*(Y-3);
  }//endF03

  public  double F04(double X,double Y){
/*
 * (4)
 *
 *	Minimizar
 *
 *	Z --> 0
 *	N=200
 *	E=4
 *	D=40
 *	Pc=1
 *	Pm=.005
 *	G=1000
 *
 *	(X,Y) = (1.198231,0.706163); Z=0.00000
 *
 *
 */
	double Z;
	Z=Math.pow(X,2)-2*X*Y+Math.pow(Y,3)+X+Y-2; 
	if (Z<0){
		if (Z>-100)
  			return 100000;
  		else
  			if (Z>-1000)
  				return 1000000;
  			else
  				if (Z>-10000)
  					return 100000000;
  				else
  					return 1000000000;
  				//endif
  			//endif
  		//endif
	}//endif
	return Z;
  }//endF04

  public  double F05(double X, double Y){
/*
 * (5)
 *
 */
//  2X+3Y-13=0
//   X-2Y+ 4=0
//  ----------
//  3X+ Y- 9=0
//  X=2; Y=3
/*
 *	USAR:
 *		  2 BITS ENTEROS
 *		 40 BITS DECIMALES
 *		200 INDIVIDUOS
 *		400 GENERACIONES
 *		Pc: 0.9
 *		Pm: 0.02
 */
 	double R1,R2;
	int C=0;
	R1=2*X+3*Y-13;
	R2=X-2*Y+4;
	if (R1>=0) C=C+1;
	if (R2>=0) C=C+1;
	if (C==2) return 3*X+Y-9;
	if (C==1) return 5000000d;
	return 10000000d;
  }//endF05
} //endClass
