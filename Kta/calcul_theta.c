/*--------------------------------------------------------------------------*/
/*  Name           : calcul_theta.c                                         */
/*  Version        : 1.0                                                    */
/*  Creation       : 10/07/2024                                             */
/*  Last update    : 26/07/2024                                             */
/*  Subject        : Gaussian kernel optimization                           */
/*  Algo.          : Kernel target alignment + decomposition method         */
/*  Author         : Guillaume Pradel                                       */
/*--------------------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h> 
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>

#define minimum(a,b) ((a)<=(b)?(a):(b))
#define maximum(a,b) ((a)>=(b)?(a):(b))

#define true 1
#define false 0
#define taille 81
#define very_small 1e-3
#define very_large 1e3
#define step 100
#define len_theta 7

FILE *fs, *fc, *f_theta;

int status, r=0;

long jump=false,
i, j, k, l, m, dim_input, iter, nb_data, *y, row,
y_i, y_j, Q=0, chunk_size=0, rec=100000,random_num,*table_chunk, *out_of_chunk, *in_chunk;

char fichier_data[taille], fichier_fichcom[taille],
choice='y',fichier_theta[taille];

double **X,*mean, *variance, *st_dev,lr;

double gradient_theta[len_theta];

double blosum62[25][25] = {
      { 4.91449e+00, 7.97143e-02, -1.07564e+00, -9.12302e-01, 1.31479e+00, -2.34335e-01, -1.04862e-01, 1.42341e+00, -1.03049e+00, 1.55751e-01, 1.24310e-01, -3.06781e-02, -2.15076e-01, -8.38272e-01, 3.69662e-01, 1.76739e+00, 8.59677e-01, -1.65617e+00, -1.06370e+00, 1.03310e+00, -9.82115e-01, 9.01294e-02, -1.89558e-01, -1.78863e-02, -1.36620e+00},
      { 7.97143e-02, 6.31194e+00, 1.04837e+00, -7.40637e-01, -1.42817e+00, 1.88210e+00, 1.04419e+00, -2.85762e-01, 1.15459e+00, -1.61128e+00, -6.13355e-01, 3.13549e+00, -5.61208e-02, -1.61745e+00, -3.53466e-01, -7.08956e-02, 3.00472e-02, -1.39516e+00, -8.65269e-01, -1.75771e+00, 3.16038e-01, -7.40603e-01, 1.02556e+00, 1.11650e-01, -8.61128e-01},
      { -1.07564e+00, 1.04837e+00, 7.01813e+00, 2.16109e+00, -1.68749e+00, 8.21328e-01, 9.33958e-01, 1.40044e+00, 1.97253e+00, -1.84942e+00, -1.93462e+00, 9.98793e-01, -1.24001e+00, -1.82827e+00, -6.46539e-01, 1.73846e+00, 8.43546e-01, -2.65430e+00, -1.10467e+00, -1.97592e+00, 4.83583e+00, -1.87087e+00, 7.00835e-01, -3.00013e-02, -1.38545e+00},
      { -9.12302e-01, -7.40637e-01, 2.16109e+00, 7.34256e+00, -1.44434e+00, 9.47712e-01, 3.08713e+00, 6.69362e-01, 1.49825e-01, -1.63098e+00, -2.70355e+00, 1.65769e-01, -2.09032e+00, -1.61645e+00, 6.11135e-01, 8.87485e-01, 4.01903e-03, -2.40617e+00, -1.92311e+00, -1.77973e+00, 5.07427e+00, -1.69538e+00, 1.87803e+00, 1.30667e-01, -9.00476e-01},
      { 1.31479e+00, -1.42817e+00, -1.68749e+00, -1.44434e+00, 1.09050e+01, -1.90049e+00, -2.70939e+00, -9.31509e-01, -1.59828e+00, 6.80383e-01, 6.46027e-01, -1.60678e+00, 1.34753e-01, -3.19582e-01, -1.01160e+00, 1.12548e-01, 2.44234e-01, -5.43352e-02, -6.42077e-01, 4.96560e-01, -1.48461e+00, 5.56896e-01, -1.81424e+00, 3.56831e-01, -1.99982e-01},
      { -2.34335e-01, 1.88210e+00, 8.21328e-01, 9.47712e-01, -1.90049e+00, 5.74239e+00, 2.82889e+00, -8.30210e-01, 8.11660e-01, -1.99348e+00, -1.07661e+00, 1.81283e+00, 6.50220e-01, -2.02024e+00, 1.31373e-01, 6.21433e-01, -2.94568e-01, -8.70476e-01, -2.60246e-01, -1.13945e+00, 7.43412e-01, -1.15963e+00, 4.49087e+00, -2.01755e-01, -1.82112e+00},
      { -1.04862e-01, 1.04419e+00, 9.33958e-01, 3.08713e+00, -2.70939e+00, 2.82889e+00, 5.94158e+00, -6.16253e-01, 9.51395e-01, -1.83312e+00, -1.90730e+00, 1.95047e+00, -1.23438e+00, -1.85464e+00, 3.34674e-01, 7.37458e-01, -1.65750e-01, -1.67761e+00, -1.11343e+00, -9.91452e-01, 1.93234e+00, -1.98732e+00, 4.65182e+00, -7.06421e-02, -1.43747e+00},
      { 1.42341e+00, -2.85762e-01, 1.40044e+00, 6.69362e-01, -9.31509e-01, -8.30210e-01, -6.16253e-01, 8.25310e+00, -4.79449e-01, -2.18011e+00, -2.20113e+00, -4.93309e-01, -1.76546e+00, -1.17873e+00, 1.63615e-01, 1.21374e+00, -6.45582e-01, 1.10691e-01, -1.51462e+00, -1.37485e+00, 6.91574e-01, -2.30642e+00, -6.70150e-01, 4.59015e-01, 1.26500e-01},
      { -1.03049e+00, 1.15459e+00, 1.97253e+00, 1.49825e-01, -1.59828e+00, 8.11660e-01, 9.51395e-01, -4.79449e-01, 9.03204e+00, -1.76489e+00, -1.79325e+00, 2.79031e-02, -1.16513e+00, 2.37203e-01, -5.37879e-01, -1.82080e-01, -1.08450e+00, -5.68059e-01, 2.99821e+00, -1.89900e+00, 1.10408e+00, -1.84891e+00, 8.68723e-01, 9.99438e-03, -1.20130e+00},
      { 1.55751e-01, -1.61128e+00, -1.84942e+00, -1.63098e+00, 6.80383e-01, -1.99348e+00, -1.83312e+00, -2.18011e+00, -1.76489e+00, 5.50860e+00, 3.47003e+00, -1.78382e+00, 2.00637e+00, 1.48297e+00, -1.24862e+00, -1.01885e+00, 9.38481e-02, -1.28100e+00, 1.87342e-01, 4.32504e+00, -1.65868e+00, 4.30595e+00, -2.01908e+00, 1.83678e-01, -6.56562e-01},
      { 1.24310e-01, -6.13355e-01, -1.93462e+00, -2.70355e+00, 6.46027e-01, -1.07661e+00, -1.90730e+00, -2.20113e+00, -1.79325e+00, 3.47003e+00, 5.48386e+00, -8.28872e-01, 2.99756e+00, 1.44482e+00, -1.27239e+00, -1.01832e+00, 7.97753e-02, -3.18583e-01, 1.94341e-01, 2.30877e+00, -2.56526e+00, 4.26386e+00, -1.92344e+00, 1.51121e-01, -7.14770e-01},
      { -3.06781e-02, 3.13549e+00, 9.98793e-01, 1.65769e-01, -1.60678e+00, 1.81283e+00, 1.95047e+00, -4.93309e-01, 2.79031e-02, -1.78382e+00, -8.28872e-01, 6.03749e+00, -1.78803e-01, -1.76705e+00, 4.49489e-01, 8.04103e-01, -9.21421e-02, -1.57696e+00, -1.01591e+00, -9.12131e-01, 1.03999e+00, -1.81273e+00, 1.85061e+00, 1.30049e-02, -1.21824e+00},
      { -2.15076e-01, -5.61208e-02, -1.24001e+00, -2.09032e+00, 1.34753e-01, 6.50220e-01, -1.23438e+00, -1.76546e+00, -1.16513e+00, 2.00637e+00, 2.99756e+00, -1.78803e-01, 5.68958e+00, 9.98317e-01, -8.13547e-01, -3.26174e-01, -2.52882e-01, 1.59223e-01, -1.78909e-01, 1.89655e+00, -2.05112e+00, 2.90536e+00, -2.77644e-01, -1.40793e-01, -1.72058e+00},
      { -8.38272e-01, -1.61745e+00, -1.82827e+00, -1.61645e+00, -3.19582e-01, -2.02024e+00, -1.85464e+00, -1.17873e+00, 2.37203e-01, 1.48297e+00, 1.44482e+00, -1.76705e+00, 9.98317e-01, 7.48414e+00, -2.24819e+00, -1.02220e+00, -9.04588e-01, 2.71715e+00, 4.19203e+00, 3.19784e-01, -1.68941e+00, 1.37592e+00, -1.97206e+00, 1.98772e-01, -6.49329e-01},
      { 3.69662e-01, -3.53466e-01, -6.46539e-01, 6.11135e-01, -1.01160e+00, 1.31373e-01, 3.34674e-01, 1.63615e-01, -5.37879e-01, -1.24862e+00, -1.27239e+00, 4.49489e-01, -8.13547e-01, -2.24819e+00, 9.07857e+00, 1.65622e-01, 3.00802e-01, -1.97019e+00, -1.57512e+00, -4.36657e-01, -3.87055e-01, -1.37586e+00, 2.66711e-01, 4.08010e-01, -3.28494e-02},
      { 1.76739e+00, -7.08956e-02, 1.73846e+00, 8.87485e-01, 1.12548e-01, 6.21433e-01, 7.37458e-01, 1.21374e+00, -1.82080e-01, -1.01885e+00, -1.01832e+00, 8.04103e-01, -3.26174e-01, -1.02220e+00, 1.65622e-01, 4.66183e+00, 1.73263e+00, -1.86427e+00, -1.19094e+00, -1.12029e+00, 9.42915e-01, -1.10411e+00, 7.35260e-01, -1.76282e-01, -1.76758e+00},
      { 8.59677e-01, 3.00472e-02, 8.43546e-01, 4.01903e-03, 2.44234e-01, -2.94568e-01, -1.65750e-01, -6.45582e-01, -1.08450e+00, 9.38481e-02, 7.97753e-02, -9.21421e-02, -2.52882e-01, -9.04588e-01, 3.00802e-01, 1.73263e+00, 5.81727e+00, -7.30273e-01, -1.10334e+00, 9.77587e-01, 1.93788e-02, 2.27560e-02, -1.99612e-01, -8.79052e-02, -1.50867e+00},
      { -1.65617e+00, -1.39516e+00, -2.65430e+00, -2.40617e+00, -5.43352e-02, -8.70476e-01, -1.67761e+00, 1.10691e-01, -5.68059e-01, -1.28100e+00, -3.18583e-01, -1.57696e+00, 1.59223e-01, 2.71715e+00, -1.97019e+00, -1.86427e+00, -7.30273e-01, 1.29881e+01, 3.38440e+00, -1.46994e+00, -2.46175e+00, -4.16666e-01, -7.99884e-01, 3.90477e-01, -1.18510e-01},
      { -1.06370e+00, -8.65269e-01, -1.10467e+00, -1.92311e+00, -6.42077e-01, -2.60246e-01, -1.11343e+00, -1.51462e+00, 2.99821e+00, 1.87342e-01, 1.94341e-01, -1.01591e+00, -1.78909e-01, 4.19203e+00, -1.57512e+00, -1.19094e+00, -1.10334e+00, 3.38440e+00, 7.99398e+00, 7.08361e-02, -1.83576e+00, 1.18456e-01, -1.07002e+00, 3.09899e-03, -1.27319e+00},
      { 1.03310e+00, -1.75771e+00, -1.97592e+00, -1.77973e+00, 4.96560e-01, -1.13945e+00, -9.91452e-01, -1.37485e+00, -1.89900e+00, 4.32504e+00, 2.30877e+00, -9.12131e-01, 1.89655e+00, 3.19784e-01, -4.36657e-01, -1.12029e+00, 9.77587e-01, -1.46994e+00, 7.08361e-02, 5.18364e+00, -1.79370e+00, 3.19922e+00, -1.05907e+00, 7.83166e-02, -1.00963e+00},
      { -9.82115e-01, 3.16038e-01, 4.83583e+00, 5.07427e+00, -1.48461e+00, 7.43412e-01, 1.93234e+00, 6.91574e-01, 1.10408e+00, -1.65868e+00, -2.56526e+00, 1.03999e+00, -2.05112e+00, -1.68941e+00, -3.87055e-01, 9.42915e-01, 1.93788e-02, -2.46175e+00, -1.83576e+00, -1.79370e+00, 5.59633e+00, -1.84780e+00, 1.18799e+00, 6.85294e-02, -9.53277e-01},
      { 9.01294e-02, -7.40603e-01, -1.87087e+00, -1.69538e+00, 5.56896e-01, -1.15963e+00, -1.98732e+00, -2.30642e+00, -1.84891e+00, 4.30595e+00, 4.26386e+00, -1.81273e+00, 2.90536e+00, 1.37592e+00, -1.37586e+00, -1.10411e+00, 2.27560e-02, -4.16666e-01, 1.18456e-01, 3.19922e+00, -1.84780e+00, 4.47162e+00, -1.93141e+00, 1.51254e-01, -8.78017e-01},
      { -1.89558e-01, 1.02556e+00, 7.00835e-01, 1.87803e+00, -1.81424e+00, 4.49087e+00, 4.65182e+00, -6.70150e-01, 8.68723e-01, -2.01908e+00, -1.92344e+00, 1.85061e+00, -2.77644e-01, -1.97206e+00, 2.66711e-01, 7.35260e-01, -1.99612e-01, -7.99884e-01, -1.07002e+00, -1.05907e+00, 1.18799e+00, -1.93141e+00, 5.09145e+00, -1.35928e-01, -1.59558e+00},
      { -1.78863e-02, 1.11650e-01, -3.00013e-02, 1.30667e-01, 3.56831e-01, -2.01755e-01, -7.06421e-02, 4.59015e-01, 9.99438e-03, 1.83678e-01, 1.51121e-01, 1.30049e-02, -1.40793e-01, 1.98772e-01, 4.08010e-01, -1.76282e-01, -8.79052e-02, 3.90477e-01, 3.09899e-03, 7.83166e-02, 6.85294e-02, 1.51254e-01, -1.35928e-01, 3.85228e-01, -1.17991e+00},
      { -1.36620e+00, -8.61128e-01, -1.38545e+00, -9.00476e-01, -1.99982e-01, -1.82112e+00, -1.43747e+00, 1.26500e-01, -1.20130e+00, -6.56562e-01, -7.14770e-01, -1.21824e+00, -1.72058e+00, -6.49329e-01, -3.28494e-02, -1.76758e+00, -1.50867e+00, -1.18510e-01, -1.27319e+00, -1.00963e+00, -9.53277e-01, -8.78017e-01, -1.59558e+00, -1.17991e+00, 8.61426e+00},
};

double deux_sigma_carre = 8*7;
double theta[] = {0.2, 0.2, 0.2, 1 , 0.2, 0.2, 0.2};

/* Functions included in this program */

void caract_db();
void alloc_memory();
void read_data();
void standardize_data();
void gradient_norme_ktheta_ktarget();
void gradient_moins_cosinus();
void pas_de_gradient();
void descente_en_gradient();
void display_theta();
void write_theta();double **matrix(int nrow, int ncol);
double **matrix(int nrow, int ncol);
double moins_cosinus();


int main(int argc, char *argv[])

{

strcpy(fichier_fichcom, argv[1]);
srand(time(NULL));

caract_db();
read_data();

printf("\nDo you want data to be standardized [y/n]? ");
choice = getc(stdin);
if((choice == 'y') || (choice == 'Y'))
  standardize_data();


printf("\n\n*** Calcul de theta: \n");

descente_en_gradient();
display_theta();
write_theta();


}

double **matrix(int nrow, int ncol)

/* Allocation d'une matrice Cf. Numerical Recipes in C */

{

int ind1, ind2;

double **m;

m = (double **) malloc((size_t)((nrow+1)*sizeof(double*)));
if(!m)
  {
  printf("\nallocation failure 1 in matrix()");
  exit(0);
  }

m[1] = (double *) malloc((size_t)((nrow*(ncol+1))*sizeof(double)));
if(!m[1])
  {
  printf("\nallocation failure 2 in matrix()");
  exit(0);
  }

for(ind1=2;ind1<=nrow;ind1++)
  m[ind1]=m[ind1-1]+ncol;

return m;

}


void caract_db()

/* Lecture des variables dans Fichcom */

{

if((fs=fopen(fichier_fichcom, "r"))==NULL)
  {
  printf("\nFile of parameters %s: cannot be open...\n", fichier_fichcom);
  exit(0);
  }

status = fscanf(fs, "%lf", &lr);
status = fscanf(fs, "%ld", &rec);
status = fscanf(fs, "%ld", &chunk_size);
status = fscanf(fs, "%s", fichier_data);
status = fscanf(fs, "%s", fichier_theta);
printf("\nThe data file is: %s\n", fichier_data);

fclose(fs);



}

void alloc_memory()

{

out_of_chunk = (long *) calloc(nb_data+1, sizeof(long));
in_chunk = (long *) calloc(nb_data+1, sizeof(long));
table_chunk = (long *) calloc(chunk_size+1, sizeof(long));
X = matrix(nb_data, dim_input);
mean = (double *) calloc(dim_input+1, sizeof(double));
variance = (double *) calloc(dim_input+1, sizeof(double));
st_dev = (double *) calloc(dim_input+1, sizeof(double));
y = (long *) calloc(nb_data+1, sizeof(long));

}

void read_data()

/* Lecture des données pour l'apprentissage de theta */

{

long min_y, max_y;
double value_y;

if((fs=fopen(fichier_data, "r"))==NULL)
  {
  printf("\nData file %s: cannot be open...\n", fichier_data);
  exit(0);
  }

status = fscanf(fs, "%ld", &nb_data);
status = fscanf(fs, "%ld", &dim_input);
status = fscanf(fs, "%ld", &Q);

min_y = nb_data;
max_y = 0;

alloc_memory();

for(i=1; i<=nb_data; i++)
  {
  for(j=1; j<=dim_input; j++)
    status = fscanf(fs, "%lf", &X[i][j]);
  status = fscanf(fs, "%lf", &value_y);
  y[i] = (long) value_y;
  if(y[i] < min_y)
     min_y = y[i];
  if(y[i] > max_y)
    max_y = y[i];
  }

fclose(fs);

if((min_y != 1) || (max_y!= Q))
  {
  printf("\nWrong numbering of the categories\n");
  exit(0);
  }

}

void standardize_data()

/* Standardize the data of the training set per predictor */

{

for(j=1; j<=dim_input; j++)
  {
  mean[j] = 0.0;
  variance[j] = 0.0;
  st_dev[j] = 0.0;
  }

for(i=1; i<=nb_data; i++)
  for(j=1; j<=dim_input; j++)
    mean[j] += X[i][j];

for(j=1; j<=dim_input; j++)
  {
  mean[j] /= (double) nb_data;
  }

for(i=1; i<=nb_data; i++)
  for(j=1; j<=dim_input; j++)
    variance[j] += (X[i][j] - mean[j]) * (X[i][j] - mean[j]);

for(j=1; j<=dim_input; j++)
  {
  variance[j] /= (double) nb_data;
  st_dev[j] = sqrt(variance[j]);
  }

for(i=1; i<=nb_data; i++)
  for(j=1; j<=dim_input; j++)
    {
    X[i][j] -= mean[j];
    if(st_dev[j] > very_small)
      X[i][j] /= st_dev[j];
    }

}


void gradient_norme_ktheta_ktarget() 

/* Calcul du gradient de la fonction : || Ko - Kt || ** 2   */

{
  // Initialisation du gradient à zéro
  for (int k = 0; k < len_theta; k++) {
      gradient_theta[k] = 0;
  }

  // Boucle principale pour chaque dérivée partielle
  for (int k = 0; k < len_theta; k++) {
    if (k != len_theta / 2) {

      double somme = 0;

      // Boucles pour parcourir les données dans le chunk
      for (int i = 1; i <= chunk_size; i++) {
        for (int j = 1; j <= chunk_size; j++)  {

          double somme_exp = 0;
          long ij, ii, jj;

          // Calcul de la somme dans l'exponentielle
          for (int l = 0; l < len_theta; l++) {

          ii = X[table_chunk[i]][l];
          jj = X[table_chunk[j]][l];

            somme_exp -= theta[l] * theta[l] * (blosum62[ii][ii]  + blosum62[jj][jj]  - 2 * blosum62[ii][jj]);
          }                 

          double expo = exp(somme_exp / (deux_sigma_carre)); // Vaut ko(x^i, x^j)


          ii = X[table_chunk[i]][k];
          jj = X[table_chunk[j]][k];

          double derive_expo = ((-2) * theta[k] / deux_sigma_carre ) * (blosum62[ii][ii] + blosum62[jj][jj] - 2 * blosum62[ii][jj]) * expo;

          if (y[table_chunk[i]] == y[table_chunk[j]]) {
            somme += (expo - 1) * derive_expo;
          } else {
            somme += (expo) * derive_expo;
          }
        }
      }
      gradient_theta[k] = 2*somme;
    }
  }
}

void gradient_moins_cosinus() 

/* Calcul du gradient de la fonction objectif 'cosinus'

On Kernel-Target Alignment
N. Cristianini, J. Shawe-Taylor, A. Elisseeff, and J. Kandola.  */

{
  // Initialisation du gradient à zéro
  for (int k = 0; k < len_theta; k++) {
      gradient_theta[k] = 0;
  }

  // Boucle principale pour le calcul de chaque dérivée partielle
  for (int k = 0; k < len_theta; k++) {
    if (k != len_theta / 2) {

      double ps = 0;
      double norme_kt = 0;
      double norme_ko = 0;
      double derive_norme_ko = 0;
      double derive_ps = 0;
      double tmp = 0;

      // Boucles pour parcourir les données dans le chunk
      for (int i = 1; i <= chunk_size; i++) {
        for (int j = 1; j <= chunk_size; j++) {

          double somme_exp = 0;
          long ii, jj;


          // Boucle pour parcourir les composantes des données
          for (int l = 0; l < len_theta; l++) {

            ii = X[table_chunk[i]][l];
            jj = X[table_chunk[j]][l];


            somme_exp -= theta[l] * theta[l] * (blosum62[ii][ii]  + blosum62[jj][jj]  - 2.0 * blosum62[ii][jj]);
          }

          double expo = exp(somme_exp / deux_sigma_carre); // Vaut ko(x^i, x^j)

          norme_ko += (expo * expo);

          ii = X[table_chunk[i]][k];
          jj = X[table_chunk[j]][k];

          tmp += ((-4.0) * theta[k] / deux_sigma_carre ) * (blosum62[ii][ii] + blosum62[jj][jj] - 2.0 * blosum62[ii][jj]) * expo * expo;

          if (y[table_chunk[i]] == y[table_chunk[j]]) {
            ps += expo;
            norme_kt += 1;
            derive_ps += ((-2.0) * theta[k] / deux_sigma_carre ) * (blosum62[ii][ii] + blosum62[jj][jj] - 2.0 * blosum62[ii][jj]) * expo;
          }
        }
      }

      norme_kt = sqrt(norme_kt);
      norme_ko = sqrt(norme_ko);
      derive_norme_ko =  (1.0 / 2.0) * (1.0 / norme_ko) * tmp;
      gradient_theta[k] = -( (derive_ps * norme_ko) - (ps * derive_norme_ko) ) / ((norme_ko * norme_ko) * (norme_kt));
    }
  }
}



void pas_de_gradient() 

/* Itération d'un pas de gradient sur un chunk */

{
  for (int j = 0; j < len_theta; j++) {
      gradient_theta[j] = 1.0;
  }
  gradient_theta[len_theta / 2] = 0.0;

  gradient_moins_cosinus();
  for (int j = 0; j < len_theta; j++) {
      theta[j] = theta[j] - lr * gradient_theta[j];
  }
}



double moins_cosinus() 

/* Calcul de la fonction objectif de l'inverse du cosinus */

{
  double ps = 0;
  double norme_kt = 0;
  double norme_ko = 0;

  for (int i = 1; i <= nb_data; i++) {
    for (int j = 1; j <= nb_data; j++) {

          double somme_exp = 0;
          long ii, jj;


          // Calcul de somme_exp
          for (int l = 0; l < len_theta; l++) {

              ii = X[i][l];
              jj = X[j][l];

              somme_exp -= theta[l] * theta[l] * (blosum62[ii][ii]  + blosum62[jj][jj]  - 2 * blosum62[ii][jj]);

          }

          double expo = exp(somme_exp / deux_sigma_carre);

          norme_ko += expo * expo;

          if (y[i] == y[j]) {
            ps += expo;
            norme_kt += 1;
          }
    }
  }
  norme_kt = sqrt(norme_kt);
  norme_ko = sqrt(norme_ko);
  return -( ps ) / (norme_ko * norme_kt);
}

double norme_ktheta_ko() 

/* Calcul de la fonction objectif : || Ko - Kt || ** 2 */

{
  double ps = 0;
  double norme_kt = 0;
  double norme_ko = 0;

  for (int i = 1; i <= nb_data; i++) {
    for (int j = 1; j <= nb_data; j++) {

          double somme_exp = 0;
          long ii, jj;
          // Calcul de somme_exp
          for (int l = 0; l < len_theta; l++) {

              ii = X[i][l];
              jj = X[j][l];

              somme_exp -= theta[l] * theta[l] * (blosum62[ii][ii]  + blosum62[jj][jj]  - 2 * blosum62[ii][jj]);

          }

          double expo = exp(somme_exp / deux_sigma_carre);

          norme_ko += expo * expo;

          if (y[i] == y[j]) {
            ps += expo;
            norme_kt += 1;
          }
    }
  }

  return norme_ko + norme_kt - 2 * ps;


}


void compute_table_chunk()

{
  for(i=1; i<=nb_data; i++)
    {
    out_of_chunk[i] = i;
    in_chunk[i] = 0;
    }

  for(i=1; i<=chunk_size; i++)
    {
    if(i == 1)
      {
      row = iter % nb_data;
      if(row == 0)
        row = nb_data;
      }
    else
      {
      random_num = rand();
      row = (random_num % (nb_data-i+1))+1;
      }
    table_chunk[i] = out_of_chunk[row];
    for(j=row; j<=nb_data-i; j++)
      out_of_chunk[j] = out_of_chunk[j+1];
    in_chunk[table_chunk[i]] = 1;
    }
}



void descente_en_gradient() 

/* Calcul de la descente en gradient en effectuant 'rec' pas de gradient */

{
  for (iter = 1; iter <= rec; iter++) {

    // if( iter % 1000 == 0) {
    //   lr = 0.995*lr;
    // }
    

    if ( iter == rec ) {
      printf("La valeur de la fonction objectif est : %lf \n", moins_cosinus());
      display_theta();
    }

    compute_table_chunk();
    pas_de_gradient();
  }
}



void display_theta() 

/* Fonction permettant d'afficher la valeur de theta */

{
  int i = 0;
  printf(" Affichage de theta : ");
  for(i = 0; i < len_theta; i++) {
    printf("%lf  ", theta[i]);
  }
  printf("\n\n");
}



void write_theta()

/* Ecriture de theta dans le fichier tableau_theta.txt */

{
  f_theta = fopen(fichier_theta,"w+");

  if(f_theta == NULL)
    {
    printf("\nFile of parameters fichier_theta : cannot be open...\n");
    exit(0);
    }

  for (int i = 0; i < len_theta; i++) {
      fprintf(f_theta, "%lf", theta[i]);
      if (i < len_theta - 1) {
          fprintf(f_theta, "\n");
      }
  }

  fclose(f_theta);
}


