#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

#define MAX(x,y) (x>y)?x:y

/**
 * Print a message before exit
 */
void exit_msg(const char *msg) {
	perror(msg);
	exit(EXIT_FAILURE);
}


/**
 * Simulate the reading of an image
 */
int **getImage(const int height, const int width) {
	int **matrix = malloc(sizeof(int *)*height);
	
	if (matrix == NULL)
		exit_msg("getImage: cannot allocate memory (1)");
	
	int i,j;
	for (i = 0; i < height; i++) {
		if ((matrix[i] = malloc(sizeof(int) * width)) == NULL)
			exit_msg("getImage: cannot allocate memory (2)");
		for (j = 0; j < width; j++)
			matrix[i][j] = rand() % 256; //grey level
	}
	
	return matrix;
}

/**
 * only print
 */
void print_image(int **matrix, const int height, const int width) {
	int i, j;
	
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++)
			printf("%3d ", matrix[i][j]);
		printf("\n");
	}
		
}

/**
 * Create a new image (original's copy)
 * Delay: Mirror effect. es el trozo que pone por los bordes de mas (rellenandolo efecto espejo)
 */
int ** getMirror(int** matrix, const int height, const int width, const int delay) {
	int **tmp = malloc(sizeof(int *)*(height+2*delay));
	
	if (tmp == NULL)
		exit_msg("getImage: cannot allocate memory (1)");
	
	int i,j;
	for (i = 0; i < height+2*delay; i++) {
		if ((tmp[i] = malloc(sizeof(int) * (width+2*delay))) == NULL)
			exit_msg("getImage: cannot allocate memory (2)");
	}
	
	//cpy image
	for (i = 0; i < height; i++)
		for (j = 0; j < width; j++)
			tmp[i+delay][j+delay] = matrix[i][j];
	//put mirror up
	for (i = 0; i < delay; i++)
		for (j = 0; j < width; j++)
			tmp[delay-i-1][j+delay] = matrix[i][j];
	//put mirror down
	for (i = 0; i < delay; i++)
		for (j = 0; j < width; j++)
			tmp[height+delay+i][j+delay] = matrix[height-i-1][j];
	//put mirror left
	for (i = 0; i < height+2*delay; i++)
		for (j = 0; j < delay; j++)
			tmp[i][delay-j-1] = tmp[i][j+delay];
	//put mirror right
	for (i = 0; i < height+2*delay; i++)
		for (j = 0; j < delay; j++)
			tmp[i][width+delay+j] = tmp[i][width+delay-j-1];
	
	return tmp;
}

float Ck(const float _k){
	float k = MAX(0, _k);
	
	return k;
}

float Pk (const float k) {
	
	return (1.0f/6.0f) * (powf(Ck(k+2),3.0) - 4 * powf(Ck(k+1),3.0) + 6 * powf(Ck(k),3.0) - 4 * powf(Ck(k-1),3.0));
}

/**
 * Lo que va a hacer es devolver la imagen escalada
 * 
 */
int ** getScale(int **mirror, const int height, const int width, const int delay, const float ix, const float iy) {
	int **tmp = malloc(sizeof(int *)*(height*iy));
	
	if (tmp == NULL)
		exit_msg("getImage: cannot allocate memory (1)");
	
	int i,j;
	for (i = 0; i < height*iy; i++) {
		if ((tmp[i] = malloc(sizeof(int) * (width*ix))) == NULL)
			exit_msg("getImage: cannot allocate memory (2)");
	}

	int cnt = 0;
	int n, m;
	float sum;
	float pn, pm;
	float a, b;
	//Iniciar el proceso de convolucion
	for (i = 0; i < height*iy; i++) {//MPI
		for (j = 0; j < width*ix; j++) {
			//tmp[i][j] = MAX(0, cnt++);
			sum = 0.0f;
			a = ((float) i)/ix  - ((int) i/ix);
			b = ((float) j)/iy  - ((int) j/iy);
			for (n = -1; n < 3; n++) {
				for (m = -1; m < 3; m++) {
					pn = Pk(n - a);
					pm = Pk(b - m);
					sum += ((float) mirror[(int) (i/ix+delay+n)][(int) (j/iy+delay+m)])*pn*pm;
				}
			}
			tmp[i][j] = (int) sum;
			
		}
	}


	return tmp;
}

int main(int argc, char **argv) {
	if (argc != 5) {
		printf("Invalid arguments\n");
		printf("\tprogram height width ix iy\n");
		exit(EXIT_FAILURE);	
	}

	int height = atoi(argv[1]),
		width = atoi(argv[2]);
	float ix = atof(argv[3]), //zoom en x
		  iy = atof(argv[4]); //zoom en y
	int delay = 2;

	struct timeval t0, t1, t;	
	
	int **originalImage = getImage(height, width);	
	assert (gettimeofday (&t0, NULL) == 0);
//	El mirrorImage se cuenta como entrada del procesamiento, asi no contamos el tiempo que se tarda crearlo
	int **mirrorImage   = getMirror(originalImage, height, width, delay);
	int **scaleImage    = getScale(mirrorImage, height, width, delay, ix, iy);
	
	assert (gettimeofday (&t1, NULL) == 0);
	timersub(&t1, &t0, &t);
	
//	printf("Original\n");
//	print_image(originalImage, height, width);

//	printf("Espejo\n");
//	print_image(mirrorImage, height+4, width+4);

//	printf("Resultado\n");
//	print_image(scaleImage, height*iy, width*ix);
	
	printf ("Tiempo      = %ld:%ld(seg:mseg)\n", t.tv_sec, t.tv_usec/1000);

	exit(0);

}
