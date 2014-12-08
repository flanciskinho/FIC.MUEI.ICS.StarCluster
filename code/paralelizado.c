#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

//#include <omp.h>

#include <mpi.h> // MPI header file

#define MAX(x,y) ((x>y)?x:y)
#define MIN(x,y) ((x<y)?x:y)

#define INDEX(i,j,ncols) [(i)*(ncols) + (j)]

/**
 * Print a message before exit
 */
void exit_msg(char *msg, int nproc) {
	printf("\tPROC: %d > %s\n", nproc, msg);
	fflush(stdout);
	MPI_Finalize();
	exit(EXIT_FAILURE);
}

/**
 * Simulate the reading of an image
 */
void getImage(short *matrix, int height, int width) {
	int i,j;
	
	for (i = 0; i < height; i++)
		for (j = 0; j < width; j++)
			matrix INDEX(i, j, width) = rand() % 256;//gray level
}

/**
 * Create a new image (original's copy)
 * Delay: Mirror effect. es el trozo que pone por los bordes de mas (rellenandolo efecto espejo)
 */
void getMirror(short* matrix, short *mirror, int height, int width, int delay) {
	int i, j;
	
	int depth = width+2*delay;
	//cpy image
#pragma omp parallel private(i, j)
{
#pragma omp for
	for (i = 0; i < height; i++)
		for (j = 0; j < width; j++)
		mirror INDEX(i+delay, j+delay, depth) = matrix INDEX(i, j, width);
}

	//put mirror up
	for (i = 0; i < delay; i++)
		for (j = 0; j < width; j++)
			mirror INDEX(delay-i-1, j+delay, depth) = matrix INDEX(i, j, width);

	//put mirror down
	for (i = 0; i < delay; i++)
		for (j = 0; j < width; j++)
			mirror INDEX(height+delay+i, j+delay, depth) = matrix INDEX((height-i-1), j, width);

	//put mirror left
	for (i = 0; i < height+2*delay; i++)
		for (j = 0; j < delay; j++)
			mirror INDEX(i, delay-j-1, depth) = mirror INDEX(i, j+delay, depth);

	//put mirror right
	for (i = 0; i < height+2*delay; i++)
		for (j = 0; j < delay; j++)
			mirror INDEX(i, width+delay+j, depth) = mirror INDEX(i, width+delay-j-1, depth);

}

float Ck(float _k){
	float k = MAX(0, _k);
	
	return k;
}

float Pk (float k) {	
	return (1.0f/6.0f) * (powf(Ck(k+2),3.0) - 4 * powf(Ck(k+1),3.0) + 6 * powf(Ck(k),3.0) - 4 * powf(Ck(k-1),3.0));
}


void getScale(short *mirror, short *result, int height, int width, int delay, float ix, float iy) {
	int i, j;
	int n, m;
	float a, b;
	float pn, pm;
	int size = height*iy;
	int depth = width+2*delay;

	float sum;

#pragma omp parallel private(i, j, sum, a, b, n, m, pn, pm)
{
#pragma omp for
	for (i = 0; i < size; i++) {//MPI
		for (j = 0; j < width*ix; j++) {//OMP
			sum = 0.0f;
			a = ((float) i)/ix  - ((int) i/ix);//solo cogemos la parte decimal
			b = ((float) j)/iy  - ((int) j/iy);
			for (n = -1; n < 3; n++) { //Ambos bucles con SSE/AVX
				for (m = -1; m < 3; m++) {
					pn = Pk(n - a);
					pm = Pk(b - m);
					sum += ((float) mirror INDEX ((int) (i/ix+delay+n), (int) (j/iy+delay+m), depth))*pn*pm;
				}
			}
			result INDEX(i, j, (int) (width*ix)) = (int) sum;
			
		}
	}

}
	

}


/**
 * Imprime la matriz
 */
void print_image(short *matrix, int height, int width) {
	int i, j;
	
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++)
			printf("%3d ", matrix INDEX(i, j, width));
		printf("\n");
	}
}

/**************************************************************
 ** Funciones para inicializar variables para poder usar MPI **
 **************************************************************/
//Función para indicar las filas con las que trabajará cada proceso
void fill_aux_rows(int *array, int size, int height, int delay) {
	int rows = height / size;	
	if (height % size != 0)
		rows++;
	
	int i;
	for (i = 0; i < size; i++) {
		array[i] = MIN(rows, height-i*rows) + 2*delay;
//printf("rows(%d) = %d\n", i,array[i]);
	}
}

//Función para conocer el número de elementos a enviar a cada proceso
void fill_size_send_child(int *fill, int size, int *rows, int width, int delay) {
	int i;
	
	for (i = 0; i < size; i++) {
		fill[i] = rows[i]* (width + 2*delay);
//printf("size(%d) = %d\n", i,fill[i]);
	}

}

//Función para conocer desde donde se empieza a enviar a cada proceso
void fill_aux_split(int *fill, int size, int stride, int width, int delay) {
	int i;
	
	for (i = 0; i < size; i++) {
		fill[i] = i * (stride - (width + delay*2)*4);
//printf("begin(%d) = %d\n", i, fill[i]);
	}

}

//Funcion para conocer el tamaño de la respuesta que enviara cada proceso
void fill_size_send_dady(int *fill, int size, int *rows, int width, float ix, float iy, int delay) {
	int i;
	
	for (i = 0; i < size; i++) {
		fill[i] = (rows[i] -(2*delay))*iy * (width*ix);
//printf("SIZE(%d) = %d\n", i, fill[i]);
	}
}

//Funcion que indica al proceso que almacena el resultado donde tiene que colocar los resultados parciales
void fill_aux_begin(int *fill, int size, int bytes) {
	int i;
	
	for (i = 0; i < size; i++) {
		fill[i] = bytes*i;
//printf("BEGIN(%d) = %d\n", i, fill[i]);
	}
}
/**************************************************************
 **             Fin de Funciones para MPI                    **
 **************************************************************/
void
diff_time(struct timeval *result, struct timeval *a, struct timeval *b)
{
	int oflg, usec, sec;

	oflg = 0;

	usec = a->tv_usec - b->tv_usec;

	if (usec < 0) {
		usec += 1000000;
		oflg = 1;
	}

	sec = a->tv_sec - b->tv_sec - oflg;

	if (a->tv_sec < b->tv_sec)
	sec += 60 * 60 * 24;

	result->tv_usec = usec;
	result->tv_sec = sec;
}

int main(int argc, char **argv) {
//	_MM_ROUND_NEAREST
//	_MM_ROUND_DOWN
//	_MM_ROUND_UP
//	_MM_ROUND_TOWARD_ZERO
//	_MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);
	
	if (argc != 5) {
		printf("Invalid arguments\n");
		printf("\tprogram height width ix iy\n");
		exit(EXIT_FAILURE);	
	}

	struct timeval t, t2;
	int height = atoi(argv[1]),
		width = atoi(argv[2]);
	float ix = atof(argv[3]), //zoom en x
		  iy = atof(argv[4]); //zoom en y
	int delay = 2;
	
	short *mirror = NULL; //donde dira el padre que empieza el espejo global
	short *scale = NULL; //donde dira el padre que hay que guardar las cosas
	
	//Las dos siguientes las modifica MPI
	int nprocs, myid;

	MPI_Init(&argc, &argv);//inicializamos mpi
	
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	int aux_rows[nprocs];//Las filas con las que trabajara cada proceso
	int size_send_child[nprocs];//El numero de elementos a enviar a cada hijo
	int aux_split[nprocs];//Es para saber donde empezar a contar (este es para enviar a los hijos)
	int size_send_dad[nprocs]; //para saber el tamano que se le va a enviar al padre
	int aux_begin[nprocs]; //Es el para saber donde empezar a contar (este es para enviar al padre)
	
	//Con esto indicamos que todos saben lo que tienen todos
	fill_aux_rows(aux_rows, nprocs, height, delay);
	fill_size_send_child(size_send_child, nprocs, aux_rows, width, delay);
	fill_aux_split(aux_split, nprocs, size_send_child[0], width, delay);
	fill_size_send_dady(size_send_dad, nprocs, aux_rows, width, ix, iy, delay);
	fill_aux_begin(aux_begin, nprocs, size_send_dad[0]);

	if (myid == 0) {//Parent process
		short *originalImage = malloc(sizeof(short) * height * width); //Imagen original
		mirror = malloc(sizeof(short) * (height+delay*2) * (width+delay*2)); //Imagen en espejo
		if (originalImage == NULL || mirror == NULL) {
			exit_msg("main: cannot allocate memory", myid);	
		}
		
		getImage(originalImage, height, width);
		getMirror(originalImage, mirror, height, width, delay);
		free(originalImage);

		gettimeofday(&t, NULL);//Empezamos a contar el tiempo
	}

	short *work = malloc(sizeof(short) * (aux_rows[myid]) * (width+2*delay) ); //con la que trabajara cada proceso
	if (work == NULL)
		exit_msg("main: cannot allocate memory (for work)", myid);

	//Enviamos a los procesos la informacion
	if (0 != MPI_Scatterv(mirror, size_send_child, aux_split, MPI_SHORT,
			work, size_send_child[myid], MPI_SHORT,
			0, MPI_COMM_WORLD))
		exit_msg("main: MPI_Scatterv", myid);


	short *result = malloc(sizeof(short) * ((int) ((aux_rows[myid]-2*delay)*iy)) * (int) (width*ix));
	if (result == NULL)
		exit_msg("main: cannot allocate memory (for result)", myid);

	getScale(work, result, aux_rows[myid]-2*delay, width, delay, ix, iy);

	free(work);
	if (myid == 0) {
		free(mirror);
		scale  = malloc(sizeof(short) * ((int) (height*iy)) * ((int) (width*ix))); // donde se guardara el resultado final
	}

	//Los procesos envian la informacion al padre
	if (0 != MPI_Gatherv(result, size_send_dad[myid], MPI_SHORT,
 				scale, size_send_dad, aux_begin,
 				MPI_SHORT, 0, MPI_COMM_WORLD))
		exit_msg("main: MPI_Gatherv", myid);
	
	if (myid == 0) {
		gettimeofday(&t2, NULL);

		struct timeval diff;
		diff_time(&diff, &t2, &t);
//		print_image(scale, height*iy, width*ix);
		printf("Tiempo      = %ld:%ld(seg:mseg)\n\n", diff.tv_sec, diff.tv_usec/1000 );
	}

	MPI_Finalize(); //Clean UP for MPI

	exit(EXIT_SUCCESS);

}
