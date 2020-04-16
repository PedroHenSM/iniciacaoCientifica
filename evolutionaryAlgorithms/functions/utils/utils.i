/* cec20.i */
%module utils
%include "carrays.i" /* %array_functions(type,name) */
%array_functions (double, doubleArray);
%array_functions (long double, longDoubleArray);
/* Creates an void array | array_functions(type,name) */
/* %array_class(void,voidArray) Creates an void array for classes | array_class(type,name) */
/* %typemap(in) void* = double*; */

%{

/* Helper functions */

double **new_doubleddArray(int rows){
	double **arr = new double *[rows];
	return arr;
}

double **castToDouble(void *b){
	return (double**)b;
}

double **new_doubleddArray(int rows, int cols) {
    int i;
    double **arr = new double *[rows];
    for (i=0; i<rows; i++)
		arr[i] = new double[cols];
    return arr;
}

void delete_doubleddArray (double **arr, int rows, int cols){
	int i;
	for (i=0; i<rows; i++)
		delete[] arr[i];
	delete[] arr;
}

void doubleddArray_setitem(double **array, int row, int col, double value) {
    array[row][col] = value;
}

double doubleddArray_getitem(double **array, int row, int col) {
    return array[row][col];
}

void printDoubleArray(double *array, int size){
    int i;
    for (i = 0; i < size; i++){
        printf("%f ", array[i]);
    }
    printf("\n");
}

void printLongDoubleArray(long double *array, int size){
    int i;
    for (i = 0; i < size; i++){
        printf("%Lf ", array[i]);
    }
    printf("\n");
}

%}

double **new_doubleddArray(int rows);

double **castToDouble(void *b);

void delete_doubleddArray (double **arr, int rows, int cols);

void doubleddArray_setitem(double **array, int row, int col, double value);

double doubleddArray_getitem(double **array, int row, int col);

void printDoubleArray(double *array, int size);

void printLongDoubleArray(long double *array, int size);
