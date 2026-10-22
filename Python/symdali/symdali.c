#define PY_SSIZE_T_CLEAN //You need the usual python syntaxe + some extra calls for numpy
#include <Python.h>
#include <numpy/arrayobject.h>
#include <numpy/ndarraytypes.h>



long long binomial(long long n, long long m) {
    if (m > n) return 0;
    if (m == 0 || m == n) return 1;
    if (m > n - m) m = n - m;  // symmetry
    long long res = 1;
    for (long long i = 1; i <= m; i++) {
        res = res * (n - m + i) / i;
    }
    return res;
}



/* DeltaP_Taylor:
 * DeltaP: array of length dim
 * dim: dimension of DeltaP
 * rank: highest Taylor rank
 * out_size_out: returns the size of the output array
 */

void DeltaP_Taylor(double *DeltaP, double *output, int dim, int Taylor_Order)
{
    /* ----- Compute DeltaP2 ----- */
    long long size_2 = (long long) dim*(dim + 1)/2;
    double* DeltaP2 = (double*) malloc(size_2 * sizeof(double));

    int idx = 0;
    for (int i = 0; i < dim; i++) {
        for (int j = i; j < dim; j++) {
            DeltaP2[idx++] = DeltaP[i]*DeltaP[j];
        }
    }

    //Copy DeltaP2 to output:
    for(int i=0; i<size_2; i++) {
      output[i] = DeltaP2[i];
    }

    //If taylor order > 2 there is more work to be done:
    if(Taylor_Order>2) {
      long long current_idx = size_2; //First idx with no elements in output
      long long start_pr = 0; /* start index of previous rank block */

      for (int j = 2; j < Taylor_Order; j++) {

        long long i_start_pr = start_pr;

        for (int i = 0; i < dim; i++) {
            /*number of LI componenents to multiply DeltaP[i]. It always starts with all LI components 
            of the previous rank and as you move foward it is all LI components of the previous rank in 
            lower dimensions*/
            long long LI = binomial(j + (dim - i) - 1, j);

            for (long long k = 0; k < LI; k++) { //Execute the multiplication.
                output[current_idx++] = DeltaP[i] * output[i_start_pr + k];
            }

            i_start_pr += binomial((j - 1) + (dim - i) - 1, j - 1);
        }

        start_pr += binomial(j + dim - 1, j);
    }
  }
  
  free(DeltaP2);
}



static PyObject *
test_binomial(PyObject *self, PyObject *args) {
  int n, m;

  if (!PyArg_ParseTuple(args, "ii", &n, &m)) 
        return NULL;

  long long res = binomial(n, m); // pass the input to the C-function

  return PyLong_FromLongLong(res); 
}


static PyObject *
Taylor_DeltaPs (PyObject *self, PyObject *args) {
  PyObject *input; 
  npy_intp dims[1]; //To create the output array
  PyArrayObject *in_arr = NULL;
  PyArrayObject *output = NULL;
  int DALI_order;//Derivative approximation order you are working with
  int taylor_order; //the corresponding taylor order

  /* Parse arguments */
  if (!PyArg_ParseTuple(args, "Oi", &input, &DALI_order))
    return NULL;
  taylor_order = 2*DALI_order;
  
  in_arr = (PyArrayObject *) PyArray_FROM_OTF 
    (input, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY); // Converting the array into contiguous

  if (in_arr == NULL)
    goto fail;

  //int ndim = 1; //Number of dimensions of the input vector
  int len = (int) PyArray_SIZE(in_arr); //len of the vector
  //length of the output vector:
  long long out_size = 0;
  for(int rank=2; rank<=taylor_order; rank++) {
    out_size += binomial(len+rank-1, rank);
  }
  dims[0] = (npy_intp) (out_size);
  /* Create output array */
  output = (PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_DOUBLE);
  if (output==NULL) 
    goto fail;

  /* Get pointers to input and output data */ 
  double * in_data = (double *) PyArray_DATA(in_arr);
  double * out_data = (double *) PyArray_DATA(output);
  
  //calculate output vector
  DeltaP_Taylor(in_data, out_data, len, taylor_order);

  /* Decrement reference count of in_arr since it won't be returned */
  Py_XDECREF(in_arr);
  return (PyObject *)output; //Return the output

  fail: 
  /* Decrement reference counts and return NULL */ //IF THE FUNCTION FAILS DELETES INPUT AND OUTPUT ARRAYS and return null
  Py_XDECREF(in_arr);
  Py_XDECREF(output);
  return NULL;
}


// CREATES MODULE, PUTS THE function inside the MODULE and initialize it

static PyMethodDef symdaliMethods[] = {
  {"Test_binomial",  test_binomial, METH_VARARGS,  "test_binomial implementation"},
  {"Taylor_DeltaPs",  Taylor_DeltaPs, METH_VARARGS,  "Taylor_DeltaPs(DeltaP, DALI_order): calculates all DeltaPs to be contracted with the tensors in the approximation of order DALI_order"},
  {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef symdaliModule = {
  PyModuleDef_HEAD_INIT,
  "symdali",   /* name of module */
  "Square elements of arrays", /* module documentation, may be NULL */
  -1,       /* size of per-interpreter state of the module,
	       or -1 if the module keeps state in global variables. */
  symdaliMethods
};

PyMODINIT_FUNC
PyInit_symdali(void)
{
  import_array();
  return PyModule_Create(&symdaliModule);
}
