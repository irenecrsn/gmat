int gram_schmidt_sel (double *mort, int *madj, double *mcov, 
		unsigned int *dim); 

int gram_schmidt (double **span_ort, double **span, 
		unsigned int *nvec, unsigned int *dim); 

int proj_ort (double *v_proj_u, double *v, double *u, 
		unsigned int *dim);
