%module sp
%{
struct sparse_matrix_t;
%}

extern struct sparse_matrix_t*
sp_load (const char* path, const char* fmt);

extern int
sp_save (struct sparse_matrix_t* A, const char* path, 
      const char* fmt);

extern void
sp_format (struct sparse_matrix_t* A);

extern int
sp_convert (struct sparse_matrix_t* A, const char* type);

extern struct sparse_matrix_t* 
sp_mult (struct sparse_matrix_t* B, struct sparse_matrix_t* A);

