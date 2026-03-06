#include "../base.h"
#include <stdarg.h>

// DEBUGING MACRO FROM numpy/linalg/umath_linalg.cpp
#define TRACE_TXT(...) do { fprintf (stderr, __VA_ARGS__); } while (0)

typedef union {
    struct {f128 r, i;};
    f128 v[2];
} complex_data;

typedef struct{
    u32 row;
    u32 col;
    complex_data *data;
} matrix_cplx;

matrix_cplx *mat_c_create(mem_arena* arena, u32 row, u32 col);

b32 mat_r_to_c(matrix_cplx* dst, const matrix* src);
b32 mat_c_to_r(matrix* dst, const matrix_cplx* src);


void print_matc(const matrix_cplx* mat);
void mat_clearc(matrix_cplx* mat);
b32 mat_copyc(matrix_cplx* dst, const matrix_cplx* src);
void mat_fillc(matrix_cplx* mat, f128 r, f128 i);
void mat_fillc_rand(
    matrix_cplx* mat, f128 lower_r, f128 upper_r, 
    f128 lower_i, f128 upper_i);

u64 mat_argmaxc(matrix_cplx* mat, b8 use_i);

complex_data mat_sumc(matrix_cplx* mat);
b32 mat_addc(matrix_cplx* out, const matrix_cplx* a, const matrix_cplx* b);
b32 mat_subc(matrix_cplx* out, const matrix_cplx* a, const matrix_cplx* b);
void mat_scalec(matrix_cplx* mat, f128 scale);
b32 mat_mulc(
    matrix_cplx* out, const matrix_cplx* a, const matrix_cplx* b, 
    b8 zero_out, b8 transpose_a, b8 transpose_b
);