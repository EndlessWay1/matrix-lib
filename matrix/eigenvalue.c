#include "mat_cplx.h"
#include "mat_cplx.c"

// b32 mat_QR(matrix_cplx* Q, matrix_cplx* R, matrix* A);
// b32 mat_QRc(matrix_cplx* Q, matrix_cplx* R, matrix_cplx* A);
/*
We Know QQ.T have determinant +- 1, thus we construct R

*/

int main (void){
    mem_arena* perm_arena = arena_create(GiB(1), MiB(1));

    matrix* Q = mat_create(perm_arena, 3, 3);
    Q->data[0] = 1/sqrt(2);
    Q->data[1] = 1/sqrt(6);
    Q->data[2] = 1/sqrt(3);
    Q->data[3] = -1/sqrt(2);
    Q->data[4] = 1/sqrt(6);
    Q->data[5] = 1/sqrt(3);
    Q->data[6] = 0;
    Q->data[7] = -2/sqrt(6);
    Q->data[8] = 1/sqrt(3);

    TRACE_TXT("+++++++++++++++++++++++++\n");
    print_mat(Q);
    TRACE_TXT("+++++++++++++++++++++++++\n");
    
    TRACE_TXT("DETERMINANT: %Lf\n", mat_det(Q));
    // mat_ref(Q, NULL, Q, NULL);
    mat_rref(Q,Q,Q);
    TRACE_TXT("+++++++++++++++++++++++++\n");
    print_mat(Q);
    TRACE_TXT("+++++++++++++++++++++++++\n");
    
    matrix* P = mat_create(perm_arena, 3, 4);
    P->data[0] = 1;
    P->data[1] = 3;
    P->data[2] = 0;
    P->data[3] = 0;

    P->data[4] = 0;
    P->data[5] = 0;
    P->data[6] = 1;
    P->data[7] = 2;
    
    P->data[8] = 0;
    P->data[9] = 0;
    P->data[10] = 0;
    P->data[11] = 1;

    TRACE_TXT("+++++++++++++++++++++++++\n");
    print_mat(P);
    TRACE_TXT("+++++++++++++++++++++++++\n");
    mat_rref(P, P, P);
    print_mat(P);
    TRACE_TXT("+++++++++++++++++++++++++\n");

    arena_destroy(perm_arena);

    return 1;
}

// b32 mat_QR(matrix_cplx* Q, matrix_cplx* R, matrix* A){
//     if (mat_rank(A) != A->row) return false;
//     return false;
// }

// b32 mat_QRc(matrix_cplx* Q, matrix_cplx* R, matrix_cplx* A);

