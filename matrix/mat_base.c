
matrix* mat_create(mem_arena* arena, u32 row, u32 col){
    matrix* mat = PUSH_STRUCT(arena, matrix);
    mat->row = row;
    mat->col = col;
    mat->data = PUSH_ARRAY(arena, f128, (u64)(row)*(col));

    return mat;
}

matrix* mat_load(mem_arena* arena, u32 row, u32 col, const char* filename){
    matrix* mat = mat_create(arena, row, col);

    FILE* f = fopen(filename, "rb");

    fseek(f, 0, SEEK_END);
    u64 size = ftell(f);
    fseek(f, 0, SEEK_SET);

    size = MIN(size, sizeof(f128) * row * col);

    fread(mat->data, 1, size, f);

    fclose(f);

    return mat;
}


void print_mat(const matrix* mat){
    for(u64 i= 0; i < mat->row; i++){
        printf("[ %.12Lf", mat->data[i*mat->col]);
        for(u64 j = 1; j < mat->col; j++){
            printf(", %.12Lf", mat->data[i*mat->col + j]);
        }
        printf("]\n");
    }
}

b32 mat_copy(matrix* dst, const matrix* src){
    if (dst->row != src->row || dst->col != src->col){
        return false;
    }

    memcpy(dst->data, src->data, sizeof(f128)*(u64) dst->row*dst->col);

    return true;
}

void mat_clear(matrix* mat){
    memset(mat->data, 0, sizeof(f128)*(u64) mat->col*mat->row);
}

void mat_fill(matrix* mat, f128 x){
    u64 size = (u64) mat->row * mat->col;

    for (u64 i = 0; i < size; i++){
        mat->data[i] = x;
    }
}

void mat_fill_rand(matrix* mat, f128 lower, f128 upper){
    u64 size = (u64) mat->row * mat->col;
    f128 lower_r = MIN(upper, lower);
    f128 upper_r = MAX(upper, lower); 
    f128 delta = upper_r - lower_r;
    for (u64 i = 0; i < size; i++){
        mat->data[i] = lower_r + (f128)(rand()/ (f128)(RAND_MAX/delta));
    }
}


f128 mat_sum(matrix* mat){
    u64 size = (u64) mat->row * mat->col;
    f128 sum = 0.0f;

    for (u64 i = 0; i < size; i++){
        sum += mat->data[i];
    }
    return sum;
}

void mat_scale(matrix* mat, f128 scale){
    u64 size = (u64) mat->row * mat->col;

    for (u64 i = 0; i < size; i++){
        mat->data[i] *= scale;
    }
}

u64 mat_argmax(matrix* mat){
    u64 size = (u64) mat->col* mat->row;

    u64 max_i = 0;
    for (u64 i = 0; i < size; i++){
        if (mat->data[i] > mat->data[max_i]){
            max_i = i;
        }
    }

    return max_i;

}


b32 mat_add(matrix* out, const matrix* a, const matrix* b){
    if (a->row != b->row || a->col != b->col) return false;
    if (a->row != out->row || a->col != out->col) return false;

    u64 size = (u64) a->row * a->col;

    for (u64 i = 0; i < size; i++){
        out->data[i] = a->data[i] + b->data[i];
    }

    return true;
    
}
b32 mat_sub(matrix* out, const matrix* a, const matrix* b){
    if (a->row != b->row || a->col != b->col) return false;
    if (a->row != out->row || a->col != out->col) return false;


    u64 size = (u64) a->row * a->col;

    for (u64 i = 0; i < size; i++){
        out->data[i] = a->data[i] - b->data[i];
    }

    return true;
}

void _mat_mul_nn(matrix* out, const matrix* a, const matrix* b){
    for (u64 i = 0; i < out->row; i++){
        for (u64 k = 0; k < a->col; k++){
            for (u64 j = 0; j < out->col; j++){
                out->data[i*out->col + j] += 
                a->data[i*a->col + k]*
                b->data[k*out->col + j];
            }
        }
    }
}
void _mat_mul_tn(matrix* out, const matrix* a, const matrix* b){
    for (u64 k = 0; k < a->col; k++){
        for (u64 i = 0; i < out->row; i++){
            for (u64 j = 0; j < out->col; j++){
                out->data[i*out->col + j] += 
                a->data[k*a->col + i]*
                b->data[k*out->col + j];
            }
        }
    }
}
void _mat_mul_nt(matrix* out, const matrix* a, const matrix* b){
    for (u64 i = 0; i < out->row; i++){
        for (u64 j = 0; j < out->col; j++){
            for (u64 k = 0; k < a->col; k++){
                out->data[i*out->col + j] += 
                a->data[i*a->col + k]*
                b->data[j*b->col + k];
            }
        }
    }
}
void _mat_mul_tt(matrix* out, const matrix* a, const matrix* b){
    for (u64 i = 0; i < out->row; i++){
        for (u64 j = 0; j < out->col; j++){
            for (u64 k = 0; k < a->col; k++){
                out->data[i*out->col + j] += 
                a->data[k*a->col + i]*
                b->data[j*out->col + k];
            }
        }
    }
}

b32 mat_mul(
    matrix* out, const matrix* a, const matrix* b, 
    b8 zero_out, b8 transpose_a, b8 transpose_b
){
    mem_arena_temp scratch = arena_scratch_get(NULL, 0);
    matrix* copy_a = mat_create(scratch.arena, a->row, a->col);
    matrix* copy_b = mat_create(scratch.arena, b->row, b->col);

    mat_copy(copy_a, a);
    mat_copy(copy_b, b);
    a = copy_a;
    b = copy_b;

    u32 a_row = transpose_a ? a->col : a->row;
    u32 a_col = transpose_a ? a->row : a->col;
    u32 b_row = transpose_b ? b->col : b->row;
    u32 b_col = transpose_b ? b->row : b->col;
    if (a_col != b_row || out->row != a_row || out->col != b_col) return false;
    

    if (zero_out){
        mat_clear(out);
    }

    u32 transpose = (transpose_a << 1) | transpose_b;

    switch (transpose){
        case 0b00: { _mat_mul_nn(out, a, b); } break; 
        case 0b10: { _mat_mul_tn(out, a, b); } break; 
        case 0b01: { _mat_mul_nt(out, a, b); } break; 
        case 0b11: { _mat_mul_tt(out, a, b); } break; 
    }

    arena_scratch_release(scratch);
    
    return true;
}
