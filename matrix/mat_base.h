typedef struct {
    u32 row, col;
    f128* data;
} matrix;

matrix* mat_create(mem_arena* arena, u32 row, u32 col);
matrix* mat_load(mem_arena* arena, u32 row, u32 col, const char* filename);

void print_mat(const matrix* mat);
void mat_clear(matrix* mat);
b32 mat_copy(matrix* dst, const matrix* src);
void mat_fill(matrix* mat, f128 x);
void mat_fill_rand(matrix* mat, f128 lower, f128 upper);

u64 mat_argmax(matrix* mat);

f128 mat_sum(matrix* mat);
b32 mat_add(matrix* out, const matrix* a, const matrix* b);
b32 mat_sub(matrix* out, const matrix* a, const matrix* b);
void mat_scale(matrix* mat, f128 scale);
b32 mat_mul(
    matrix* out, const matrix* a, const matrix* b, 
    b8 zero_out, b8 transpose_a, b8 transpose_b
);

