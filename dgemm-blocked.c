/*
    Please include compiler name below (you may also include any other modules you would like to be loaded)

COMPILER= gnu

    Please include All compiler flags and libraries as you want them run. You can simply copy this over from the Makefile's first few lines

CC = cc
OPT = -O3
CFLAGS = -Wall -std=gnu99 $(OPT)
MKLROOT = /opt/intel/composer_xe_2013.1.117/mkl
LDLIBS = -lrt -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm

*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>


const char* dgemm_desc = "My Simple blocked dgemm.";

#if !defined(BLOCK_SIZE)
#define BLOCK_SIZE 41
#endif

#define min(a,b) (((a)<(b))?(a):(b))
#define calculateOffset(row, col, numRow) (col * numRow + row)

//utility functions
void pack(double* dest, double* src, int offset, int width, int height, int lda){
  int numBytes = height * sizeof(double);
  double* pointer_src = src + offset;
  double* pointer_dest = dest;
  for (int j = 0; j < width; j++){
    memcpy(pointer_dest, pointer_src, numBytes);
    pointer_dest  += height;
    pointer_src   += lda;
  }
}

/*
Each block
A = m_c * k_c
B = k_c * n_r
C = m_c * n_r
*/
void do_block(int lda, int m_c, int n_r, int k_c, double* packed_A, double* packed_B, double* C){
  for (int i = 0; i < m_c; ++i){
      for (int j = 0; j < n_r; ++j){
        double val = C[calculateOffset(i, j, lda)];
        for (int p = 0; p < k_c; ++p){
          val += packed_A[calculateOffset(i, p, m_c)] * packed_B[calculateOffset(p, j, k_c)];
        }
        C[calculateOffset(i, j, lda)] = val;
      }
  }
}

void GEBP(int lda, int i, int p, double* A, double* B, double* C, double* packed_A, double* packed_B, int k_c, int m_c){
  //pack A into packed_A
  pack(packed_A, A, calculateOffset(i, p, lda), k_c, m_c, lda);

  for (int j = 0; j < lda; j += BLOCK_SIZE){
    int n_r = min(BLOCK_SIZE, lda - j);
    //call work horse here
    do_block(lda, m_c, n_r, k_c,
      packed_A,
      packed_B + calculateOffset(0, j, k_c),
      C + calculateOffset(i, j, lda));
  }
}

void GEPP(int lda, int p, double* A, double* B, double* C, double* packed_A, double* packed_B, int k_c){
  //pack B into packed_B
  pack(packed_B, B, calculateOffset(p, 0, lda), lda, k_c, lda);

  //break each column of A into multiple blocks
  for (int i = 0; i < lda; i += BLOCK_SIZE){
    int m_c = min(BLOCK_SIZE, lda - i);
    //call GEBP here
    GEBP(lda, i, p, A, B, C, packed_A, packed_B, k_c, m_c);
  }
}


void square_dgemm(int lda, double* A, double* B, double* C){
  //prealocaate memory here
  double* packed_A = malloc((BLOCK_SIZE * BLOCK_SIZE + 1) * sizeof(double));
  double* packed_B = malloc((BLOCK_SIZE * lda + 1) * sizeof(double));

  //break B into multiple rows of size: k_c * lda
  //break A into multiple cols of size: lda * k_c
  for (int p = 0; p < lda; p += BLOCK_SIZE){
    int k_c = min(BLOCK_SIZE, lda - p);
    //call GEPP here
    GEPP(lda, p, A, B, C, packed_A, packed_B, k_c);
  }

  free(packed_A);
  free(packed_B);
}





// #include <stdlib.h>
// #include <string.h>
// #include <stdio.h>
//
//
// const char* dgemm_desc = "Simple blocked dgemm.";
//
// #if !defined(BLOCK_SIZE)
// #define BLOCK_SIZE 41
// #endif
//
// #define min(a,b) (((a)<(b))?(a):(b))
// #define calculateOffset(row, col, numRow) (col * numRow + row)
//
// //utility functions
// double* pack(double* src, int offset, int width, int height, int lda){
//   //printf("offset %d \n", offset);
//   //malloc the length of need to copied arr
//   double* dest = malloc((width * height + 1) * sizeof(double));
//
//   int numBytes = height * sizeof(double);
//   double* pointer_src = src + offset;
//   double* pointer_dest = dest;
//   for (int j = 0; j < width; j++){
//     memcpy(pointer_dest, pointer_src, numBytes);
//     pointer_dest  += height;
//     pointer_src   += lda;
//   }
//   return dest;
// }
//
//
// /*
// Each block
// A = m_c * k_c
// B = k_c * n_r
// C = m_c * n_r
// */
// void do_block(int lda, int m_c, int n_r, int k_c, double* packed_A, double* packed_B, double* C){
//   for (int i = 0; i < m_c; i++){
//       for (int j = 0; j < n_r; j++){
//         double val = C[calculateOffset(i, j, lda)];
//         for (int p = 0; p < k_c; p++){
//           val += packed_A[calculateOffset(i, p, m_c)] * packed_B[calculateOffset(p, j, k_c)];
//           //val += packed_A[calculateOffset(i, p, lda)] * packed_B[calculateOffset(p, j, k_c)];
//         }
//         C[calculateOffset(i, j, lda)] = val;
//       }
//   }
// }
//
// void GEBP(int lda, int i, int p, double* A, double* packed_B, double* C, int k_c, int m_c){
//   //pack A into packed_A
//   double* packed_A = pack(A, calculateOffset(i, p, lda), k_c, m_c, lda);
//
//   for (int j = 0; j < lda; j += BLOCK_SIZE){
//     int n_r = min(BLOCK_SIZE, lda - j);
//     //printf("lda %d, m_c %d, n_r %d, k_c %d \n", lda, m_c, n_r, k_c);
//     //call work horse here
//     do_block(lda, m_c, n_r, k_c,
//       packed_A,
//       //A + calculateOffset(i, p, lda),
//       packed_B + calculateOffset(0, j, k_c),
//       C + calculateOffset(i, j, lda));
//   }
//   //free packed_A
//   free(packed_A);
// }
//
// void GEPP(int lda, int p, double* A, double* B, double* C, int k_c){
//   //pack B into packed_B
//   double* packed_B = pack(B, calculateOffset(p, 0, lda), lda, k_c, lda);
//
//   //break each column of A into multiple blocks
//   for (int i = 0; i < lda; i += BLOCK_SIZE){
//     int m_c = min(BLOCK_SIZE, lda - i);
//     //call GEBP here
//     GEBP(lda, i, p, A, packed_B, C, k_c, m_c);
//   }
//   //free packed_B
//   free(packed_B);
// }
//
//
// void square_dgemm(int lda, double* A, double* B, double* C){
//   //break B into multiple rows of size: k_c * lda
//   //break A into multiple cols of size: lda * k_c
//   for (int p = 0; p < lda; p += BLOCK_SIZE){
//     int k_c = min(BLOCK_SIZE, lda - p);
//     //call GEPP here
//     GEPP(lda, p, A, B, C, k_c);
//   }
// }










// /*
//     Please include compiler name below (you may also include any other modules you would like to be loaded)
//
// COMPILER= gnu
//
//     Please include All compiler flags and libraries as you want them run. You can simply copy this over from the Makefile's first few lines
//
// CC = cc
// OPT = -O3
// CFLAGS = -Wall -std=gnu99 $(OPT)
// MKLROOT = /opt/intel/composer_xe_2013.1.117/mkl
// LDLIBS = -lrt -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm
//
// */
//
// #include <stdlib.h>
// #include <string.h>
//
// const char* dgemm_desc = "Simple blocked dgemm.";
//
// #if !defined(BLOCK_SIZE)
// #define BLOCK_SIZE 41
// #endif
//
// #define min(a,b) (((a)<(b))?(a):(b))
// #define calculateOffset(row, col, col_length) (col * col_length + row)
//
// //utility functions
// double* pack(double* src, int lda, int width, int height){
//   //malloc the length of need to copied arr
//   double* dest = (double*) malloc(width * height * sizeof(double));
//
//   int height_in_bytes = width * sizeof(double);
//   for (int j = 0; j < width; j++){
//     memcpy(dest, src+j, height_in_bytes);
//   }
//   return packed_arr;
// }
//
//
// /*
// Each block
// A = m_c * k_c
// B = k_c * n_r
// C = m_c * n_r
// */
// void do_block(int lda, int m_c, int n_r, int k_c, double* A, double* B, double* C){
//   for (int i = 0; i < m_c; i++){
//       for (int j = 0; j < n_r; j++){
//
//         double val = C[calculateOffset(i, j, lda)];
//
//         for (int p = 0; p < k_c; p++){
//
//           val += A[calculateOffset(i, p, lda)] * B[calculateOffset(p, j, lda)];
//         }
//
//         C[calculateOffset(i, j, lda)] = val;
//       }
//   }
// }
//
// void GEBP(int lda, int i, int p, double* A, double* B, double* C, int k_c, int m_c){
//   for (int j = 0; j < lda; j += BLOCK_SIZE){
//     //calculate the small block here
//     int n_r = min(BLOCK_SIZE, lda - j);
//     do_block(lda, m_c, n_r, k_c,
//       A + calculateOffset(i, p, lda),
//       B + calculateOffset(p, j, lda),
//       C + calculateOffset(i, j, lda));
//   }
// }
//
// void GEPP(int lda, int p, double* A, double* B, double* C, int k_c){
//   for (int i = 0; i < lda; i += BLOCK_SIZE){
//     //call GEBP here
//     int m_c = min(BLOCK_SIZE, lda - i);
//     GEBP(lda, i, p, A, B, C, k_c, m_c);
//   }
// }
//
// void square_dgemm(int lda, double* A, double* B, double* C){
//   //break B into multiple rows of size: k_c * lda
//   for (int p = 0; p < lda; p += BLOCK_SIZE){
//     //call GEPP here
//     int k_c = min(BLOCK_SIZE, lda - p);
//     GEPP(lda, p, A, B, C, k_c);
//   }
// }
