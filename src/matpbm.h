/**
  * MatPBM - simple I/O utility to read/write bitmaps into sparse matrices.
  *
  * MatPBM uses Compressed Row Storage to load sparse matrices in memory,
  * and matrices can be converted from and to binary PBM.
  *
  * Copyright 2018 Micky Faas <micky@edukitty.org>
  */

#ifndef MATPBM_H
#define MATPBM_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#ifndef MIN
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#endif

typedef unsigned int idx_t;

typedef struct {
    idx_t p, q; // First row and column of the submatrix
    idx_t s; // Width of the border
} submat_t;

typedef struct {
    idx_t* rowptr;      // rowptr[i] gives the start of the i-th row in the colind array
    idx_t* rowlen;      // rowlen[i] gives the number of non-zero elements in row i
    idx_t* colind;      // colind[rowptr[i]+j] gives the index of the j-th column in the i-th row
    size_t m;           // Number of rows/columns (matrix is square)
    size_t nz;          // Number of non-zeroes (size of the colind array)

    idx_t* roworder;    // Order of the rows, roworder[0] givens number of the row etc.
    idx_t* colorder;    // Order of the columns

    submat_t active;
} bmat_t;

bmat_t*
matpbm_loadFromStream( FILE* f );

bool
matpbm_writeToStream( FILE* f, bmat_t* mat );

void
matpbm_printDense( bmat_t* mat );

bool
matpbm_isNZ( bmat_t* mat, idx_t row, idx_t col );

idx_t
matpbm_minNZPerRow( bmat_t* mat );

idx_t
matpbm_rowNNZ( bmat_t* mat, idx_t row, idx_t q, idx_t s );

void
matpbm_swapRows( bmat_t* mat, idx_t row1, idx_t row2 );

void
matpbm_swapCols( bmat_t* mat, idx_t col1, idx_t col2 );

void
matpbm_moveRow( bmat_t* mat, idx_t from, idx_t to );

void
matpbm_moveCol( bmat_t* mat, idx_t from, idx_t to );

void
matpbm_free( bmat_t* mat );

#endif
