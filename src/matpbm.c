/**
  * MatPBM - simple I/O utility to read/write bitmaps into sparse matrices.
  *
  * MatPBM uses Compressed Row Storage to load sparse matrices in memory,
  * and matrices can be converted from and to binary PBM.
  *
  * Copyright 2018 Micky Faas <micky@edukitty.org>
  */
#include "matpbm.h"
#include <stdbool.h>
#include <string.h>
#include <limits.h>


#define MAX_LINE 4096 // Maximum length of the line buffer
#define BLOCK_SIZE 256 // Allocate memory in blocks of this size
#define PBM_HEADER_COMMENT "CREATOR: Rareman matrix transform demo, https://github.com/mickymuis/rareman"

bmat_t*
matpbm_loadFromStream( FILE* f ) {
    if( f == NULL ) return NULL;

    char line[MAX_LINE];
    bool has_magic = false;
    int m,n;

    /* Read the header while skipping the comments */
    while( 1 ) {
        if( fgets( line, MAX_LINE, f ) == NULL ) {
            fprintf( stderr, "(e) Unexpected EOF while parsing header.\n" );
            return NULL;
        }
        if( line[0] == '#' ) continue;
        if( strcmp( line, "P4\n" ) == 0 )
            has_magic =true;
        else if( has_magic ) {
            sscanf( line, "%d %d\n", &m, &n );  
            break;
        }
    }
    if( !has_magic ) {
        fprintf( stderr, "(e) Unsupported file format (expected binary PBM).\n" );
        return NULL;
    } else if( m != n ) {
        fprintf( stderr, "(e) Image/matrix is not square.\n" );
        return NULL;
    }

    /* Allocate the space for the matrix, only mat->colind is undetermined until we've read all data */
    const int b = m / 8 + ( m % 8 != 0 ? 1 : 0 ); // Number of 8 bit blocks
    bmat_t* mat = malloc( sizeof( bmat_t ) );
    memset( mat, 0, sizeof( bmat_t ) );
    mat->rowptr = malloc( sizeof( idx_t ) * m );
    mat->rowlen = malloc( sizeof( idx_t ) * m );
    mat->colind = malloc( sizeof( idx_t ) * BLOCK_SIZE );
    size_t capacity = BLOCK_SIZE;
    mat->m = m;

    idx_t pos = 0; // Pointer to the current element in colind

    /* Read the binary data in blocks of a byte */
    for( int i =0; i < m; i++ ) {
        mat->rowptr[i] = pos; // Begin of this row
        idx_t len = 0; // Length of the row 

        for( int bj=0; bj < b; bj++ ) {
            int block = fgetc( f );
            //printf( "%x ", (char)block );
            if( block == EOF ) {
                fprintf( stderr, "(e) Unexpected EOF.\n" );
                goto err;
            }
            int mask = 1<<7;
            for( int j = bj*8; j < n && j < bj*8+8; j++ ) {
                bool x = block & mask;
                //printf( "%c ", x ? 'W' : ' ' );
                if( x ) { // Non-zero entry
                    mat->colind[pos++] = j;
                    if( pos == capacity ) { 
                        capacity += BLOCK_SIZE;
                        mat->colind = realloc( mat->colind, sizeof( idx_t ) * capacity );
                        if( !mat->colind ) {
                            fprintf( stderr, "(e) Unable to allocate more memory.\n" );
                            goto err;
                        }
                    }
                    len++;
                }
                mask = mask >> 1;
            }
        }
        mat->rowlen[i] = len;
        //printf( "\n" );
    }

    /* Initialize the trivial reordering arrays */
    mat->roworder = malloc( sizeof( idx_t ) * m );
    mat->colorder = malloc( sizeof( idx_t ) * m );
    for( idx_t i =0; i < m ; i++ ) {
        mat->roworder[i] = i;
        mat->colorder[i] = i;
    }
    mat->nz =pos;

    fprintf( stderr, "(i) Loaded matpbm matrix with density %.5f\n", (double)pos / (double)(m*m) );

    return mat;
err:
    matpbm_free( mat );
    return NULL;
}

bool
matpbm_writeToStream( FILE* f, bmat_t* mat ) {
    if( !f ) return false;

    // Let's write a simple PBM header
    fprintf( f, "P4\n# %s\n%ld %ld\n", PBM_HEADER_COMMENT, mat->m, mat->m );

    const int blockLen =8;
    int blockPos =0;
    char block =0;

    // Iterate over all logical rows
    for( int ii =0; ii < mat->m; ii++ ) {
        // Convert logical to physical rows number
        idx_t i = mat->roworder[ii];
        idx_t ptr = mat->rowptr[i];

        blockPos =0; block =0;
        
        // Iterate over all logical columns
        for( int jj=0; jj < mat->m; jj++ ) {
            // Convert logical to physical column number
            idx_t j =mat->colorder[jj];
    
            // Unfortunately, due to reordering, the colind[] array is not sorted.
            // We have to do a linear search to know if this column is non-zero.
            for( int k =ptr; k < ptr+ mat->rowlen[i]; k++ ) {
                if( mat->colind[k] == j ) { 
                    block |= (1 << (blockLen-1)) >> blockPos;
                    break;
                }
            }
            // Completed one block, write it to disk
            if( ++blockPos == blockLen ) {
                fprintf( f, "%c", block );
                block =0;
                blockPos =0;
            }

        }
        // Write the last (partial) block
        if( blockPos != 0 ) {
            fprintf( f, "%c", block );
        }

    }

    return true;
}

void
matpbm_printDense( bmat_t* mat ) {
    for( int ii =0; ii < mat->m; ii++ ) {
        idx_t i = mat->roworder[ii];
        idx_t ptr = mat->rowptr[i];
        char c = '+';
        for( int jj=0; jj < mat->m; jj++ ) {
            if( jj >= mat->active.p && ii >= mat->active.p )
                c = 'X';
            if( jj >= mat->m - mat->active.q )
                c = 'B';

            idx_t j =mat->colorder[jj];
            bool x =false;
            for( int k =ptr; k < ptr+ mat->rowlen[i]; k++ )
                if( mat->colind[k] == j ) x =true;
            printf( "%c ", x ? c : ' ' );
        }

        printf( "\n" );
    }
}

bool
matpbm_isNZ( bmat_t* mat, idx_t row, idx_t col ) {
    idx_t i = mat->roworder[row];
    idx_t ptr = mat->rowptr[i];
    idx_t j =mat->colorder[col];
    for( int k =ptr; k < ptr + mat->rowlen[i]; k++ )
        if( mat->colind[k] == j ) return true;
    return false; // No non-zero at (row,col)
}

idx_t
matpbm_minNZPerRow( bmat_t* mat ) {
    idx_t min = INT_MAX;
    for( int i =0; i < mat->m; i++ ) {
        min = MIN( min, mat->rowlen[i] );
    }
    return min;
}

idx_t
matpbm_rowNNZ( bmat_t* mat, idx_t row, idx_t q, idx_t s ) {
    idx_t i = mat->roworder[row];
    idx_t nnz = mat->rowlen[i];
    idx_t base = mat->rowptr[i];
    for( int j =0; j < mat->rowlen[i]; j++ ) {
        if( mat->colind[base+j] < q ) nnz--;
        else if( mat->colind[base+j] >= s ) nnz--;
    }
    return nnz;
}

void
matpbm_swapRows( bmat_t* mat, idx_t row1, idx_t row2 ) {
    idx_t tmp = mat->roworder[row1];
    mat->roworder[row1] = mat->roworder[row2];
    mat->roworder[row2] = tmp;
}

void
matpbm_swapCols( bmat_t* mat, idx_t col1, idx_t col2 ) {
    idx_t tmp = mat->colorder[col1];
    mat->roworder[col1] = mat->roworder[col2];
    mat->roworder[col2] = tmp;
}

void
arrayMove( idx_t* A, idx_t from, idx_t to ) {
    if( from == to ) return;

    idx_t value = A[from];

    if( from < to ) {
        memmove( A + from, A + from + 1, sizeof( idx_t ) * (to - from) );    
    } else {
        memmove( A + to + 1, A + to, sizeof( idx_t ) * (from - to) );    
    }
    A[to] = value;
}

void
matpbm_moveRow( bmat_t* mat, idx_t from, idx_t to ) {
    arrayMove( mat->roworder, from, to );
}

void
matpbm_moveCol( bmat_t* mat, idx_t from, idx_t to ) {
    arrayMove( mat->colorder, from, to );
}

void
matpbm_free( bmat_t* mat ) {
    free( mat->rowptr );
    free( mat->rowlen );
    free( mat->colind );
    free( mat->roworder );
    free( mat->colorder );
    free( mat );
}
