/**
  * Rareman - simple C implementation of Hellerman-Rarick transform
  *
  * usage: rareman <input pbm file> [output pbm file]
  * Outfile is optional and ascii art is generated if none is specified.
  *
  * Copyright 2018 Micky Faas <micky@edukitty.org>
  */

#include <stdio.h>
#define _GNU_SOURCE // for qsort_r()
#include <stdlib.h>
#include <string.h>

#include "matpbm.h"

#define IDX_LAST (1 << 31)

int
compareByMagnitude( const void* a, const void* b, void* arg ) {
    idx_t elemA = *(idx_t*)a;
    idx_t elemB = *(idx_t*)b;
    idx_t *magnitude = (idx_t*)arg;

    return magnitude[elemA] - magnitude[elemB];
}

int
isIn( idx_t value, idx_t* set, size_t n ) {
    for( int i =0; i < n; i++ ) {
        if( set[i] == IDX_LAST ) break;
        if( set[i] == value ) return i;
    }

    return -1;
}

void
hr_transform( bmat_t* mat ) {

    // Start of the active submatrix
    mat->active.p = mat->active.q =0;
    // Array of active rows and columns
    idx_t a_rows[mat->m];
    idx_t a_cols[mat->m];
    idx_t n_intersect[mat->m];
    idx_t n_nz[mat->m];
    idx_t pivot_cols[mat->m];
    size_t numPivotCols =0;
    int m =0;

    while( mat->active.p + mat->active.q < mat->m ) {

        int progress = mat->active.p + mat->active.q;
        fprintf( stderr, "\b\b\b\b\b%4g%%", (double)progress / (double)mat->m * 100.0 );

        numPivotCols =0;
        m =0;
        while( 1 ) {
            
            // Active columns that are not pivot columns
            int numCols =0;
            for( int j =mat->active.p; j < mat->m - mat->active.q; j++ ) {
                if( isIn( j, pivot_cols, numPivotCols ) != -1 ) continue;
                a_cols[numCols++] = j;
            }
            
            // Number of rows in the active matrix
            int numRows = mat->m - mat->active.p;

            // Initialize an array of all rows in the active matrix
            for( int i =0; i < numRows; i++ ) {
                idx_t row =mat->active.p + i;
                a_rows[i] = row;
                n_nz[row] = 0;

               // matpbm_rowNNZ( mat, row, mat->active.q, mat->m - mat->active.q );
                for( int j =0; j < numCols; j++ ) {
                    if( isIn( mat->colorder[a_cols[j]], 
                              &mat->colind[mat->rowptr[mat->roworder[row]]], 
                              mat->rowlen[mat->roworder[row]] ) != -1 )
                        n_nz[row]++;
                }
                if( n_nz[row] == 0 ) n_nz[row] = mat->m;
            }

            // Sort the active rows by non-zero count, ascending
            qsort_r( a_rows, numRows, sizeof( idx_t ), compareByMagnitude, (void*)n_nz );

            // The smallest number of non-zeroes in an active row
            idx_t min = n_nz[a_rows[0]];
            
            // Determine the maximum size of this block at the first iteration of the while()
            if( m == 0 ) 
                m = MIN( min, mat->m - mat->active.p - mat->active.q );

            /*printf( "p=%d, q=%d, m=%d, \nrow: ", mat->active.p, mat->active.q, m );
            for( int i =0; i < numRows; i++ )
                printf( "%d, ", a_rows[i] );
            printf ("\nnnz: " );
            for( int i =0; i < numRows; i++ )
                printf( "%d, ", n_nz[a_rows[i]] );
            printf( "\n" );*/

            // Count the intersections with the active rows containing (min) non-zeroes,
            // and the active non-pivot columns
            for( int j =0; j < numCols; j++ ) {
                int col = a_cols[j];
                // Count the number intersections between this column, 
                // and rows with (min) non-zeroes
                n_intersect[col] =0;
                for( int i =0; i < numRows; i++ ) { 
                    int k =a_rows[i];
                    if( n_nz[k] > min ) break;
                    if( isIn( mat->colorder[col], 
                              &mat->colind[mat->rowptr[mat->roworder[k]]], 
                              mat->rowlen[mat->roworder[k]] ) != -1 )
                        n_intersect[col]++;
                }
            }

            // Sort the active columns by number of intersections. ascending
            qsort_r( a_cols, numCols, sizeof( idx_t ), compareByMagnitude, (void*)n_intersect );
            
           /* printf( "col: " );
            for( int i =0; i < numCols; i++ )
                printf( "%d, ", a_cols[i] );
            printf ("\nits: " );
            for( int i =0; i < numCols; i++ )
                printf( "%d, ", n_intersect[a_cols[i]] );
            printf( "\n" );*/

            // The column with the most intersections, add it to the list of pivot columns
            int col =a_cols[numCols -1];
            pivot_cols[numPivotCols++] = col;


            if( min == 1 || numPivotCols == m ) { // We've completed a block

                int s =0;
                if( min == 1 )
                    s = n_intersect[col]; // Number of singleton rows = size of the block
                if( m < s ) s = m;

                //printf( "s=%d, m=%d, numPivotCols=%d\n", s, m, numPivotCols );

                for( int i =0, k =0; i < numRows; i++ ) {
                    int row = a_rows[i];
                    if( n_nz[row] != min ) break;
                    if( isIn( mat->colorder[col], 
                              &mat->colind[mat->rowptr[mat->roworder[row]]], 
                              mat->rowlen[mat->roworder[row]] ) != -1 ) {
                        a_rows[i] = a_rows[k];
                        a_rows[k++] = row;
                    }
                }

                // Convert logical row numbers to physical row numbers
                for( int i =0; i < s; i++ ) {
                    a_rows[i] = mat->roworder[a_rows[i]];
                }
                // Convert logical column numbers to physical column numbers
                for( int i =0; i < numPivotCols; i++ ) {
                    pivot_cols[i] = mat->colorder[pivot_cols[i]];
                }

                // Permute the s singleton rows to the front
                for( int i =0; i < s; i++ ) {
                    int logical = isIn( a_rows[i], mat->roworder, mat->m );
                    matpbm_moveRow( mat, logical, mat->active.p );
                }

                // Permute the last s pivot columns to the front
                for( int j =0; j < s; j++ ) {
                    int logical = isIn( pivot_cols[numPivotCols-j-1], mat->colorder, mat->m );
                    matpbm_moveCol( mat, logical, mat->active.p );
                }

                int spikes =m - s;

                // Permute the spikes to the border
                for( int j =0; j < spikes; j++ ) {
                    int logical = isIn( pivot_cols[j], mat->colorder, mat->m );
                    matpbm_moveCol( mat, logical, mat->m - mat->active.q - 1 );
                }

                mat->active.p += s;
                mat->active.q += spikes;
                break;
            }   

        }
       
       /* matpbm_printDense( mat );
        char c;
        scanf( "%c", &c );*/

    }
    fprintf( stderr, "\b\b\b\b\bDone.\n" );

}

int
main( int argc, const char** argv ) {
    FILE *infile, *outfile;
    if( argc < 2 ) {
        infile = stdin;
    } else {
        infile =fopen( argv[1], "r" );
        if( infile == NULL ) {
            fprintf( stderr, "(e) Could not open %s for reading.\n", argv[1] );
        }
    }
    bmat_t* mat = matpbm_loadFromStream( infile );
    if( infile != stdin ) fclose( infile );
    if( !mat ) return -1;
    
    if( argc < 3 ) {
        outfile = stdout;
    } else {
        outfile =fopen( argv[2], "w" );
        if( outfile == NULL ) {
            fprintf( stderr, "(e) Could not open %s for writing.\n", argv[2] );
        }
    }

    //matpbm_printDense( mat );

    hr_transform( mat );
    
    if( outfile == stdout )
        matpbm_printDense( mat );
    else
        matpbm_writeToStream( outfile, mat );

    matpbm_free( mat );
    if( outfile != stdin ) fclose( outfile );

    return 0;
}
