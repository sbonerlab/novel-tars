#include <bios/format.h>
#include <bios/log.h>
#include <bios/hlrmisc.h>
#include <bios/bgrParser.h>
#include <bios/intervalFind.h>
#include <bios/common.h>
#include <math.h>
#include <bios/linestream.h>
#include <gsl/gsl_statistics.h>

typedef struct {
  char* name;
  char* coordinates;
  Array values; // of double
} Expression;


void readExpression( Array exprs, char* cmd ) {
  LineStream ls;
  char* line;

  ls = ls_createFromPipe( cmd );
  while( line = ls_nextLine( ls ) ) {
    if( strStartsWith( line, "ID") )
	continue;
    WordIter w = wordIterCreate( line, "\t", 0);
    Expression* currEl = arrayp( exprs, arrayMax( exprs ), Expression );
    currEl->name = hlr_strdup( wordNext( w ) );
    currEl->coordinates = hlr_strdup ( wordNext( w ) );
    currEl->values = arrayCreate( 60, double );
    char* strValue;
    while ( strValue = wordNext( w ) ) {
        array( currEl->values, arrayMax( currEl->values ), double) = atof( strValue );
    }
    wordIterDestroy( w );
  }
  ls_destroy( ls );
}

int sortExpressionByName( Expression* exprA, Expression* exprB) {
  return strcmp( exprB->name, exprA->name );
}


int correlation( Expression* exprA, Expression* exprB ) {
  int i;
  int n=arrayMax(exprA->values);
  if( arrayMax( exprB->values) != n ) die("Not the same size");
  double valueA[n] ;
  double valueB[n];
  for( i=0; i<n; i++ ) {
    valueA[i] = *arrp( exprA->values, i, double );
    valueB[i] = *arrp( exprB->values, i, double );
  }
  double corr = gsl_stats_correlation( valueA, 1, valueB, 1, n );
  double meanA =  gsl_stats_mean( valueA, 1, n);
  double meanB = gsl_stats_mean( valueB, 1, n );
  double meandiff = fabs( meanA - meanB );
  double sdA = gsl_stats_sd_m( valueA, 1, n, meanA );
  double sdB = gsl_stats_sd_m( valueB, 1, n, meanB );
  double sddiff = sqrt( sdA*sdA + sdB*sdB );
  
  if( corr > 0 && (sddiff > meandiff) )
    return 1;
  return 0;
}

void writeMergedInterval(Interval* leftInterval, Interval* rightInterval ) {
  if( strEqual( leftInterval->name, rightInterval->name ) ) {
    printf( "%s\n", intervalFind_writeInterval( leftInterval ) );
    return;
  }
  Interval* newInterval;
  AllocVar( newInterval );
  newInterval->chromosome= hlr_strdup( leftInterval->chromosome );
  newInterval->strand = leftInterval->strand;
  newInterval->start = leftInterval->start;
  newInterval->end = rightInterval->end;
  newInterval->name = hlr_strdup(  leftInterval->name );
  newInterval->subIntervalCount = 1;
  newInterval->subIntervals = arrayCreate( newInterval->subIntervalCount, SubInterval);
  SubInterval* subInterval = arrayp( newInterval->subIntervals, newInterval->subIntervalCount - 1, SubInterval );
  subInterval->start = leftInterval->start;
  subInterval->end = rightInterval->end;
  printf("%s\n",  intervalFind_writeInterval( newInterval ) );
  freeMem( newInterval );
}


int main( int argc, char** argv) 
{
  Array intervals; // of Interval
  Array merge; // of short int
  int i, j;
  short int found;
  Stringa buffer = stringCreate( 20 );

  if( argc != 4 ) {
    usage( "%s <tar.interval> <tar.expression.gz> <tarIntronic.expression.gz>", argv[0] );
  }

  Array tars = arrayCreate( 5000, Expression ); // of Expression
  Array introns = arrayCreate( 5000, Expression ); // of Expression

 // reading interval file
  intervals = intervalFind_parseFile( argv[1], 0 );
  merge = arrayCreate( arrayMax(intervals)-1, short int );

  stringPrintf( buffer ,  "zcat %s", argv[2] ); //all.novel.filtered.merged.expression_2011.02.23.csv.gz" );
  readExpression( tars, string(buffer) );
  stringPrintf( buffer,  "zcat %s", argv[3]); //all.tarintronic.expression_2011.10.19.csv.gz" );
  readExpression( introns, string(buffer) );
  stringDestroy( buffer );

  arraySort( tars, (ARRAYORDERF) sortExpressionByName ); 
  arraySort( introns,  (ARRAYORDERF) sortExpressionByName ); 

  for ( i = 0; i<arrayMax( intervals)-1; i++ ) {
    if( strEqual( arrp(intervals, i, Interval)->chromosome, arrp(intervals, i+1, Interval)->chromosome) &&
	(arrp(intervals, i+1, Interval)->start - arrp(intervals, i, Interval)->end) < 100000 ) { // same chr and < 100Kb
      Expression* testFind;
      AllocVar( testFind );
      testFind->name = hlr_strdup( arrp(intervals, i, Interval)->name );
      found = arrayFind( tars, testFind, &j, (ARRAYORDERF) sortExpressionByName );
      if( !found )
	die( "TAR not found %s", testFind->name);
      Expression* exprTARL = arrp( tars, j, Expression );
      
      testFind->name = hlr_strdup( arrp(intervals, i+1, Interval)->name );
      found = arrayFind( tars, testFind, &j, (ARRAYORDERF) sortExpressionByName );
      if( !found )
	die( "TAR not found %s", testFind->name);
      Expression* exprTARR = arrp( tars, j, Expression );
      freeMem (testFind );
      
      if ( correlation( exprTARL, exprTARR) ) {
	found = arrayFind( introns, exprTARL, &j, (ARRAYORDERF) sortExpressionByName );
	if( !found ) {
	  arru( merge, i, short int) = 0;
	  warn("TAR not found %s", exprTARL->name );
	  continue;
	}
	Expression* currIntron;
	AllocVar( currIntron );
	currIntron = arrp( introns, j, Expression );
	if( correlation( exprTARL, currIntron ) && correlation( currIntron, exprTARR) ) 
	  arru( merge, i, short int) = 1;
	else 
	  arru( merge, i, short int) = 0;
      }
    } else
      arru( merge, i, short int) = 0;
  }
  i = 0;
  while( i < (arrayMax( intervals )-1) ) {
    Interval* leftInterval = arrp( intervals, i, Interval );
    while( arru( merge, i, short int) && (i < (arrayMax( intervals )-1)) ) i++;
    Interval* rightInterval = arrp( intervals, i, Interval );
    writeMergedInterval( leftInterval, rightInterval ); // it takes care of non-merged intervals
    i++;
  }
	
  arrayDestroy( introns );
  arrayDestroy( tars );
  arrayDestroy( intervals );
  arrayDestroy( merge );
  return 0;
}

 
