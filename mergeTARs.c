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
  Array intervals; Array merge;
  int i, j;
  short int found;

  Array tars = arrayCreate( 5000, Expression );
  Array introns = arrayCreate( 5000, Expression );
 // reading interval file
  intervals = intervalFind_parseFile( argv[1], 0 );
  merge = arrayCreate( arrayMax(intervals)-1, short int );

  char* buffer = hlr_strdup( "zcat all.novel.filtered.merged.expression_2011.02.23.csv.gz" );
  readExpression( tars, buffer );
  buffer = hlr_strdup( "zcat all.tarintronic.expression_2011.10.19.csv.gz" );
  readExpression( introns, buffer );
  hlr_free( buffer );

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
	if( !found )
	  warn("TAR not found %s", exprTARL->name );
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

  /*  
  bgrFileName = hlr_strdup( argv[2] );
 
  
  for( i=0; i<arrayMax( intervals)-1; i++) {
    leftInterval = arrp( intervals, i, Interval );
    rightInterval = arrp( intervals, i+1, Interval );
    if( !strEqual( leftInterval->chromosome, rightInterval->chromosome ) )
      continue;
    leftStats = intervalStatistics( leftInterval );
    rightStats = intervalStatistics ( rightInterval );
    Interval* intron;
    AllocVar( intron );
    intron->chromosome = hlr_strdup( leftInterval->chromosome );
    intron->start = leftInterval->end+1;
    intron->end = rightInterval->start-1;
    intronStats = intervalStatistics ( intron );
    if( (intronStats->mean + (1.95*intronStats->stdev) ) > (leftStats->mean - (1.95*leftStats->stdev)) ) {
      printf("join:");
     
    }
    printf( "(%d):left\t%1.5e\t%1.5e\tright:\t%1.6e\t%1.6e\t\tIntron\t%e\t%e\n", i, leftStats->mean, leftStats->stdev, rightStats->mean, rightStats->stdev, intronStats->mean, intronStats->stdev);
    freeMem(intron);
  }





typedef struct {
  double mean;
  double stdev;
} Statistics;

static char* bgrFileName;

Statistics* intervalStatistics( Interval* currInterval ) {
  Array bedGraph;
  float sumValues = 0.0;
  float sumValues2 = 0.0;
  int i;
  static Stringa buffer = NULL;
  stringCreateClear(buffer, 50 );
  stringPrintf( buffer, "tabix %s %s:%d-%d", bgrFileName, currInterval->chromosome, currInterval->start, currInterval->end );
  bgrParser_initFromPipe ( string( buffer ) ); // "tabix 4240sc_P.bgr.gz.bgz chr21:36000000-38000000" 
  bedGraph = bgrParser_getAllEntries();
  bgrParser_deInit();
  for ( i=0; i<arrayMax( bedGraph ); i++ ) {
    BedGraph *currBGR = arrp( bedGraph, i, BedGraph );
    double bgrValue = (double) currBGR->value;
    sumValues  +=   bgrValue              * (double) ( currBGR->end - currBGR->start );
    sumValues2 += ( bgrValue * bgrValue ) * (double) ( currBGR->end - currBGR->start );
  }
  Statistics *stats;
  AllocVar( stats );
  float N = (float) ( currInterval->end - currInterval->start );
  stats->mean = sumValues / N;
  stats->stdev = sqrt( (N * sumValues2) - ( sumValues * sumValues ) ) / ( N * (N - 1.0) );
  //printf( "sumValues:\t%1.6e\tsumValues2:\t%1.6e\t(N*s2):%1.6e\t(s1*s1): %1.6e\t(N*s2) - (s1*s1): %1.6e\tall: %1.6e\nmean: %1.6e\tsd: %1.6e\n", sumValues, sumValues2, (N* sumValues2), (sumValues * sumValues ),  (N* sumValues2) - (sumValues * sumValues ), ( N * (N - 1.0) ), stats->mean, stats->stdev );
  return stats;
}



  Array bedGraphIntron; 
  Array bedGraphExonLeft;
  Array bedGraphExonRight;

bgrParser_initFromPipe ( "tabix 4240sc_P.bgr.gz.bgz chr21:35000000-35999999 " );
  bedGraphExonLeft = bgrParser_getAllEntries();
  bgrParser_deInit();
  bgrParser_initFromPipe ( "tabix 4240sc_P.bgr.gz.bgz chr21:38000001-39000000" );
  bedGraphExonRight = bgrParser_getAllEntries();
  bgrParser_deInit();
  float valueLeft=0.0, valueIntron=0.0, valueRight=0.0;
  for ( i=0; i<arrayMax( bedGraphExonLeft ); i++ ) {
    BedGraph *currBGR = arrp( bedGraphExonLeft, i, BedGraph );
    valueLeft += currBGR->value * ( currBGR->end - currBGR->start);
  }
 for ( i=0; i<arrayMax( bedGraphExonRight ); i++ ) {
    BedGraph *currBGR = arrp( bedGraphExonRight, i, BedGraph );
    valueRight += currBGR->value * ( currBGR->end - currBGR->start);
  }
  printf( "valueLeft:\t %f\nvalueIntron:\t%f\nvalueRight:\t%f\n", (valueLeft*1000.0)/(35999999-35000000+1), (valueIntron*1000.0)/(38000000-36000000+1), (valueRight*1000.0)/(39000000-38000001+1));
  */
    //   printf( "%s\t%d\t%d\t%f\n", currBGR->chromosome, currBGR->start, currBGR->end, currBGR->value );


