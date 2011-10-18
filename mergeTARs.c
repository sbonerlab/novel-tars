#include <bios/format.h>
#include <bios/log.h>
#include <bios/hlrmisc.h>
#include <bios/bgrParser.h>
#include <bios/intervalFind.h>
#include <bios/common.h>
#include <math.h>
#include <bios/linestream.h>

typedef struct {
  char* name;
  Array values;
} Expression;

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

int main( int argc, char** argv) 
{
  Array intervals;
  Interval* leftInterval;
  Interval* rightInterval;
  Statistics* leftStats;
  Statistics* intronStats;
  Statistics* rightStats;
  int i, j;

  File
  bgrFileName = hlr_strdup( argv[2] );
  // reading interval file
  intervals = intervalFind_parseFile( argv[1], 0 );
  
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
      Interval* newInterval;
      AllocVar( newInterval );
      newInterval->chromosome= hlr_strdup( leftInterval->chromosome );
      newInterval->strand = leftInterval->strand;
      newInterval->start = leftInterval->start;
      newInterval->end = rightInterval->end;
      newInterval->name = hlr_strdup(  leftInterval->name );
      newInterval->subIntervals = arrayCreate( 1, SubInterval);
      for( j=0; j<1; j++) {
	SubInterval* subInterval = arrayp( newInterval->subIntervals, arrayMax( newInterval->subIntervals ), SubInterval );
	subInterval->start = leftInterval->start;
	subInterval->end = rightInterval->end;
      }
      printf("new:\t%s\n",  intervalFind_writeInterval( newInterval ) );
      freeMem( newInterval );
    }
    printf( "(%d):left\t%1.5e\t%1.5e\tright:\t%1.6e\t%1.6e\t\tIntron\t%e\t%e\n", i, leftStats->mean, leftStats->stdev, rightStats->mean, rightStats->stdev, intronStats->mean, intronStats->stdev);
    freeMem(intron);
  }

  return 0;
}

  /*  

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


