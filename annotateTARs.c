#include <bios/format.h>
#include <bios/log.h>
#include <bios/intervalFind.h>
#include <bios/common.h>
#include <bios/array.h>

#define UTR5 0
#define UTR3 1
#define INTRON 2

int computeDistance( SubInterval *currSubInterval, Interval* currInterval ) {
  int currDistance = currSubInterval->start - currInterval->start ;
  if( currDistance < 0 ) 
    currDistance = abs( currSubInterval->end - currInterval->start );
  else 
    currDistance = abs( currSubInterval->start - currInterval->end );
  return currDistance;
}
int processOverlap(Interval *currOverlap, Interval* currInterval ) {
  int i, index;
  int minDistance, currDistance;
  for( i=0; i < arrayMax( currOverlap->subIntervals ); i++ ) {
    SubInterval *currSubInterval = arrp( currOverlap->subIntervals, i, SubInterval);
    currDistance = computeDistance( currSubInterval, currInterval);
    if( i==0 ) 
      minDistance =  currDistance ;
    else {
      if( currDistance < minDistance ) {
	minDistance = currDistance;
	index = i;
      }
    }
  }
  if( (index == 0 && currOverlap->strand =='+')  || (index == (arrayMax( currOverlap->subIntervals )-1)  && currOverlap->strand == '-') )
    return UTR5;
  else if ( (index == 0 && currOverlap->strand =='-')  || (index == (arrayMax( currOverlap->subIntervals )-1)  && currOverlap->strand == '+') )
    return UTR3;
  else
    return INTRON;
}

int main( int argc, char** argv) 
{
  Array TARs = arrayCreate( 100, Interval ); 
  Array overlaps;
  Interval *currInterval, *currOverlap;
  int i,j;
  int numIntergenic = 0;
  int numIntronic   = 0;
  int num3UTR       = 0;
  int num5UTR       = 0;
  int numMulti      = 0;
  int maxUTRdist    = 2000;
  int type;
  if( argc < 3 ) {
    usage("%s <transcriptome.interval> <TARs.interval> <maxUTRdist>", argv[0]);
  }
  // "/net/gerstein/annotations/human/hg19/GENCODE/v3c/gencode.v3c.annotation.GRCh37.transcript.gtpc.ttpc.composite.interval"
  intervalFind_addIntervalsToSearchSpace ( argv[1], 0);
  // "/net/gerstein/brainSeq/human_adult_project/withoutChrM/NOVEL/all.novel.noChrM.novel.0.0505728.50.100.filtered.interval"
  TARs = intervalFind_parseFile ( argv[2], 1);
  
  maxUTRdist = atoi ( argv[3] );
  if( maxUTRdist < 1 ) 
    die( "Max UTR distance needs to be positve (or non-numeric data): %d", maxUTRdist);
  
  AllocVar(currInterval);
  for( i=0; i < arrayMax( TARs ); i++ ) {
    currInterval = arrp( TARs, i, Interval );
    overlaps =  arrayCopy( intervalFind_getOverlappingIntervals( currInterval->chromosome, currInterval->start - maxUTRdist, currInterval->end + maxUTRdist) );
    if( arrayMax( overlaps ) == 1 ) {
      for( j=0; j<arrayMax( overlaps ); j++ ) {
	currOverlap = arru( overlaps, j, Interval* );
	printf( "%s\tTRANSCRIPT\n", intervalFind_writeInterval( currOverlap ) ); 
	type = processOverlap( currOverlap, currInterval);
      }
      if( type == UTR5) {
	num5UTR++;
	printf( "%s\t5UTR\n\n", intervalFind_writeInterval ( currInterval ));
      }
      else if (type == UTR3 ) {
	num3UTR++;
	printf( "%s\t3UTR\n\n", intervalFind_writeInterval ( currInterval ));
      }
      else {
	numIntronic++;
	printf( "%s\tINTRON\n\n", intervalFind_writeInterval ( currInterval ));
      }
    } else if ( arrayMax( overlaps ) > 1 ) {
      printf( "%s\tmulti\n\n", intervalFind_writeInterval ( currInterval ) );
      numMulti++; //type = processMultipleOverlaps( overlaps, currInterval );
    } else {
      numIntergenic++;
      printf( "%s\tintergenic\n\n", intervalFind_writeInterval ( currInterval ) );
    }
  }
  warn("%s: max distance from UTRs: %d", argv[0], maxUTRdist );
  warn("%s: num intergenic TARs %d", argv[0], numIntergenic );
  warn("%s: num 5'UTR      TARs %d", argv[0], num5UTR );
  warn("%s: num intronic   TARs %d", argv[0], numIntronic );
  warn("%s: num 3'UTR      TARs %d", argv[0], num3UTR );
  warn("%s: num multi      TARs %d", argv[0], numMulti );
  arrayDestroy( TARs );
  arrayDestroy( overlaps );

  return 0;
}



