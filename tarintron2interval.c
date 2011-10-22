#include <bios/format.h>
#include <bios/log.h>
#include <bios/intervalFind.h>
#include <bios/common.h>

int main( int argc, char** argv) 
{
  Array intervals;
  Array genes;
  Array overlaps;

  Interval* leftInterval;
  Interval* rightInterval;
  Interval* intron;
  int i,j;

   // reading tar interval file
  intervals = intervalFind_parseFile( argv[1], 0 );
  genes = intervalFind_parseFile( argv[2], 0 );

  for( i=0; i<arrayMax( intervals)-1; i++) {
    leftInterval = arrp( intervals, i, Interval );
    rightInterval = arrp( intervals, i+1, Interval );
    if( !strEqual( leftInterval->chromosome, rightInterval->chromosome ) || 
	(rightInterval->start - leftInterval->end ) > 100000 )
      continue;
    AllocVar( intron );
    intron->chromosome = hlr_strdup( leftInterval->chromosome );
    intron->strand = leftInterval->strand;
    intron->start = leftInterval->end+1;
    intron->end = rightInterval->start-1;
    intron->name = hlr_strdup( leftInterval->name );
    intron->subIntervalCount = 1;
    intron->subIntervals = arrayCreate(intron->subIntervalCount , SubInterval);
    SubInterval* subInterval = arrayp( intron->subIntervals, arrayMax( intron->subIntervals ), SubInterval );
    subInterval->start = leftInterval->end+1;
    subInterval->end = rightInterval->start-1;
    
    overlaps = intervalFind_getOverlappingIntervals( intron->chromosome, intron->start, intron->end ); 
    if( arrayMax( overlaps ) == 0 ) {
      printf("%s\n",intervalFind_writeInterval( intron ) );
    } else {
      unsigned short int contained=0;
      for( j = 0; j< arrayMax( overlaps ); j++) {
	Interval* currOverlap = arrp( overlaps, j, Interval );
	if( (intron->start <= currOverlap->start) && ( intron->end >= currOverlap->end ) ) contained=1; // checking if fully contained
      }
      if( !contained ) {
	printf("%s\n",intervalFind_writeInterval( intron ) );
      }
    }
    freeMem( intron );
  }
  return 0;
}
