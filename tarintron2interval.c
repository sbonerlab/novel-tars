#include <bios/format.h>
#include <bios/log.h>
#include <bios/intervalFind.h>
#include <bios/common.h>

int main( int argc, char** argv) 
{
  Array intervals;
  Interval* leftInterval;
  Interval* rightInterval;
  Interval* intron;
  int i,j;

   // reading interval file
  intervals = intervalFind_parseFile( argv[1], 0 );

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
    intron->subIntervals = arrayCreate( 1, SubInterval);
    intron->subIntervalCount = 1;
    for( j=0; j<1; j++) {
      SubInterval* subInterval = arrayp( intron->subIntervals, arrayMax( intron->subIntervals ), SubInterval );
      subInterval->start = leftInterval->end+1;
      subInterval->end = rightInterval->start-1;
    }
    printf("%s\n",intervalFind_writeInterval( intron ) );
    freeMem( intron );
  }
  return 0;
}
