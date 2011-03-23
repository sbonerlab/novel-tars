#include "log.h"
#include "format.h"
#include "blatParser.h"

int checkOriginal ( BlatQuery* blQ, PslEntry* entry ) 
{
  int res = 0;
  Texta  t = textFieldtokP( blQ->qName, "|");
  if( ( strEqual ( textItem( t, 1 ), entry->tName ) ) &&
      ( atoi( textItem( t, 3 ) ) == entry->tStart ) &&
      ( atoi( textItem( t, 4 ) ) == entry->tEnd ) && 
      entry->misMatches==0 && (entry->matches+entry->nCount) == entry->qSize) 
    res = 1;
  else
    res = 0;
  textDestroy( t );
  return res;
}

int processBlatQuery( BlatQuery* blQ, int *idxOrig , float cutoff) {
  int i,j;
  PslEntry *curr;
  int sizes[ arrayMax( blQ->entries ) ];
  *idxOrig = -1;
  for( i=0; i<arrayMax(blQ->entries); i++ ) {
    curr = arrp( blQ->entries, i, PslEntry );
    sizes[i]=0;
    for( j=0; j < arrayMax( curr->blockSizes ); j++) 
      sizes[i] += arru( curr->blockSizes, j, int);
    sizes[i] -=  curr->misMatches;
    if( checkOriginal ( blQ, curr ) == 1 )  
      *idxOrig = i;
  }
  if( *idxOrig < 0 ) die("Cannot find exact match: %s", blQ->qName);
  int sizeOrig = sizes[ *idxOrig ];
  for( i=0; i< arrayMax( blQ->entries ); i++ ) {
    curr = arrp( blQ->entries, i, PslEntry );
    warn( "%s\t%s\t%d\t%d\t[ %d, %d - %f]\t%d\t--\t%d\t%d\t%d", blQ->qName, curr->tName, curr->tStart, curr->tEnd, sizes[i], sizeOrig, ( (float)sizes[i] / (float)sizeOrig ) ,curr->blockCount, curr->misMatches, curr->qNumInsert, curr->tNumInsert); 
    if( ( (float)sizes[i] / (float)sizeOrig ) > cutoff && (i != *idxOrig) )
      return 1;
  }
  return 0;
}

int main (int argc, char *argv[])
{
  if( argc < 2 ) {
    usage("%s <overlap> < file.psl", argv[0]);
  }
  blatParser_initFromFile( "-" );
  BlatQuery* blQ = NULL;
  PslEntry* pslE = NULL;
  int discard = 0;
  int idxOrig ;
  while( blQ = blatParser_nextQuery() ) {
    idxOrig = 0;
    if( arrayMax( blQ->entries ) > 1 ) {
      idxOrig = -1;
      discard = processBlatQuery( blQ, &idxOrig, atof( argv[1] ) );
    }
    if( discard == 1 ) {
      discard = 0;
      continue;
    } else {
      if( idxOrig == -1 )
	die( "Error");
      pslE = arrp( blQ->entries, idxOrig, PslEntry );
      printf( "%s\t%d\t%d\t%s\n", pslE->tName, pslE->tStart, pslE->tEnd, blQ->qName);
    }
  }
  blatParser_deInit();
  return 0;
}
