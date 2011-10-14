#include <bios/format.h>
#include <bios/log.h>
#include <bios/hlrmisc.h>
#include <bios/bgrParser.h>

int main( int argc, char** argv) 
{
  Array values;
  Array bedGraph; 
  int i;
  bgrParser_initFromPipe ( "tabix 4240sc_P.bgr.gz.bgz chr21:36000000-38000000" );
  bedGraph = bgrParser_getAllEntries();
  for ( i=0; i<arrayMax( bedGraph ); i++ ) {
    BedGraph *currBGR = arrp( bedGraph, i, BedGraph );
    printf( "%s\t%d\t%d\t%f\n", currBGR->chromosome, currBGR->start, currBGR->end, currBGR->value );
  }
  bgrParser_deInit();
  return 0;
}



/*

int c, skip = -1, meta = -1, list_chrms = 0, force = 0, print_header = 0, bed_reg = 0;
  ti_conf_t conf = ti_conf_gff;
  struct stat stat_tbi,stat_vcf;
  char *fnidx = calloc(strlen(argv[optind]) + 5, 1);
  strcat(strcpy(fnidx, argv[optind]), ".tbi");
 
  // retrieve
  tabix_t *t;
  // Common source of errors: new VCF is used with an old index
  stat(fnidx, &stat_tbi);
  /*stat(argv[optind], &stat_vcf);
  if ( stat_vcf.st_mtime > stat_tbi.st_mtime )
    {
      fprintf(stderr, "[tabix] the index file is older than the vcf file. Please use '-f' to overwrite or reindex.\n");
      free(fnidx);
      return 1;
      }//*/
/*  free(fnidx);	
  if ((t = ti_open(argv[optind], 0)) == 0) 
    die(  "[main] fail to open the data file: %s\n", argv[1] );

  if (strcmp(argv[optind+1], ".") == 0) { // retrieve all
    ti_iter_t iter;
    const char *s;
    int len;
    iter = ti_query(t, 0, 0, 0);
    while ((s = ti_read(t, iter, &len)) != 0) {
      fputs(s, stdout); fputc('\n', stdout);
    }
    ti_iter_destroy(iter);
  } else { // retrieve from specified regions
    int i, len;
    ti_iter_t iter;
    const char *s;
    const ti_conf_t *idxconf;
    
    if (ti_lazy_index_load(t) < 0 && bed_reg == 0) {
      die("[tabix] failed to load the index file: %s.tbi\n", argv[1]);
    }
    idxconf = ti_get_conf(t->idx);
    
    /*if ( print_header ) {
	// If requested, print the header lines here
	iter = ti_query(t, 0, 0, 0);
	while ((s = ti_read(t, iter, &len)) != 0) {
	  if ((int)(*s) != idxconf->meta_char) break;
	  fputs(s, stdout); fputc('\n', stdout);
	}
	ti_iter_destroy(iter);
	} //*/
/*    if (bed_reg) {
      extern int bed_overlap(const void *_h, const char *chr, int beg, int end);
      extern void *bed_read(const char *fn);
      extern void bed_destroy(void *_h);
      
      const ti_conf_t *conf_ = idxconf ? idxconf : &conf; // use the index file if available
      void *bed = bed_read(argv[optind+1]); // load the BED file
      ti_interval_t intv;
      
      if (bed == 0) 
	die( "[main] fail to read the BED file: %s\n", argv[1] );
      
      iter = ti_query(t, 0, 0, 0);
      while ((s = ti_read(t, iter, &len)) != 0) {
	int c;
	ti_get_intv(conf_, len, (char*)s, &intv);
	c = *intv.se; *intv.se = '\0';
	if (bed_overlap(bed, intv.ss, intv.beg, intv.end)) {
	  *intv.se = c;
	  puts(s);
	}
	
	*intv.se = c;
      }
      ti_iter_destroy(iter);
      bed_destroy(bed);
    } else {
      for (i = optind + 1; i < argc; ++i) {
	int tid, beg, end;
	if (ti_parse_region(t->idx, argv[i], &tid, &beg, &end) == 0) {
	  iter = ti_queryi(t, tid, beg, end);
	  while ((s = ti_read(t, iter, &len)) != 0) {
	    fputs(s, stdout); fputc('\n', stdout);
	  }
	  ti_iter_destroy(iter);
	} 
	// else fprintf(stderr, "[main] invalid region: unknown target name or minus interval.\n");
      }
    }
  }
  ti_close(t);
   //*/
